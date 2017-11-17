#ifndef STREAMCONTAINER_H_
#define STREAMCONTAINER_H_

#include "../third_party/zlib/zconf.h"
#include "../third_party/zlib/zlib.h"
#include "../io/BasicBuffer.h"
#include "../algorithm/OpenHashTable.h"

namespace Tomahawk{
namespace Core{

enum CORE_TYPE{TYPE_8B, TYPE_16B, TYPE_32B, TYPE_64B,
			   TYPE_FLOAT, TYPE_DOUBLE};

enum CORE_COMPRESSION{ENCODE_NONE, ENCODE_DEFLATE, ENCODE_FSE};

// Controller type for stream container
struct StreamContainerHeaderController{
	typedef StreamContainerHeaderController self_type;

public:
	StreamContainerHeaderController() :
		signedness(0),
		mixedStride(0),
		type(0),
		encoder(0),
		uniform(0),
		unused(0)
	{}
	~StreamContainerHeaderController(){}

	inline void clear(){ memset(this, 0, sizeof(U16)); }

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& controller){
		const U16* c = reinterpret_cast<const U16* const>(&controller);
		stream.write(reinterpret_cast<const char*>(&c), sizeof(U16));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& controller){
		U16* c = reinterpret_cast<U16*>(&controller);
		stream.read(reinterpret_cast<char*>(c), sizeof(U16));
		return(stream);
	}

public:
	// 6 base values (4 integers + 2 floats)
	U16 signedness: 1,
		mixedStride: 1,
		type: 6,    // base typing (extra bits saved for future use)
		encoder: 5, // encoder bits (0 = uncompressed)
		uniform: 1, // triggered if all values are the same
		unused: 2;
};

/*
 Stream header data structure
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+----+
 |   2   |   2   |       4       |        4      |       4       |       4       | X ~
 +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
 ^       ^       ^               ^               ^               ^               ^
 |       |       |               |               |               |               |
 CNT     STRIDE  |               COMP LENGTH     UNCOMP LENGTH   CRC             POSSIBLE PARAMETERS
                 OFFSET
*/
struct StreamContainerHeader{
	typedef StreamContainerHeader self_type;
	typedef StreamContainerHeaderController controller_type;

	StreamContainerHeader() :
		stride(-1),
		offset(0),
		cLength(0),
		uLength(0),
		crc(0),
		n_extra(0),
		extra(nullptr)
	{}

	virtual ~StreamContainerHeader(){
		delete [] this->extra;
	}

	inline void reset(void){
		this->stride = -1;
		this->offset = 0;
		this->cLength = 0;
		this->uLength = 0;
		this->crc = 0;
		this->n_extra = 0;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.controller;
		stream.write(reinterpret_cast<const char*>(&entry.stride), sizeof(S16));
		stream.write(reinterpret_cast<const char*>(&entry.offset), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.cLength),sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uLength),sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.crc),    sizeof(U32));
		if(entry.n_extra > 0)
			stream.write(entry.extra, entry.n_extra);

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream >> entry.controller;
		stream.read(reinterpret_cast<char*>(&entry.stride),  sizeof(S16));
		stream.read(reinterpret_cast<char*>(&entry.offset),  sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.cLength), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uLength), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.crc),     sizeof(U32));
		return(stream);
	}

public:
	controller_type controller;
	S16 stride;
	U32 offset;
	U32 cLength;
	U32 uLength;
	U32 crc;
	U32 n_extra; // not written; used internally only
	char* extra; // extra length is encoder specific
};

/*
 If stride is mixed then use this
 secondary structure. Stride is always
 a single unsigned integer
 +---+---+---+---+---+---+---+---+---+
 |   2   |   2   |       4       | X ~
 +---+---+---+---+---+---+---+---+---+
 ^       ^       ^               ^
 |       |       |               |
 CNT     CLENGTH |               POSSIBLE PARAMETERS
                 ULENGTH
*/
struct StreamContainerHeaderStride{
	typedef StreamContainerHeaderStride self_type;
	typedef StreamContainerHeaderController controller_type;

	StreamContainerHeaderStride() :
		cLength(0),
		uLength(0),
		crc(0),
		n_extra(0),
		extra(nullptr)
	{}
	~StreamContainerHeaderStride(){
		delete [] this->extra;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.controller;
		stream.write(reinterpret_cast<const char*>(&entry.cLength), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.uLength), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.crc),     sizeof(U32));
		if(entry.n_extra > 0)
			stream.write(entry.extra, entry.n_extra);

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream >> entry.controller;
		stream.read(reinterpret_cast<char*>(&entry.cLength), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.uLength), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.crc),     sizeof(U32));
		return(stream);
	}

public:
	controller_type controller;
	U32 cLength;
	U32 uLength;
	U32 crc;
	U32 n_extra; // not written; used internally only
	char* extra; // extra length is encoder specific
};

// Stream container for importing
class StreamContainer{
	typedef StreamContainer self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef StreamContainerHeader header_type;
	typedef StreamContainerHeaderStride header_stride_type;

public:
	StreamContainer() :
		n_entries(0),
		n_additions(0)
	{}

	StreamContainer(const U32 start_size) :
		n_entries(0),
		n_additions(0),
		buffer_data(start_size),
		buffer_strides(start_size)
	{}

	~StreamContainer(){ this->buffer_data.deleteAll(); this->buffer_strides.deleteAll(); }

	inline void setType(const CORE_TYPE value){ this->header.controller.type = value; }
	inline void setStrideSize(const S32 value){ this->header.stride = value; }
	inline const bool checkStrideSize(const S32& value) const{ return this->header.stride == value; }
	inline void setMixedStrides(void){ this->header.stride = -1; }

	inline void operator++(void){ ++this->n_additions; }

	inline void addStride(const U32& value){ this->buffer_strides += value; }

	inline void operator+=(const SBYTE& value){
		//assert(this->header.controller.type == 0);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const BYTE& value){
		//assert(this->header.controller.type == 1);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const S16& value){
		//assert(this->header.controller.type == 2);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const U16& value){
		//assert(this->header.controller.type == 3);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const S32& value){
		//assert(this->header.controller.type == 4);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const U32& value){
		//assert(this->header.controller.type == 5);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const U64& value){
		//assert(this->header.controller.type == 6);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const float& value){
		//assert(this->header.controller.type == 7);
		this->buffer_data += value;
		++this->n_entries;
	}

	inline void operator+=(const double& value){
		//assert(this->header.controller.type == 8);
		this->buffer_data += value;
		++this->n_entries;
	}

	void reset(void){
		this->n_entries = 0;
		this->n_additions = 0;
		this->buffer_data.reset();
		this->buffer_strides.reset();
		this->header.reset();
	}

	inline void resize(const U32 size){
		this->buffer_data.resize(size);
		this->buffer_strides.resize(size);
	}

	bool checkSum(bool both = false){
		if(this->buffer_data.size() == 0)
			return false;

		// Checksum for main buffer
		U32 crc = crc32(0, NULL, 0);
		crc = crc32(crc, (Bytef*)this->buffer_data.data, this->buffer_data.pointer);
		this->header.crc = crc;

		// Checksum for strides
		if(both){
			if(this->buffer_strides.size() > 0){
				crc = crc32(0, NULL, 0);
				crc = crc32(crc, (Bytef*)this->buffer_strides.data, this->buffer_strides.pointer);
				this->header_stride.crc = crc;
			}
		}
		return true;
	}

	bool checkUniformity(void){
		if(this->n_entries == 0)
			return false;

		const U32 stride_size = this->header.stride;
		if(stride_size == -1)
			return false;

		U32 stride_update = stride_size;

		switch(this->header.controller.type){
		case CORE_TYPE::TYPE_32B: stride_update *= sizeof(S32);   break;
		case CORE_TYPE::TYPE_FLOAT: stride_update *= sizeof(float); break;
		case CORE_TYPE::TYPE_8B: stride_update *= sizeof(char);  break;
		default: return false; break;
		}

		const U64 first_hash = XXH64(&this->buffer_data.data[0], stride_update, 2147483647);

		for(U32 i = 1; i < this->n_entries; ++i){
			if(XXH64(&this->buffer_data.data[i*stride_update], stride_update, 2147483647) != first_hash){
				std::cerr << "not uniform" << std::endl;
				return(false);
			}
		}
		std::cerr << "is uniform" << std::endl;

		this->n_entries = 1;
		this->buffer_data.pointer = stride_size;
		this->header.controller.uniform = true;
		this->header.controller.mixedStride = false;
		this->header.controller.encoder = 0;

		return(true);
	}

	/*
	 This function is called during import to shrink each
	 word-type to fit min(x) and max(x) in the worst case.
	 At this stage all integer values in the stream is of
	 type S32. No other values can be shrunk

	 */
	void reformatStream(void){
		// Recode integer types
		if(this->header.controller.type != TYPE_32B || this->header.controller.signedness != 1){
			return;
		}

		const S32* dat = reinterpret_cast<const S32*>(this->buffer_data.data);
		S32 min = dat[0];
		S32 max = dat[0];

		for(U32 j = 1; j < this->n_entries; ++j){
			if(dat[j] < min) min = dat[j];
			if(dat[j] > max) max = dat[j];
		}

		BYTE byte_width = 0;
		if(min < 0) byte_width = ceil((ceil(log2(abs(min) + 1))+1)/8);  // One bit is used for sign
		else byte_width = ceil(ceil(log2(max + 1))/8);

		if(byte_width >= 3 && byte_width <= 4) byte_width = 4;
		else if(byte_width > 4) byte_width = 8;
		if(byte_width == 0) byte_width = 1;

		// Phase 2
		// Here we re-encode values using the smallest possible
		// word-size
		if(this->header.controller.uniform){
			// Non-negative
			this->buffer_data.pointer = 0;
			this->n_entries = 1;
			this->n_additions = 1;
			if(min >= 0){
				this->header.controller.signedness = 0;
				switch(byte_width){
				case 1: this->buffer_data += (BYTE)min; this->header.controller.type = TYPE_8B;  break;
				case 2: this->buffer_data += (U16)min;  this->header.controller.type = TYPE_16B; break;
				case 4: this->buffer_data += (U32)min;  this->header.controller.type = TYPE_32B; break;
				case 8: this->buffer_data += (U64)min;  this->header.controller.type = TYPE_64B; break;
				default: std::cerr << "illegal: " << std::endl; exit(1);
				}
			} else {
				this->header.controller.signedness = 1;
				switch(byte_width){
				case 1: this->buffer_data += (SBYTE)min; this->header.controller.type = TYPE_8B;  break;
				case 2: this->buffer_data += (S16)min;   this->header.controller.type = TYPE_16B; break;
				case 4: this->buffer_data += (S32)min;   this->header.controller.type = TYPE_32B; break;
				default: std::cerr << "illegal" << std::endl; exit(1);
				}
			}
			// done
			std::cerr << "swap: uniform" << std::endl;
			return;
		}
		// Not unfirom
		else {
			buffer_type buffer(this->buffer_data.pointer);
			dat = reinterpret_cast<const S32*>(this->buffer_data.data);

			// Is non-negative
			if(min >= 0){
				this->header.controller.signedness = 0;

				if(byte_width == 1){
					this->header.controller.type = TYPE_8B;

					for(U32 j = 0; j < this->n_entries; ++j)
						buffer += (BYTE)*(dat++);
				} else if(byte_width == 2){
					this->header.controller.type = TYPE_16B;

					for(U32 j = 0; j < this->n_entries; ++j)
						buffer += (U16)*(dat++);
				} else if(byte_width == 4){
					this->header.controller.type = TYPE_32B;

					for(U32 j = 0; j < this->n_entries; ++j)
						buffer += (U32)*(dat++);
				} else if(byte_width == 8){
					this->header.controller.type = TYPE_64B;

					for(U32 j = 0; j < this->n_entries; ++j)
						buffer += (U64)*(dat++);
				} else {
					std::cerr << "illegal" << std::endl;
					exit(1);
				}
			}
			// Is negative
			else {
				this->header.controller.signedness = 1;

				if(byte_width == 1){
					this->header.controller.type = TYPE_8B;

					for(U32 j = 0; j < this->n_entries; ++j)
						buffer += (SBYTE)*(dat++);
				} else if(byte_width == 2){
					this->header.controller.type = TYPE_16B;

					for(U32 j = 0; j < this->n_entries; ++j)
						buffer += (S16)*(dat++);
				} else if(byte_width == 4){
					this->header.controller.type = TYPE_32B;

					for(U32 j = 0; j < this->n_entries; ++j)
						buffer += (S32)*(dat++);
				} else {
					std::cerr << "illegal" << std::endl;
					exit(1);
				}
			}
			std::cerr << "swap: " << this->buffer_data.pointer << '\t' << buffer.pointer << std::endl;
			memcpy(this->buffer_data.data, buffer.data, buffer.pointer);
			this->buffer_data.pointer = buffer.pointer;
			buffer.deleteAll();
		}
	}

	/////////////////////
	// Loading a container
	/////////////////////
	// Sketch of algorithm
	// 1) Load header
	// 2) Read compressed data into temp buffer
	// 3) Decompress into local buffer
	// 4) Finished if stride is not equal to -1
	// 5) Load stride header
	// 6) Read stride data into temp buffer
	// 7) Decompress into local buffer
	bool read(std::istream& stream, buffer_type& temp_buffer){
		stream >> this->header;
		// Todo: currently no encoding-specific fields
		stream.read(temp_buffer.data, this->header.cLength);
		if(!stream.good()){
			std::cerr << Helpers::timestamp("ERROR","CONTAINER") << "Stream is not good!" << std::endl;
			return false;
		}
		// Uncompress into local buffer

		if(this->header.controller.mixedStride){
			stream >> this->header_stride;
			stream.read(temp_buffer.data, this->header_stride.cLength);
			if(!stream.good()){
				std::cerr << Helpers::timestamp("ERROR","CONTAINER") << "Stream is not good!" << std::endl;
				return false;
			}
			// Uncompress into local buffer
		} // end stride extra

		// check crc checksum

		return true;
	}

public:
	// Headers
	header_type header;
	header_stride_type header_stride;
	// Not written - used internally only
	U32 n_entries;   // number of entries
	U32 n_additions; // number of times added
	// Buffers
	buffer_type buffer_data;
	buffer_type buffer_strides;
};

}
}

#endif /* STREAMCONTAINER_H_ */
