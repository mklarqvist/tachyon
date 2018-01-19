#ifndef CORE_STREAMCONTAINER_H_
#define CORE_STREAMCONTAINER_H_

#include <bitset>

#include "../io/BasicBuffer.h"
#include "ContainerHeaderController.h"
#include "ContainerHeader.h"

namespace Tachyon{
namespace Core{

enum TACHYON_CORE_TYPE{
	YON_TYPE_8B,
	YON_TYPE_16B,
	YON_TYPE_32B,
	YON_TYPE_64B,
	YON_TYPE_FLOAT,
	YON_TYPE_DOUBLE,
	YON_TYPE_BOOLEAN,
	YON_TYPE_CHAR,
	YON_TYPE_STRUCT
};

enum TACHYON_CORE_COMPRESSION{
	YON_ENCODE_NONE,
	YON_ENCODE_ZSTD
};


// Stream container for importing
class Container{
	typedef Container self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef ContainerHeader header_type;
	typedef ContainerHeaderStride header_stride_type;

public:
	Container();
	Container(const U32 start_size);
	~Container();

	/**<
	 *
	 * @param value
	 */
	inline void setType(const TACHYON_CORE_TYPE value){ this->header.controller.type = value; }

	/**<
	 *
	 * @param value
	 */
	inline void setStrideSize(const S32 value){ this->header.stride = value; }

	/**<
	 *
	 * @param value
	 * @return
	 */
	inline const bool checkStrideSize(const S32& value) const{ return this->header.stride == value; }

	/**<
	 *
	 */
	inline void setMixedStrides(void){ this->header.stride = -1; this->header.controller.mixedStride = true; }

	// Operators
	inline void operator++(void){ ++this->n_entries; }
	inline void addStride(const U32& value){ this->buffer_strides_uncompressed += (U32)value; }

	inline void operator+=(const SBYTE& value){
		//assert(this->header.controller.type == 0);
		this->buffer_data_uncompressed += value;
		++this->n_additions;
	}

	inline void operator+=(const BYTE& value){
		//assert(this->header.controller.type == 1);
		this->buffer_data_uncompressed += value;
		++this->n_additions;
	}

	inline void operator+=(const S16& value){
		//assert(this->header.controller.type == 2);
		this->buffer_data_uncompressed += value;
		++this->n_additions;
	}

	inline void operator+=(const U16& value){
		//assert(this->header.controller.type == 3);
		this->buffer_data_uncompressed += value;
		++this->n_additions;
	}

	inline void operator+=(const S32& value){
		//assert(this->header.controller.type == 4);
		this->buffer_data_uncompressed += value;
		++this->n_additions;
	}

	inline void operator+=(const U32& value){
		//assert(this->header.controller.type == 5);
		this->buffer_data_uncompressed += value;
		++this->n_additions;
	}

	inline void operator+=(const U64& value){
		//assert(this->header.controller.type == 6);
		this->buffer_data_uncompressed += value;
		++this->n_additions;
	}

	inline void operator+=(const float& value){
		//assert(this->header.controller.type == 7);
		this->buffer_data_uncompressed += value;
		++this->n_additions;
	}

	inline void operator+=(const double& value){
		//assert(this->header.controller.type == 8);
		this->buffer_data_uncompressed += value;
		++this->n_additions;
	}

	void reset(void);
	void resize(const U32 size);

	inline const U64& getSizeUncompressed(void) const{ return(this->buffer_data_uncompressed.size()); }
	inline const U64& getSizeCompressed(void) const{ return(this->buffer_data.size()); }


	/**<
	 * Generates a CRC32 checksum of the uncompressed
	 * data and, if set, the uncompressed strides data.
	 * CRC32 checksums are stored in the header
	 * @return
	 */
	const bool generateCRC(void);

	/**<
	 *
	 * Targets:
	 * 0: Uncompressed data
	 * 1: Uncompressed strides data
	 * 2: Compressed data
	 * 3: Compressed strides data
	 *
	 * @param target Target buffer stream
	 * @return       Returns TRUE if passing checks or FALSE otherwise
	 */
	bool checkCRC(int target = 0);

	/**<
	 *
	 * @return
	 */
	bool checkUniformity(void);

	/**<
	 *
	 */
	void reformat(void);

	/**<
	 *
	 */
	void reformatStride(void);

	inline const U32 getObjectSize(void) const{
		U32 total_size = 0;
		total_size += header.getObjectSize();
		if(this->header.controller.mixedStride)
			total_size += header_stride.getObjectSize();

		total_size += this->buffer_data.pointer;
		if(this->header.controller.mixedStride)
			total_size += this->buffer_strides.pointer;

		return(total_size);
	}

	inline const TACHYON_CORE_TYPE getDataPrimitiveType(void) const{ return(TACHYON_CORE_TYPE(this->header.controller.type)); }
	inline const TACHYON_CORE_TYPE getStridePrimitiveType(void) const{ return(TACHYON_CORE_TYPE(this->header_stride.controller.type)); }

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.header;
		if(entry.header.controller.mixedStride)
			stream << entry.header_stride;

		stream << entry.buffer_data;
		if(entry.header.controller.mixedStride)
			stream << entry.buffer_strides;

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.header;
		if(entry.header.controller.mixedStride)
			stream >> entry.header_stride;

		entry.buffer_data.resize(entry.header.uLength);
		stream.read(entry.buffer_data.data, entry.header.cLength);
		entry.buffer_data.pointer = entry.header.cLength;
		if(entry.header.controller.mixedStride){
			entry.buffer_strides.resize(entry.header_stride.uLength);
			stream.read(entry.buffer_strides.data, entry.header_stride.cLength);
			entry.buffer_strides.pointer = entry.header_stride.cLength;
		}
		return(stream);
	}

public:
	header_type header;
	header_stride_type header_stride;
	// Not written - used internally only during import
	U32 n_entries;   // number of container entries
	U32 n_additions; // number of times an addition operation was executed
	// Buffers - only bit that are written to disk
	// from here
	buffer_type buffer_data;
	buffer_type buffer_strides;
	// These buffers are for internal use only
	// They are used during decompression and are
	// not written to disk
	buffer_type buffer_data_uncompressed;
	buffer_type buffer_strides_uncompressed;
};

}
}

#endif /* CORE_STREAMCONTAINER_H_ */
