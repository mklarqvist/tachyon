#ifndef CORE_STREAMCONTAINER_H_
#define CORE_STREAMCONTAINER_H_

#include <cassert>

#include "../support/enums.h"
#include "../io/basic_buffer.h"
#include "components/datacontainer_header.h"
#include "components/datacontainer_header_controller.h"

namespace tachyon{
namespace containers{

/**<
 * Primary data container in Tachyon. Stores the data itself
 * in compressed/uncompressed form and possibly the data stride
 * size in compressed/uncompressed form
 */
class DataContainer{
	typedef DataContainer       self_type;
	typedef DataContainerHeader header_type;
	typedef io::BasicBuffer     buffer_type;

public:
	DataContainer();
	DataContainer(const U32 start_size);
	~DataContainer();

	/**<
	 * Set the primitive (or higher-order primitive) type for
	 * the objects in this container.
	 * @param value Data primitive type
	 */
	inline void setType(const TACHYON_CORE_TYPE value){ this->header.data_header.controller.type = value; }

	/**<
	 * Set the stride size of this container to some value.
	 * The value -1 is reserved for the case when the controller
	 * flag for mixed strides is set to FALSE. Valid numbers are
	 * [0,..inf)
	 * @param value Stride size
	 */
	inline void setStrideSize(const S32 value){ this->header.data_header.stride = value; }

	/**<
	 * Check if the stride size of this container matches the
	 * supplied stride size. Always returns FALSE if the container
	 * has mixed stride sizes
	 * @param value Stride size to compare against
	 * @return      Returns TRUE if they are the same or FALSE otherwise
	 */
	inline const bool checkStrideSize(const S32 value) const{
		if(this->header.data_header.hasMixedStride() == false)
			return false;

		return(this->header.data_header.stride == value);
	}

	/**<
	 * Triggers the necessary switches to set this container
	 * as having mixed strides
	 */
	inline void triggerMixedStride(void){
		this->header.data_header.stride = -1;
		this->header.data_header.controller.mixedStride = true;
	}

	// Operators
	inline void operator++(void){ ++this->header.n_entries; }

	/**<
	 * Adds a stride value to the uncompressed buffer. At this
	 * point all stride values added must be of type U32. This
	 * function internally checks if stride sizes is mixed or
	 * not.
	 * @param value Stride value to add
	 */
	inline void addStride(const U32 value){
		// If this is the first stride set
		if(this->header.n_strides == 0){
			this->header.stride_header.controller.type = YON_TYPE_32B;
			this->header.stride_header.controller.signedness = false;
			this->setStrideSize(value);
		}

		// Check if there are different strides
		if(!this->checkStrideSize(value)){
			this->triggerMixedStride();
		}

		// Add value
		this->buffer_strides_uncompressed += (U32)value;
		++this->header.n_strides;
	}

	// Supportive
	inline const U64& getSizeUncompressed(void) const{ return(this->buffer_data_uncompressed.size()); }
	inline const U64& getSizeCompressed(void) const{ return(this->buffer_data.size()); }
	inline const U32& size(void) const{ return(this->header.n_entries); }

	inline bool __checkInteger(void){
		if(this->header.data_header.controller.encoder == 0 && this->header.n_entries == 0){
			this->header.data_header.setType(YON_TYPE_32B);
			this->header.data_header.controller.signedness = true;
		}

		// Make checks
		if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_32B, true)){
			std::cerr << "Illegal primitive type match integer!" << std::endl;
			exit(1);
			return false;
		}
		return true;
	}

	inline bool Add(const BYTE& value){
		if(!this->__checkInteger()) return false;
		this->buffer_data_uncompressed += (S32)value;
		++this->header.n_additions;
		//++this->n_entries;
		return(true);
	}

	inline bool Add(const U16& value){
		if(!this->__checkInteger()) return false;
		this->buffer_data_uncompressed += (S32)value;
		++this->header.n_additions;
		//++this->n_entries;
		return(true);
	}

	inline bool Add(const U32& value){
		if(!this->__checkInteger()) return false;
		this->buffer_data_uncompressed += (S32)value;
		++this->header.n_additions;
		//++this->n_entries;
		return(true);
	}

	inline bool Add(const SBYTE& value){
		if(!this->__checkInteger()) return false;
		this->buffer_data_uncompressed += (S32)value;
		++this->header.n_additions;
		//++this->n_entries;
		return(true);
	}

	inline bool Add(const S16& value){
		if(!this->__checkInteger()) return false;
		this->buffer_data_uncompressed += (S32)value;
		++this->header.n_additions;
		//++this->n_entries;
		return(true);
	}

	inline bool Add(const S32& value){
		if(!this->__checkInteger()) return false;
		this->buffer_data_uncompressed += (S32)value;
		++this->header.n_additions;
		//++this->n_entries;
		return(true);
	}

	inline bool Add(const U64& value){
		if(this->header.data_header.controller.encoder == 0 && this->header.n_entries == 0){
			this->header.data_header.setType(YON_TYPE_64B);
			this->header.data_header.controller.signedness = false;
			//std::cerr << "triggering: U64" << std::endl;
		}

		// Make checks
		if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_64B, false)){
			std::cerr << "Illegal primitive type match u64!" << std::endl;
			exit(1);
			return false;
		}

		this->buffer_data_uncompressed += (U64)value;
		++this->header.n_additions;
		//++this->n_entries;
		return(true);
	}

	inline bool Add(const S64& value){
		if(this->header.data_header.controller.encoder == 0 && this->header.n_entries == 0){
			this->header.data_header.setType(YON_TYPE_64B);
			this->header.data_header.controller.signedness = true;
			//std::cerr << "triggering: S64" << std::endl;
		}


		// Make checks
		if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_64B, true)){
			std::cerr << "Illegal primitive type match s64!" << std::endl;
			exit(1);
			return false;
		}

		this->buffer_data_uncompressed += (U64)value;
		++this->header.n_additions;
		//++this->n_entries;
		return(true);
	}

	inline bool Add(const float& value){
		if(this->header.data_header.controller.encoder == 0 && this->header.n_entries == 0){
			this->header.data_header.setType(YON_TYPE_FLOAT);
			this->header.data_header.controller.signedness = true;
			//std::cerr << "triggering: float" << std::endl;
		}

		// Make checks
		if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_FLOAT, true)){
			std::cerr << "Illegal primitive type match float!" << std::endl;
			exit(1);
			return false;
		}

		this->buffer_data_uncompressed += (float)value;
		++this->header.n_additions;
		//++this->n_entries;
		return(true);
	}

	inline bool Add(const double& value){
		if(this->header.data_header.controller.encoder == 0 && this->header.n_entries == 0){
			this->header.data_header.setType(YON_TYPE_DOUBLE);
			this->header.data_header.controller.signedness = true;
			//std::cerr << "triggering: float" << std::endl;
		}

		// Make checks
		if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_DOUBLE, true)){
			std::cerr << "Illegal primitive type match double!" << std::endl;
			exit(1);
			return false;
		}

		this->buffer_data_uncompressed += (double)value;
		++this->header.n_additions;
		//++this->n_entries;
		return(true);
	}

	inline bool AddCharacter(const char& value){
		if(this->header.data_header.controller.encoder == 0 && this->header.n_entries == 0){
			this->header.data_header.setType(YON_TYPE_CHAR);
			this->header.data_header.controller.signedness = true;
			std::cerr << "triggering: char" << std::endl;
		}

		// Make checks
		if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_CHAR, true)){
			std::cerr << "Illegal primitive type match char!" << std::endl;
			exit(1);
			return false;
		}

		this->buffer_data_uncompressed += (char)value;
		++this->header.n_additions;
		//++this->n_entries;
		return(true);
	}

	inline bool AddCharacter(const char* const string, const U32 l_string){
		if(this->header.data_header.controller.encoder == 0 && this->header.n_entries == 0){
			this->header.data_header.setType(YON_TYPE_CHAR);
			this->header.data_header.controller.signedness = true;
			//std::cerr << "triggering: string" << std::endl;
		}

		// Make checks
		if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_CHAR, true)){
			std::cerr << "Illegal primitive type match string!" << std::endl;
			exit(1);
			return false;
		}

		this->buffer_data_uncompressed.Add(string, l_string);
		this->header.n_additions += l_string;
		//++this->n_entries;
		return(true);
	}
	// Aliases
	inline bool AddString(const char* const string, const U32 l_string){ return(this->AddCharacter(string, l_string)); }
	inline bool AddString(const std::string& string){ return(this->AddCharacter(&string[0], string.size())); }
	inline bool AddCharacter(const std::string& string){ return(this->AddCharacter(&string[0], string.size())); }
	inline bool Add(const std::string& string){ return(this->AddCharacter(&string[0], string.size())); }

	/**<
	 *
	 * @param value
	 */
	template <class T>
	inline void AddLiteral(const T& value){
		this->buffer_data_uncompressed += (T)value;
		++this->header.n_additions;
	}

	inline void AddLiteral(const char* const string, const U32 l_string){
		this->buffer_data_uncompressed.Add(string, l_string);
		this->header.n_additions += l_string;
	}

	void reset(void);
	void resize(const U32 size);

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
	 * Checks if the current data is uniform given the provided
	 * stride size
	 * @return Returns TRUE if the data is uniform or FALSE otherwise
	 */
	bool checkUniformity(void);

	/**<
	 * This function is called during import to shrink each
	 * word-type to fit min(x) and max(x) in the worst case.
	 * At this stage all integer values in the stream is of
	 * type S32. No other values can be shrunk
	 */
	void reformat(void);

	/**<
	 * This function is caled during import to shrink each
	 * stride size element to the smallest possible primitive
	 * type to describe it without losing precision.
	 */
	void reformatStride(void);

	/**<
	 * Utility function that calculates the space this
	 * object would take on disk if written out
	 * @return Total size in bytes
	 */
	const U32 getObjectSize(void) const;

	/**< @brief Update base container header data and evaluate output byte streams
	 * Internal use only (import): Collectively updates base
	 * container offsets and checks/builds
	 * 1) If the byte stream is uniform
	 * 2) Generates CRC checksums for both data and strides
	 * 3) Reformat (change used word-size) for strides and data; if possible
	 *
	 * @param container Data container
	 * @param reormat   Reformat boolean
	 */
	void updateContainer(bool reformat = true);

	/**<
	 *
	 */
	void deltaEncode();

	inline const TACHYON_CORE_TYPE getDataPrimitiveType(void) const{ return(TACHYON_CORE_TYPE(this->header.data_header.controller.type)); }
	inline const TACHYON_CORE_TYPE getStridePrimitiveType(void) const{ return(TACHYON_CORE_TYPE(this->header.stride_header.controller.type)); }

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.buffer_data;
		if(entry.header.data_header.hasMixedStride())
			stream << entry.buffer_strides;

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		entry.buffer_data.resize(entry.header.data_header.cLength);
		stream.read(entry.buffer_data.buffer, entry.header.data_header.cLength);
		entry.buffer_data.n_chars = entry.header.data_header.cLength;

		if(entry.header.data_header.hasMixedStride()){
			entry.buffer_strides.resize(entry.header.stride_header.cLength);
			stream.read(entry.buffer_strides.buffer, entry.header.stride_header.cLength);
			entry.buffer_strides.n_chars = entry.header.stride_header.cLength;
		}
		return(stream);
	}

public:
	// Not written
	header_type header; // usually written elsewhere

	// Buffers - only bit that are written to disk from here
	buffer_type buffer_data;
	buffer_type buffer_strides;

	// These buffers are for internal use only
	buffer_type buffer_data_uncompressed;
	buffer_type buffer_strides_uncompressed;
};

}
}

#endif /* CORE_STREAMCONTAINER_H_ */
