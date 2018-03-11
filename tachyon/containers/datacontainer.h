#ifndef CORE_STREAMCONTAINER_H_
#define CORE_STREAMCONTAINER_H_

#include <bitset>

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
	typedef DataContainer             self_type;
	typedef io::BasicBuffer           buffer_type;
	typedef DataContainerHeader       header_type;

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
	inline const bool checkStrideSize(const U32& value) const{
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
	inline void operator++(void){ ++this->n_entries; }
	inline void addStride(const U32 value){ this->buffer_strides_uncompressed += (U32)value; }
	inline const U64& getSizeUncompressed(void) const{ return(this->buffer_data_uncompressed.size()); }
	inline const U64& getSizeCompressed(void) const{ return(this->buffer_data.size()); }
	inline const U32& size(void) const{ return(this->n_entries); }

	/**<
	 *
	 * @param value
	 */
	template <class T>
	inline void operator+=(const T& value){
		//if(this->header.controller.encoder == 0 && this->n_entries == 0)
		//	this->header.setType(value);

		/*
		// Make checks
		if(this->header.controller.checkTypeMatch(value))
		*/

		this->buffer_data_uncompressed += (T)value;
		++this->n_additions;
	}

	/**<
	 *
	 * @param string
	 */
	inline void operator+=(const std::string& string){
		if(string.size() == 0)
			return;

		++this->n_entries;
		this->buffer_data_uncompressed.Add(&string[0], string.size());
		this->n_additions += string.size();
	}

	/**
	 *
	 * @param vector
	 */
	template <class T>
	inline void operator+=(const std::vector<T>& vector){
		if(vector.size() == 0)
			return;

		++this->n_entries;
		for(U32 i = 0; i < vector.size(); ++i)
			this->buffer_data_uncompressed += vector[i];
		this->n_additions += vector.size();
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
		entry.buffer_data.resize(entry.header.data_header.uLength);
		stream.read(entry.buffer_data.buffer, entry.header.data_header.cLength);
		entry.buffer_data.n_chars = entry.header.data_header.cLength;

		if(entry.header.data_header.hasMixedStride()){
			entry.buffer_strides.resize(entry.header.stride_header.uLength);
			stream.read(entry.buffer_strides.buffer, entry.header.stride_header.cLength);
			entry.buffer_strides.n_chars = entry.header.stride_header.cLength;
		}
		return(stream);
	}

public:
	header_type header;

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
