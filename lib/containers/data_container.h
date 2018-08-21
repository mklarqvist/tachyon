#ifndef CORE_STREAMCONTAINER_H_
#define CORE_STREAMCONTAINER_H_

#include <cassert>

#include "components/data_container_header.h"
#include "components/data_container_header_controller.h"
#include "support/enums.h"
#include "io/basic_buffer.h"

namespace tachyon{
namespace containers{

/**<
 * Primary data container in Tachyon. Store actual byte
 * streams and its associated data required to restore
 * the input object.
 */
class DataContainer {
public:
	typedef DataContainer       self_type;
	typedef DataContainerHeader header_type;
	typedef io::BasicBuffer     buffer_type;

public:
	DataContainer();
	DataContainer(const uint32_t start_size);
	~DataContainer();
	DataContainer(self_type&& other) noexcept;
	DataContainer(const self_type& other);
	DataContainer& operator=(const self_type& other);
	DataContainer& operator=(self_type&& other) noexcept;

	/**<
	 * Set the primitive (or higher-order primitive) type for
	 * the objects in this container.
	 * @param value Data primitive type
	 */
	inline void SetType(const TACHYON_CORE_TYPE value){ this->header.data_header.controller.type = value; }

	/**<
	 * Set the stride size of this container to some value.
	 * The value -1 is reserved for the case when the controller
	 * flag for mixed strides is set to FALSE. Valid numbers are
	 * [0,..inf)
	 * @param value Stride size
	 */
	inline void SetStrideSize(const int32_t value){ this->header.data_header.stride = value; }

	inline TACHYON_CORE_TYPE GetDataPrimitiveType(void) const{ return(TACHYON_CORE_TYPE(this->header.data_header.controller.type)); }
	inline TACHYON_CORE_TYPE GetStridePrimitiveType(void) const{ return(TACHYON_CORE_TYPE(this->header.stride_header.controller.type)); }

	inline int32_t& GetIdx(void){ return(this->header.data_header.global_key); }
	inline const int32_t& GetIdx(void) const{ return(this->header.data_header.global_key); }


	/**<
	 * Check if the stride size of this container matches the
	 * supplied stride size. Always returns FALSE if the container
	 * has mixed stride sizes
	 * @param value Stride size to compare against
	 * @return      Returns TRUE if they are the same or FALSE otherwise
	 */
	inline bool CheckStrideSize(const int32_t value) const{
		if(this->header.data_header.HasMixedStride() == false)
			return false;

		return(this->header.data_header.stride == value);
	}

	/**<
	 * Triggers the necessary switches to set this container
	 * as having mixed strides
	 */
	inline void TriggerMixedStride(void){
		this->header.data_header.stride = -1;
		this->header.data_header.controller.mixedStride = true;
	}

	// Operators
	inline void operator++(void){ ++this->header.n_entries; }

	self_type& operator+=(const self_type& other){
		// Add buffers together
		this->buffer_data                 += other.buffer_data;
		this->buffer_data_uncompressed    += other.buffer_data_uncompressed;
		this->buffer_strides              += other.buffer_strides;
		this->buffer_strides_uncompressed += other.buffer_strides_uncompressed;
		// Add header counters together
		this->header += other.header;
		return(*this);
	}

	// Utility functions.
	inline const uint64_t& GetSizeUncompressed(void) const{ return(this->buffer_data_uncompressed.size()); }
	inline const uint64_t& GetSizeCompressed(void) const{ return(this->buffer_data.size()); }
	inline const uint32_t& size(void) const{ return(this->header.n_entries); }

	/**<
	 * Adds a stride value to the uncompressed buffer. At this
	 * point all stride values added must be of type uint32_t. This
	 * function internally checks if stride sizes is mixed or
	 * not.
	 * @param value Stride value to add
	 */
	void AddStride(const uint32_t value);

	bool Add(const uint8_t& value);
	bool Add(const uint16_t& value);
	bool Add(const uint32_t& value);
	bool Add(const int8_t& value);
	bool Add(const int16_t& value);
	bool Add(const int32_t& value);
	bool Add(const uint64_t& value);
	bool Add(const int64_t& value);
	bool Add(const float& value);
	bool Add(const double& value);
	bool AddCharacter(const char& value);
	bool AddCharacter(const char* const string, const uint32_t l_string);
	// Aliases
	inline bool AddString(const char* const string, const uint32_t l_string){ return(this->AddCharacter(string, l_string)); }
	inline bool AddString(const std::string& string){ return(this->AddCharacter(&string[0], string.size())); }
	inline bool AddCharacter(const std::string& string){ return(this->AddCharacter(&string[0], string.size())); }
	inline bool Add(const std::string& string){ return(this->AddCharacter(&string[0], string.size())); }

	/**<
	 * Add a literal primitive/string to the byte stream without conversion
	 * and without performing any checks. This is function is only used
	 * when storing explicit data structures into the byte stream and not
	 * arrays of values. This function should generally not be called by
	 * end users.
	 * @param value Src primitive type.
	 */
	template <class T>
	inline void AddLiteral(const T& value){
		this->buffer_data_uncompressed += (T)value;
		++this->header.n_additions;
	}

	inline void AddLiteral(const char* const string, const uint32_t l_string){
		this->buffer_data_uncompressed.Add(string, l_string);
		this->header.n_additions += l_string;
	}

	void reset(void);
	void resize(const uint32_t size);

	/**<
	 * Generates a MD5 checksum of the uncompressed
	 * data and, if set, the uncompressed strides data.
	 * The MD5 checksums are stored in the header.
	 *
	 * Checksums for compressed buffers are computed
	 * by the appropriate compression wrapper.
	 */
	void GenerateMd5(void);

	/**<
	 * Calculates the MD5 checksum from the target data
	 * buffer and compares that value to the expected
	 * value stored in the data header. MD5 checksums
	 * for compressed buffers are always checked during
	 * decompression to ascertain correctness.
	 *
	 * Targets:
	 * 0: Uncompressed data
	 * 1: Uncompressed strides data
	 * 2: Compressed data
	 * 3: Compressed strides data
	 *
	 * @param target Target buffer stream
	 * @return       Returns TRUE if the MD5 checksums are identical or FALSE otherwise
	 */
	bool CheckMd5(int target = 0);

	/**<
	 * Checks if the current data is uniform given the provided
	 * stride size.
	 * @return Returns TRUE if the data is uniform or FALSE otherwise
	 */
	bool CheckUniformity(void);

	/**<
	 * This function is called during import to shrink each
	 * word-type to fit min(x) and max(x) in the worst case.
	 * At this stage all integer values in the stream is of
	 * type int32_t. No other values can be shrunk.
	 */
	void ReformatInteger(void);

	/**<
	 * This function is caled during import to shrink each
	 * stride size element to the smallest possible primitive
	 * type to describe it without losing precision.
	 */
	void ReformatStride(void);

	/**<
	 * Utility function that calculates the number of bytes this
	 * object would occupy if written to a stream.
	 * @return Total number of bytes.
	 */
	uint32_t GetObjectSize(void) const;

	/**<
	 * Utility function that calculates the number of actual bytes
	 * this object occupies internally when uncompressed.
	 * @return Number of uncompressed bytes.
	 */
	uint64_t GetObjectSizeUncompressed(void) const;

	/**<
	 * Update base container header data and evaluate output byte streams.
	 * Collectively updates base container offsets and checks/builds:
	 * 1) If the byte stream is uniform;
	 * 2) Generates CRC checksums for both data and strides;
	 * 3) Reformat (change used word-size) for strides and data; if possible.
	 *
	 * @param reformat_data   Flag for whether data should be reformatted.
	 * @param reformat_stride Flag for whether strides should be reformatted.
	 */
	void UpdateContainer(bool reformat_data = true, bool reformat_stride = true);

private:
	/**<
	 * Predicate for uniformity of primitive family type. If integers are
	 * stored in a byte stream then only other integers and not for example
	 * floats, doubles, or strings may be stored in the same stream.
	 * @return Returns TRUE if acceptable addition or FALSE otherwise.
	 */
	inline bool CheckInteger(void){
		if(this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0){
			this->header.data_header.SetType(YON_TYPE_32B);
			this->header.data_header.controller.signedness = true;
		}

		// Make checks
		if(!this->header.data_header.controller.compareTypeSign(YON_TYPE_32B, true)){
			std::cerr << utility::timestamp("ERROR") << "Illegal primitive type mismatch (integer)!" << std::endl;
			return false;
		}
		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry);
	friend std::istream& operator>>(std::istream& stream, self_type& entry);

public:
	header_type header;
	buffer_type buffer_data;
	buffer_type buffer_strides;
	buffer_type buffer_data_uncompressed;
	buffer_type buffer_strides_uncompressed;
};

}
}

#endif /* CORE_STREAMCONTAINER_H_ */
