/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TACHYON_DATA_CONTAINER_H_
#define TACHYON_DATA_CONTAINER_H_

#ifndef MD5_DIGEST_LENGTH
#define MD5_DIGEST_LENGTH 16
#endif

#include <vector>
#include <string>
#include <cstring>
#include <cassert>
#include <unordered_map>
#include <cmath>

#include "support/magic_constants.h"
#include "buffer.h"
#include "core/data_block_settings.h"

namespace tachyon {

struct yon_blk_bv_pair {
	yon_blk_bv_pair();
	~yon_blk_bv_pair();
	yon_blk_bv_pair& operator=(const yon_blk_bv_pair& other);
	yon_blk_bv_pair& operator=(yon_blk_bv_pair&& other) noexcept;

	void clear(void);

	/**<
	 * Predicate for a target local idx field in this pattern. Returns
	 * TRUE if the bit is set or FALSE otherwise. This function performs
	 * no checks for the target bit position being out-of-bounds.
	 * @param position Target bit position.
	 * @return         Returns TRUE if the bit is set or FALSE otherwise.
	 */
	inline bool operator[](const uint32_t position) const{ return((this->bit_bytes[position / 8] & (1 << (position % 8))) >> (position % 8)); }

	/**<
	 * Construct the lookup bit-vector for this object. This function needs
	 * to know the total number of fields that are set in the parent
	 * yon_vb_ftr structure as this will determine that byte-width
	 * of the bit-vector. Additionally, this function needs to be given a
	 * pointer to the map from global idx to local idx as the bit-vector
	 * bits corresponds to local idx predicates.
	 *
	 * Internally the bit-vector for all objects in a yon_vb_ftr
	 * structure has ceil(n_total_fields/8) bytes allocated for the base
	 * array.
	 * @param n_footer_total_fields Total number of fields set in the parent yon_vb_ftr.
	 * @param local_map             Pointer to map from global idx to local idx.
	 */
	void Build(const uint32_t n_footer_total_fields,
	           const std::unordered_map<uint32_t, uint32_t>* local_map);

	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const yon_blk_bv_pair& entry);
	friend yon_buffer_t& operator>>(yon_buffer_t& buffer, yon_blk_bv_pair& entry);

public:
	std::vector<int> pattern; // vector of global idx values.
	uint8_t  l_bytes; // number of bytes in bit-vector.
	uint8_t* bit_bytes; // byte array interpreted as a bit-vector.
};

// Controller type for stream container
struct yon_dc_hdr_cont {
public:
	typedef yon_dc_hdr_cont self_type;
	typedef yon_buffer_t               buffer_type;

public:
	yon_dc_hdr_cont();
	~yon_dc_hdr_cont();

	inline void clear();

	inline bool IsEncrypted(void) const{ return(this->encryption > 0); }
	inline bool CompareType(const uint8_t& type) const{ return(this->type == type); }
	inline bool CompareTypeSign(const uint8_t& type, const bool& sign) const{ return(this->type == type && this->signedness == sign); }

	self_type& operator=(const self_type& other);
	bool operator==(const self_type& other) const;
	inline bool operator!=(const self_type& other) const{ return(!(*this == other)); }

private:
	friend buffer_type& operator<<(buffer_type& buffer,const self_type& controller);
	friend std::ostream& operator<<(std::ostream& stream, const self_type& controller);
	friend std::istream& operator>>(std::istream& stream, self_type& controller);
	friend buffer_type& operator>>(buffer_type& buffer, self_type& controller);

public:
	uint32_t signedness:    1, // Signed type
	         mixedStride:   1, // Different stride sizes
	         type:          6, // Base typing (extra bits reserved for future use)
	         encoder:       5, // Encoder bits (see encoder for values)
	         uniform:       1, // Triggered if all values in the buffer are the same
	         encryption:    2, // Encryption type
			 preprocessor: 16; // preprocessor bits (extra reserved for future used)
};

struct yon_dc_hdr_obj {
public:
	typedef yon_dc_hdr_obj     self_type;
	typedef yon_dc_hdr_cont controller_type;

public:
	yon_dc_hdr_obj();
	yon_dc_hdr_obj(const self_type& other);
	yon_dc_hdr_obj(self_type&& other) noexcept;
	self_type& operator=(const self_type& other);
	self_type& operator=(self_type&& other) noexcept;
	~yon_dc_hdr_obj();

	void reset(void);

	bool operator==(const self_type& other) const;
	inline bool operator!=(const self_type& other) const{ return(!(*this == other)); }

	int8_t GetPrimitiveWidth(void) const;

	inline int32_t& GetStride(void){ return(this->stride); }
	inline const int32_t& GetStride(void) const{ return(this->stride); }
	inline bool IsUniform(void) const{ return(this->controller.uniform); }
	inline bool IsSigned(void) const{ return(this->controller.signedness); }
	inline bool HasMixedStride(void) const{ return(this->controller.mixedStride); }
	inline void SetUniform(const bool yes){ this->controller.uniform = yes; }
	inline void SetSignedness(const bool yes){ this->controller.signedness = yes; }
	inline void SetMixedStride(const bool yes){ this->controller.mixedStride = yes; }

	inline TACHYON_CORE_TYPE GetPrimitiveType(void) const{ return(TACHYON_CORE_TYPE(this->controller.type)); }
	inline TACHYON_CORE_COMPRESSION GetEncoder(void) const{ return(TACHYON_CORE_COMPRESSION(this->controller.encoder)); }

	// Set types
	inline void SetType(const TACHYON_CORE_TYPE& type){ this->controller.type = type; }

	// Checksum
	inline uint8_t* GetChecksum(void){ return(&this->crc[0]); }
	inline const uint8_t* GetChecksum(void) const{ return(&this->crc[0]); }
	bool CheckChecksum(const uint8_t* compare) const;

private:
	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const self_type& entry);
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry);
	friend yon_buffer_t& operator>>(yon_buffer_t& buffer, self_type& entry);
	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry);

public:
	controller_type controller; // controller bits
	int32_t  stride;  // stride size: -1 if not uniform, a non-zero positive value otherwise
	uint32_t offset; // relative file offset
	uint32_t cLength;// compressed length
	uint32_t uLength;// uncompressed length
	uint32_t eLength;// encrypted length
	uint8_t  crc[MD5_DIGEST_LENGTH];  // MD5 checksum
	int32_t  global_key;// global key
};

struct yon_dc_hdr {
private:
	typedef yon_dc_hdr       self_type;
	typedef yon_dc_hdr_obj header_type;
	typedef yon_buffer_t           buffer_type;

public:
	yon_dc_hdr();
	yon_dc_hdr(const self_type& other);
	~yon_dc_hdr();

	self_type& operator=(const self_type& other);
	self_type& operator=(self_type&& other) noexcept;

	void reset(void);

	// Comparators
	bool operator==(const self_type& other) const;
	inline bool operator!=(const self_type& other) const{ return(!(*this == other)); }

	self_type& operator+=(const self_type& other);

	// Accessors
	inline int32_t& GetGlobalKey(void){ return(this->data_header.global_key); }
	inline const int32_t& GetGlobalKey(void) const{ return(this->data_header.global_key); }
	inline bool HasMixedStride(void) const{ return(this->data_header.HasMixedStride()); }

private:
	friend buffer_type& operator<<(buffer_type& buffer, const self_type& entry);
	friend buffer_type& operator>>(buffer_type& buffer, self_type& entry);
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry);
	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry);

public:
	uint64_t identifier;
	uint32_t n_entries;      // number of container entries
	uint32_t n_additions;    // number of times an addition operation was executed
	uint32_t n_strides;      // number of stride elements
	header_type data_header;
	header_type stride_header;
};

/**<
 * Primary data container in Tachyon. Store actual byte
 * streams and its associated data required to restore
 * the input object.
 */
class yon1_dc_t {
public:
	typedef yon1_dc_t    self_type;
	typedef yon_dc_hdr   header_type;
	typedef yon_buffer_t buffer_type;

public:
	yon1_dc_t();
	yon1_dc_t(const uint32_t start_size);
	~yon1_dc_t();
	yon1_dc_t(self_type&& other) noexcept;
	yon1_dc_t(const self_type& other);
	yon1_dc_t& operator=(const self_type& other);
	yon1_dc_t& operator=(self_type&& other) noexcept;

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

	// Accessors for the primitive types.
	inline TACHYON_CORE_TYPE GetDataPrimitiveType(void) const{ return(TACHYON_CORE_TYPE(this->header.data_header.controller.type)); }
	inline TACHYON_CORE_TYPE GetStridePrimitiveType(void) const{ return(TACHYON_CORE_TYPE(this->header.stride_header.controller.type)); }

	// Accessors for the global identifier
	inline int32_t& GetIdx(void){ return(this->header.data_header.global_key); }
	inline const int32_t& GetIdx(void) const{ return(this->header.data_header.global_key); }
	inline int32_t& GetGlobalKey(void){ return(this->header.data_header.global_key); }
	inline const int32_t& GetGlobalKey(void) const{ return(this->header.data_header.global_key); }

	/**<
	 * Predicate for checking if the byte stream is encrypted or not.
	 * @return Returns TRUE if the byte streams is encrypted or FALSE otherwise.
	 */
	inline bool IsEncrypted(void) const{ return(this->header.data_header.controller.encryption != YON_ENCRYPTION_NONE); }

	/**<
	 * Check if the stride size of this container matches the
	 * supplied stride size. Always returns FALSE if the container
	 * has mixed stride sizes.
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
	 * as having mixed strides.
	 */
	inline void TriggerMixedStride(void){
		this->header.data_header.stride = -1;
		this->header.data_header.controller.mixedStride = true;
	}

	// Operators
	inline void operator++(void){ ++this->header.n_entries; }

	self_type& operator+=(const self_type& other){
		// Add buffers together
		this->data                 += other.data;
		this->data_uncompressed    += other.data_uncompressed;
		this->strides              += other.strides;
		this->strides_uncompressed += other.strides_uncompressed;
		// Add header counters together
		this->header += other.header;
		return(*this);
	}

	// Utility functions.
	inline const uint64_t& GetSizeUncompressed(void) const{ return(this->data_uncompressed.size()); }
	inline const uint64_t& GetSizeCompressed(void) const{ return(this->data.size()); }
	inline const uint32_t& size(void) const{ return(this->header.n_entries); }

	/**<
	 * Adds a stride value to the uncompressed buffer. At this
	 * point all stride values added must be of type uint32_t. This
	 * function internally checks if stride sizes is mixed or
	 * not.
	 * @param value Stride value to add
	 */
	void AddStride(const uint32_t value);

	/**<
	 * Adds a value to the data buffer. Integers are always added as int32_t
	 * types. This means that EOV sentinel node symbols and missing values
	 * have to be expanded into int32_t space.
	 * @param value Src value.
	 * @return      Returns TRUE upon success or FALSE otherwise.
	 */
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
	// Aliases for adding a string.
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
		this->data_uncompressed += (T)value;
		++this->header.n_additions;
	}

	inline void AddLiteral(const char* const string, const uint32_t l_string){
		this->data_uncompressed.Add(string, l_string);
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
	bool CheckUniformity(const uint32_t n_samples);

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
	void UpdateContainerFormat(bool reformat_data, bool reformat_stride, const uint32_t n_samples);

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
		if(!this->header.data_header.controller.CompareTypeSign(YON_TYPE_32B, true)){
			std::cerr << utility::timestamp("ERROR") << "Illegal primitive type mismatch (integer)!" << std::endl;
			return false;
		}
		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry);
	friend std::istream& operator>>(std::istream& stream, self_type& entry);

public:
	header_type header;
	buffer_type data;
	buffer_type strides;
	buffer_type data_uncompressed;
	buffer_type strides_uncompressed;
};

}

#endif /* CONTAINERS_DATA_CONTAINERS_H_ */
