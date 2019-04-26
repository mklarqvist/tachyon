#ifndef COMPRESSIONCONTAINER_H_
#define COMPRESSIONCONTAINER_H_

#include "genotypes.h"
#include "data_container.h"

namespace tachyon{
namespace algorithm{

/**< Lower bounds threshold in fold-change for compression to be kept */
#define MIN_COMPRESSION_FOLD 1.05

/**
 * Permute bits from a uint8_t-stream of uint32_t into target
 * destinations such that bit X from a
 * uint8_t [1,2,3,4,5,6,7,8] stream is permuted
 * to 1,1,1,1,1....8,8,8,8,8
 * @param data Input char* buffer
 * @param size Length of input data
 * @param destination Destination char* buffer of permuted data
 * @return TRUE if passing or FALSE otherwise
 */
inline uint32_t permuteIntBits(const char* const data,
                          const uint32_t  size,
                          char* destination)
{
	if (size == 0) return 0;
	// Balance the number of uint8_ts in the output
	// uint8_t stream to be divisible by 32. Assert
	// that this is true or the procedure fails.
	const uint32_t internal_size = size + (32-size % 32); // Balance uint8_ts
	assert(internal_size % 32 == 0);

	// Interpret the dst target as unsigned char
	// to prevent annoying signedness.
	uint8_t* dest = reinterpret_cast<uint8_t*>(destination);
	memset(dest, 0, internal_size); // Set all uint8_ts to 0
	const uint8_t* const d = reinterpret_cast<const uint8_t* const>(data); // Recast as uchar
	uint8_t* target[32]; // Bucket pointers
	const uint32_t partition_size = internal_size / 32; // Partition size

	// Assign a pointer to each bucket
	for (uint32_t i = 0; i < 32; ++i)
		target[31-i] = &dest[partition_size*i];

	uint32_t k = 0, p = 0;
	// Iterate over the data and  update position K for
	// each element. When K reaches position 7 then reset
	// to 0.
	for (uint32_t i = 0; i + 4 < internal_size; i+=4, ++k) {
		if (k == 8) { k = 0; ++p; }

		// Foreach bit in uint32_t
		// Update target T at uint8_t position P with bit J at position K
		for (uint32_t j = 0; j < 8; ++j) target[j+ 0][p] |= ((d[i]   & (1 << j)) >> j) << k;
		for (uint32_t j = 0; j < 8; ++j) target[j+ 8][p] |= ((d[i+1] & (1 << j)) >> j) << k;
		for (uint32_t j = 0; j < 8; ++j) target[j+16][p] |= ((d[i+2] & (1 << j)) >> j) << k;
		for (uint32_t j = 0; j < 8; ++j) target[j+24][p] |= ((d[i+3] & (1 << j)) >> j) << k;
	}

	return internal_size;
}

inline uint32_t unpermuteIntBits(char* data,
                            const uint32_t  size,
                            char* destination)
{
	if (size == 0) return 0;
	//uint32_t internal_size = size + (32-size%32); // Balance uint8_ts
	//assert(internal_size % 32 == 0);

	uint8_t* temp = reinterpret_cast<uint8_t*>(data); // Recast destination as uint32_t
	uint32_t* dest = reinterpret_cast<uint32_t*>(destination); // Recast destination as uint32_t
	const uint32_t n_entries = size / sizeof(uint32_t);
	memset(destination, 0, size); // Set all uint8_ts to 0
	//const uint8_t* const d = reinterpret_cast<const uint8_t* const>(data); // Recast as uchar
	uint8_t* target[32]; // Bucket pointers
	const uint32_t partition_size = size / 32; // Partition size

	// Assign a pointer to each bucket
	for (uint32_t i = 0; i < 32; ++i)
		target[31-i] = &temp[partition_size*i];

	uint32_t k = 0; uint32_t p = 0;
	// Foreach uint32_t
	// Update position K for each element
	// When K reaches position 7 then reset to 0
	for (uint32_t i = 0; i < n_entries; ++i, ++k) {
		if (k == 8) { k = 0; ++p; }

		for (uint32_t j = 0; j < 32; ++j)
			dest[i] |= ((target[j][p] & (1 << k)) >> k) << j;
	}

	//std::cerr << "out: " << internal_size << "/" << size/sizeof(uint32_t) << std::endl;
	return size;
}

/**
 * Encodes an unsigned variable-length integer using the MSB algorithm.
 * This function assumes that the value is stored as little endian.
 * @param value The input value. Any standard integer type is allowed.
 * @param output A pointer to a piece of reserved memory. Must have a minimum size dependent on the input size (32 bit = 5 bytes, 64 bit = 10 bytes).
 * @return The number of bytes used in the output memory.
 */
template<typename int_t = uint64_t>
size_t EncodeVarint(int_t value, uint8_t* output) {
    size_t outputSize = 0;
    //While more than 7 bits of data are left, occupy the last output byte
    // and set the next byte flag
    while (value > 127) {
        //|128: Set the next byte flag
        output[outputSize] = ((uint8_t)(value & 127)) | 128;
        //Remove the seven bits we just wrote
        value >>= 7;
        outputSize++;
    }
    output[outputSize++] = ((uint8_t)value) & 127;
    return outputSize;
}
/**
 * Decodes an unsigned variable-length integer using the MSB algorithm.
 * @param value A variable-length encoded integer of arbitrary size.
 * @param offset Virtual stream offset in bytes.
 */
template<typename int_t = uint64_t>
int_t DecodeVarint(uint8_t* input, size_t& offset) {
    int_t ret = 0; uint8_t its = 0;
    while(true) {
        ret |= (input[offset] & 127) << (7 * its);
        //If the next-byte flag is set
        if (!(input[offset] & 128)) {
        	++offset;
        	break;
        }
        ++offset; ++its;
    }
    return ret;
}


static uint64_t
EncodeZigzag64(int64_t input)
{
    return ((input << 1) ^ (input >> 63));
}

static int64_t
DecodeZigzag64(uint64_t input)
{
    return ((input >> 1) ^ -(input & 1));
}

static uint32_t
EncodeZigzag32(int32_t input)
{
    return ((input << 1) ^ (input >> 31));
}

static int32_t
DecodeZigzag32(uint32_t input)
{
    return ((input >> 1) ^ -(input & 1));
}

static uint16_t
EncodeZigzag16(int16_t input)
{
    return ((input << 1) ^ (input >> 15));
}

static int16_t
DecodeZigzag16(uint16_t input)
{
    return ((input >> 1) ^ -(input & 1));
}

static uint8_t
EncodeZigzag8(int8_t input)
{
    return ((input << 1) ^ (input >> 7));
}

static int8_t
DecodeZigzag8(uint8_t input)
{
    return ((input >> 1) ^ -(input & 1));
}

static size_t
EncodeInt64(int64_t input, uint8_t* output)
{
	return(EncodeVarint(EncodeZigzag64(input), output));
}

static size_t
EncodeInt32(int32_t input, uint8_t* output)
{
	return(EncodeVarint(EncodeZigzag32(input), output));
}

static size_t
EncodeInt16(int16_t input, uint8_t* output)
{
    return(EncodeVarint(EncodeZigzag16(input), output));
}



class CompressionContainer{
public:
	typedef CompressionContainer self_type;
	typedef yon1_dc_t     container_type;
	typedef yon_buffer_t  buffer_type;
	typedef yon_gt_ppa    permutation_type;

public:
	CompressionContainer() = default;
	virtual ~CompressionContainer() = default;
	virtual bool Compress(container_type& container, permutation_type& manager) =0;
	virtual bool Compress(container_type& container)          =0;
	virtual bool CompressStrides(container_type& container)   =0;
	virtual bool Decompress(container_type& container)        =0;
	virtual bool DecompressStrides(container_type& container) =0;

protected:
	buffer_type buffer;
};

}
}

#endif /* COMPRESSIONCONTAINER_H_ */
