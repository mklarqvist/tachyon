#ifndef COMPRESSIONCONTAINER_H_
#define COMPRESSIONCONTAINER_H_

#include "core/genotypes.h"
#include "containers/data_container.h"

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
	if(size == 0) return 0;
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
	for(uint32_t i = 0; i < 32; ++i)
		target[31-i] = &dest[partition_size*i];

	uint32_t k = 0, p = 0;
	// Iterate over the data and  update position K for
	// each element. When K reaches position 7 then reset
	// to 0.
	for(uint32_t i = 0; i + 4 < internal_size; i+=4, ++k){
		if(k == 8){ k = 0; ++p; }

		// Foreach bit in uint32_t
		// Update target T at uint8_t position P with bit J at position K
		for(uint32_t j = 0; j < 8; ++j) target[j+ 0][p] |= ((d[i]   & (1 << j)) >> j) << k;
		for(uint32_t j = 0; j < 8; ++j) target[j+ 8][p] |= ((d[i+1] & (1 << j)) >> j) << k;
		for(uint32_t j = 0; j < 8; ++j) target[j+16][p] |= ((d[i+2] & (1 << j)) >> j) << k;
		for(uint32_t j = 0; j < 8; ++j) target[j+24][p] |= ((d[i+3] & (1 << j)) >> j) << k;
	}

	return internal_size;
}

inline uint32_t unpermuteIntBits(char* data,
                            const uint32_t  size,
                            char* destination)
{
	if(size == 0) return 0;
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
	for(uint32_t i = 0; i < 32; ++i)
		target[31-i] = &temp[partition_size*i];

	uint32_t k = 0; uint32_t p = 0;
	// Foreach uint32_t
	// Update position K for each element
	// When K reaches position 7 then reset to 0
	for(uint32_t i = 0; i < n_entries; ++i, ++k){
		if(k == 8){ k = 0; ++p; }

		for(uint32_t j = 0; j < 32; ++j)
			dest[i] |= ((target[j][p] & (1 << k)) >> k) << j;
	}

	//std::cerr << "out: " << internal_size << "/" << size/sizeof(uint32_t) << std::endl;
	return size;
}

inline uint32_t permuteuint8_tBits(const char* const  data,
                           const uint32_t  size,
                           char* destination)
{
	if(size == 0) return 0;
	uint32_t internal_size = size + (8 - size % 8); // Balance uint8_ts
	assert(internal_size % 8 == 0);

	uint8_t* dest = reinterpret_cast<uint8_t*>(destination);
	memset(dest, 0, internal_size); // Set all uint8_ts to 0
	const uint8_t* const d = reinterpret_cast<const uint8_t* const>(data); // Recast as uchar
	uint8_t* target[8]; // Bucket pointers
	const uint32_t partition_size = internal_size / 8; // Partition size
	assert(partition_size != 0);

	// Assign a pointer to each bucket
	for(uint32_t i = 0; i < 8; ++i)
		target[7-i] = &dest[partition_size*i];

	uint32_t k = 0; uint32_t p = 0;
	// Foreach uint32_t
	// Update position K for each element
	// When K reaches position 7 then reset to 0
	for(uint32_t i = 0; i < internal_size; ++i, ++k){
		if(k == 8){ k = 0; ++p; }

		// Foreach bit in uint32_t
		// Update bucket B at uint8_t position P with bit J at position K
		for(uint32_t bucket = 0; bucket < 8; ++bucket)
			target[bucket][p] |= ((d[i] & (1 << bucket)) >> bucket) << k;
	}

	return internal_size;
}

inline uint32_t unpermuteuint8_tBits(char* data,
                             const uint32_t  size,
                             char* destination)
{
	if(size == 0) return 0;
	//uint32_t internal_size = size + (32-size%32); // Balance uint8_ts
	//assert(internal_size % 32 == 0);

	uint8_t* temp = reinterpret_cast<uint8_t*>(data); // Recast destination as uint32_t
	uint8_t* dest = reinterpret_cast<uint8_t*>(destination); // Recast destination as uint32_t
	memset(destination, 0, size); // Set all uint8_ts to 0
	//const uint8_t* const d = reinterpret_cast<const uint8_t* const>(data); // Recast as uchar
	uint8_t* target[8]; // Bucket pointers
	const uint32_t partition_size = size / 8; // Partition size

	// Assign a pointer to each bucket
	for(uint32_t i = 0; i < 8; ++i)
		target[7-i] = &temp[partition_size*i];

	/*
	for(uint32_t i = 0; i < size; ++i){
		std::cerr << (int)data[i] << ' ';
	}
	std::cerr << std::endl;
	*/

	uint32_t k = 0; uint32_t p = 0;
	// Foreach uint32_t
	// Update position K for each element
	// When K reaches position 7 then reset to 0
	for(uint32_t i = 0; i < size; ++i, ++k){
		if(k == 8){ k = 0; ++p; }

		for(uint32_t j = 0; j < 8; ++j){
			dest[i] |= ((target[j][p] & (1 << k)) >> k) << j;
		}
	}

	//std::cerr << "out: " << internal_size << "/" << size/sizeof(uint32_t) << std::endl;
	return size;
}

class CompressionContainer{
public:
	typedef CompressionContainer          self_type;
	typedef containers::DataContainer     container_type;
	typedef io::BasicBuffer               buffer_type;
	typedef yon_gt_ppa                    permutation_type;

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
