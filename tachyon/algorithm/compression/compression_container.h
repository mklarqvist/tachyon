#ifndef COMPRESSIONCONTAINER_H_
#define COMPRESSIONCONTAINER_H_

#include "../permutation/permutation_manager.h"
#include "../../containers/data_container.h"

namespace tachyon{
namespace algorithm{

/**< Lower bounds threshold in fold-change for compression to be kept */
#define MIN_COMPRESSION_FOLD 1.1

/**
 * Permute bits from a byte-stream of U32 into target
 * destinations such that bit X from a
 * byte [1,2,3,4,5,6,7,8] stream is permuted
 * to 1,1,1,1,1....8,8,8,8,8
 * @param data Input char* buffer
 * @param size Length of input data
 * @param destination Destination char* buffer of permuted data
 * @return TRUE if passing or FALSE otherwise
 */
inline const U32 permuteIntBits(const char* const  data,
                                        const U32  size,
                                             char* destination)
{
	if(size == 0) return 0;
	U32 internal_size = size + (32-size%32); // Balance bytes
	assert(internal_size % 32 == 0);

	BYTE* dest = reinterpret_cast<BYTE*>(destination);
	memset(dest, 0, internal_size); // Set all bytes to 0
	const BYTE* const d = reinterpret_cast<const BYTE* const>(data); // Recast as uchar
	BYTE* target[32]; // Bucket pointers
	const U32 partition_size = internal_size / 32; // Partition size

	// Assign a pointer to each bucket
	for(U32 i = 0; i < 32; ++i)
		target[31-i] = &dest[partition_size*i];

	U32 k = 0; U32 p = 0;
	// Foreach U32
	// Update position K for each element
	// When K reaches position 7 then reset to 0
	for(U32 i = 0; i + 4 < internal_size; i+=4, ++k){
		if(k == 8){ k = 0; ++p; }

		// Foreach bit in U32
		// Update target T at byte position P with bit J at position K
		for(U32 j = 0; j < 8; ++j)
			target[j][p] |= ((d[i] & (1 << j)) >> j) << k;

		for(U32 j = 0; j < 8; ++j)
			target[j+8][p] |= ((d[i+1] & (1 << j)) >> j) << k;

		for(U32 j = 0; j < 8; ++j)
			target[j+16][p] |= ((d[i+2] & (1 << j)) >> j) << k;

		for(U32 j = 0; j < 8; ++j)
			target[j+24][p] |= ((d[i+3] & (1 << j)) >> j) << k;
	}

	return internal_size;
}

inline const U32 unpermuteIntBits(char* data,
                             const U32  size,
                                  char* destination)
{
	if(size == 0) return 0;
	//U32 internal_size = size + (32-size%32); // Balance bytes
	//assert(internal_size % 32 == 0);

	BYTE* temp = reinterpret_cast<BYTE*>(data); // Recast destination as U32
	U32* dest = reinterpret_cast<U32*>(destination); // Recast destination as U32
	const U32 n_entries = size / sizeof(U32);
	memset(destination, 0, size); // Set all bytes to 0
	//const BYTE* const d = reinterpret_cast<const BYTE* const>(data); // Recast as uchar
	BYTE* target[32]; // Bucket pointers
	const U32 partition_size = size / 32; // Partition size

	// Assign a pointer to each bucket
	for(U32 i = 0; i < 32; ++i)
		target[31-i] = &temp[partition_size*i];

	/*
	for(U32 i = 0; i < size; ++i){
		std::cerr << (int)data[i] << ' ';
	}
	std::cerr << std::endl;
	*/

	U32 k = 0; U32 p = 0;
	// Foreach U32
	// Update position K for each element
	// When K reaches position 7 then reset to 0
	for(U32 i = 0; i < n_entries; ++i, ++k){
		if(k == 8){ k = 0; ++p; }

		for(U32 j = 0; j < 32; ++j){
			dest[i] |= ((target[j][p] & (1 << k)) >> k) << j;
		}
	}

	//std::cerr << "out: " << internal_size << "/" << size/sizeof(U32) << std::endl;
	return size;
}

inline const U32 permuteByteBits(const char* const  data,
                                        const U32  size,
                                             char* destination)
{
	if(size == 0) return 0;
	U32 internal_size = size + (8 - size % 8); // Balance bytes
	assert(internal_size % 8 == 0);

	BYTE* dest = reinterpret_cast<BYTE*>(destination);
	memset(dest, 0, internal_size); // Set all bytes to 0
	const BYTE* const d = reinterpret_cast<const BYTE* const>(data); // Recast as uchar
	BYTE* target[8]; // Bucket pointers
	const U32 partition_size = internal_size / 8; // Partition size
	assert(partition_size != 0);

	// Assign a pointer to each bucket
	for(U32 i = 0; i < 8; ++i)
		target[7-i] = &dest[partition_size*i];

	U32 k = 0; U32 p = 0;
	// Foreach U32
	// Update position K for each element
	// When K reaches position 7 then reset to 0
	for(U32 i = 0; i < internal_size; ++i, ++k){
		if(k == 8){ k = 0; ++p; }

		// Foreach bit in U32
		// Update bucket B at byte position P with bit J at position K
		for(U32 bucket = 0; bucket < 8; ++bucket)
			target[bucket][p] |= ((d[i] & (1 << bucket)) >> bucket) << k;
	}

	return internal_size;
}

inline const U32 unpermuteByteBits(char* data,
                             const U32  size,
                                  char* destination)
{
	if(size == 0) return 0;
	//U32 internal_size = size + (32-size%32); // Balance bytes
	//assert(internal_size % 32 == 0);

	BYTE* temp = reinterpret_cast<BYTE*>(data); // Recast destination as U32
	BYTE* dest = reinterpret_cast<BYTE*>(destination); // Recast destination as U32
	memset(destination, 0, size); // Set all bytes to 0
	//const BYTE* const d = reinterpret_cast<const BYTE* const>(data); // Recast as uchar
	BYTE* target[8]; // Bucket pointers
	const U32 partition_size = size / 8; // Partition size

	// Assign a pointer to each bucket
	for(U32 i = 0; i < 8; ++i)
		target[7-i] = &temp[partition_size*i];

	/*
	for(U32 i = 0; i < size; ++i){
		std::cerr << (int)data[i] << ' ';
	}
	std::cerr << std::endl;
	*/

	U32 k = 0; U32 p = 0;
	// Foreach U32
	// Update position K for each element
	// When K reaches position 7 then reset to 0
	for(U32 i = 0; i < size; ++i, ++k){
		if(k == 8){ k = 0; ++p; }

		for(U32 j = 0; j < 8; ++j){
			dest[i] |= ((target[j][p] & (1 << k)) >> k) << j;
		}
	}

	//std::cerr << "out: " << internal_size << "/" << size/sizeof(U32) << std::endl;
	return size;
}

class CompressionContainer{
private:
	typedef CompressionContainer          self_type;

protected:
	typedef containers::DataContainer     container_type;
	typedef io::BasicBuffer               buffer_type;
	typedef algorithm::PermutationManager permutation_type;

public:
	CompressionContainer(){}
	virtual ~CompressionContainer(){}
	virtual const bool compress(permutation_type& manager) =0;
	virtual const bool compress(container_type& container) =0;
	virtual const bool compressStrides(container_type& container) =0;
	virtual const bool decompress(container_type& container) =0;
	virtual const bool decompressStrides(container_type& container) =0;

protected:
	buffer_type buffer;
};

}
}

#endif /* COMPRESSIONCONTAINER_H_ */
