#include <cmath>
#include "../support/helpers.h"
#include "IndexBlockEntry.h"

namespace Tomahawk{
namespace Index{

IndexBlockEntryBase::IndexBlockEntryBase() :
	contigID(-1),
	n_variants(0),
	offset_streams_begin(0),
	l_ppa(0),
	l_meta(0),
	l_meta_complex(0),
	l_gt_rle(0),
	l_gt_simple(0),
	n_info_streams(0),
	n_format_streams(0),
	n_filter_streams(0)
{}

IndexBlockEntryBase::~IndexBlockEntryBase(){}

IndexBlockEntry::IndexBlockEntry():
	info_offsets(nullptr),
	format_offsets(nullptr),
	filter_offsets(nullptr),
	info_bit_vectors(nullptr),
	format_bit_vectors(nullptr),
	filter_bit_vectors(nullptr)
{}

IndexBlockEntry::~IndexBlockEntry(){
	delete [] this->info_offsets;
	delete [] this->format_offsets;
	delete [] this->filter_offsets;
	delete [] this->info_bit_vectors;
	delete [] this->format_bit_vectors;
	delete [] this->filter_bit_vectors;
}

void IndexBlockEntry::reset(void){
	this->IndexBlockEntryBase::reset();

	delete [] this->info_offsets;
	delete [] this->format_offsets;
	delete [] this->filter_offsets;
	this->info_offsets = nullptr;
	this->format_offsets = nullptr;
	this->filter_offsets = nullptr;

	delete [] this->info_bit_vectors;
	delete [] this->format_bit_vectors;
	delete [] this->filter_bit_vectors;
	this->info_bit_vectors = nullptr;
	this->format_bit_vectors = nullptr;
	this->filter_bit_vectors = nullptr;
}

bool IndexBlockEntry::constructBitVector(const INDEX_BLOCK_TARGET& target, hash_table& htable, const id_vector& values, const pattern_vector& patterns){
	// Determine target
	switch(target){
	case(INDEX_BLOCK_TARGET::INDEX_INFO)   : return(this->__constructBitVector(this->info_bit_vectors, htable, values, patterns));   break;
	case(INDEX_BLOCK_TARGET::INDEX_FORMAT) : return(this->__constructBitVector(this->format_bit_vectors, htable, values, patterns)); break;
	case(INDEX_BLOCK_TARGET::INDEX_FILTER) : return(this->__constructBitVector(this->filter_bit_vectors, htable, values, patterns)); break;
	}

	return false;
}

bool IndexBlockEntry::__constructBitVector(bit_vector*& target, hash_table& htable, const id_vector& values, const pattern_vector& patterns){
	const BYTE bitvector_width = ceil((float)values.size()/8);

// Clear data if present
	// Allocate new bit-vectors
	if(target != nullptr){ delete [] target; }
	target = new bit_vector[patterns.size()];

	// Allocate memory for these bit-vectors
	for(U32 i = 0; i < patterns.size(); ++i)
		target[i].allocate(bitvector_width);

	// Cycle over pattern size
	for(U32 i = 0; i < patterns.size(); ++i){
		// Dump data
		std::cerr << i << '\t';
		for(U32 j = 0; j < patterns[i].size(); ++j){
			std::cerr << patterns[i][j] << '\t';
		}
		std::cerr << std::endl;

		//
		std::cerr << i << '\t';
		for(U32 j = 0; j < patterns[i].size(); ++j){
			U32* retval = nullptr;
			if(!htable.GetItem(&patterns[i][j], retval, sizeof(U32))){
				std::cerr << "impossible" << std::endl;
				exit(1);
			}
			target[i].bit_bytes[*retval/8] ^= 1 << (*retval % 8);
			std::cerr << *retval << '\t';
		}
		std::cerr << std::endl;

		std::cerr << i << '\t';
		for(U32 j = 0; j < bitvector_width; ++j)
			std::cerr << std::bitset<8>(target[i].bit_bytes[j]);

		std::cerr << std::endl << std::endl;
	}
	std::cerr << std::endl;

	return true;
}

}
}
