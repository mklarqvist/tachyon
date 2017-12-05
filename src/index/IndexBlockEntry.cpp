#include <cmath>
#include "IndexBlockEntry.h"

namespace Tachyon{
namespace Index{

IndexBlockEntryBase::IndexBlockEntryBase() :
	offset_end_of_block(0),
	contigID(-1),
	minPosition(0),
	maxPosition(0),
	n_variants(0),
	n_info_streams(0),
	n_format_streams(0),
	n_filter_streams(0),
	n_info_patterns(0),
	n_format_patterns(0),
	n_filter_patterns(0),
	l_info_bitvector(0),
	l_format_bitvector(0),
	l_filter_bitvector(0)
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

bool IndexBlockEntry::constructBitVector(const INDEX_BLOCK_TARGET& target, hash_container_type& values, hash_vector_container_type& patterns){
	if(values.size() == 0)
		return false;

	// Determine target
	switch(target){
	case(INDEX_BLOCK_TARGET::INDEX_INFO)   :
		this->n_info_patterns = patterns.size();
		return(this->__constructBitVector(this->info_bit_vectors, values, patterns));
		break;
	case(INDEX_BLOCK_TARGET::INDEX_FORMAT) :
		this->n_format_patterns = patterns.size();
		return(this->__constructBitVector(this->format_bit_vectors, values, patterns));
		break;
	case(INDEX_BLOCK_TARGET::INDEX_FILTER) :
		this->n_filter_patterns = patterns.size();
		return(this->__constructBitVector(this->filter_bit_vectors, values, patterns));
		break;
	default: std::cerr << "unknown target type" << std::endl; exit(1);
	}

	return false;
}

bool IndexBlockEntry::__constructBitVector(bit_vector*& target, hash_container_type& values, hash_vector_container_type& patterns){
	BYTE bitvector_width = ceil((float)values.size()/8);
	if(values.size() == 1) bitvector_width = 1;

	// Clear data if present
	// Allocate new bit-vectors
	delete [] target;
	target = new bit_vector[patterns.size()];

	// Allocate memory for these bit-vectors
	for(U32 i = 0; i < patterns.size(); ++i)
		target[i].allocate(patterns[i].size(), bitvector_width);

	// Cycle over pattern size
	for(U32 i = 0; i < patterns.size(); ++i){
		/*
		// Dump data
		std::cerr << i << '\t';
		for(U32 j = 0; j < patterns[i].size(); ++j){
			std::cerr << patterns[i][j] << '\t';
		}
		std::cerr << std::endl;
		*/

		//
		//std::cerr << i << '\t';
		for(U32 j = 0; j < patterns[i].size(); ++j){
			U32 retval = 0;
			if(!values.getRaw(patterns[i][j], retval)){
				std::cerr << "impossible to get " << patterns[i][j] << std::endl;
				exit(1);
			}
			target[i].bit_bytes[retval/8] ^= 1 << (retval % 8);
			//std::cerr << retval << '\t';
			target[i].keys[j] = retval;
			//std::cerr << target[i].keys[j] << std::endl;
		}
		//std::cerr << std::endl;

		//std::cerr << i << '\t';
		//for(U32 j = 0; j < bitvector_width; ++j)
		//	std::cerr << std::bitset<8>(target[i].bit_bytes[j]);

		//std::cerr << std::endl << std::endl;
	}
	//std::cerr << std::endl;

	return true;
}

}
}
