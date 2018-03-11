#include "../components/datablock_header.h"

#include <cmath>
#include <cassert>


namespace tachyon{
namespace containers{

DataBlockHeaderBase::DataBlockHeaderBase() :
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

DataBlockHeaderBase::~DataBlockHeaderBase(){}

DataBlockHeader::DataBlockHeader():
	info_offsets(nullptr),
	format_offsets(nullptr),
	filter_offsets(nullptr),
	info_bit_vectors(nullptr),
	format_bit_vectors(nullptr),
	filter_bit_vectors(nullptr)
{}

DataBlockHeader::~DataBlockHeader(){
	delete [] this->info_offsets;
	delete [] this->format_offsets;
	delete [] this->filter_offsets;
	delete [] this->info_bit_vectors;
	delete [] this->format_bit_vectors;
	delete [] this->filter_bit_vectors;
}

void DataBlockHeader::reset(void){
	this->DataBlockHeaderBase::reset();
	this->offset_ppa.reset();
	this->offset_hot_meta.reset();
	this->offset_cold_meta.reset();
	this->offset_gt_rle.reset();
	this->offset_gt_simple.reset();
	this->offset_gt_helper.reset();
	this->offset_meta_filter_id.reset();
	this->offset_meta_format_id.reset();
	this->offset_meta_info_id.reset();

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

bool DataBlockHeader::constructBitVector(const INDEX_BLOCK_TARGET& target, hash_container_type& values, hash_vector_container_type& patterns){
	if(values.size() == 0)
		return false;

	// Determine target
	switch(target){
	case(INDEX_BLOCK_TARGET::INDEX_INFO)   :
		this->n_info_patterns = patterns.size();
		return(this->__constructBitVector(this->info_bit_vectors, this->info_offsets,  values, patterns));
		break;
	case(INDEX_BLOCK_TARGET::INDEX_FORMAT) :
		this->n_format_patterns = patterns.size();
		return(this->__constructBitVector(this->format_bit_vectors, this->format_offsets, values, patterns));
		break;
	case(INDEX_BLOCK_TARGET::INDEX_FILTER) :
		this->n_filter_patterns = patterns.size();
		return(this->__constructBitVector(this->filter_bit_vectors, this->filter_offsets, values, patterns));
		break;
	default: std::cerr << "unknown target type" << std::endl; exit(1);
	}

	return false;
}

bool DataBlockHeader::__constructBitVector(bit_vector*& target,
                                          header_type*  offset,
                                  hash_container_type&  values,
                           hash_vector_container_type& patterns)
{
	if(values.size() == 0) return false;
	BYTE bitvector_width = ceil((float)values.size()/8);
	if(bitvector_width == 0) bitvector_width = 1;

	// Create new bit-vectors
	delete [] target;
	target = new bit_vector[patterns.size()];

	// Allocate memory for these bit-vectors
	for(U32 i = 0; i < patterns.size(); ++i)
		target[i].allocate(patterns[i].size(), bitvector_width);

	// Cycle over pattern size
	for(U32 i = 0; i < patterns.size(); ++i){
		for(U32 j = 0; j < patterns[i].size(); ++j){
			U32 local_key = 0;
			// Map from absolute key to local key
			if(!values.getRaw(patterns[i][j], local_key)){
				std::cerr << "impossible to get " << patterns[i][j] << std::endl;
				exit(1);
			}

			// Flip bit at local key position
			target[i].bit_bytes[local_key/8] |= 1 << (local_key % 8);

			// Store local key in key-chain
			target[i].local_keys[j] = local_key;

			// Store absolute key
			//assert(offset[local_key].key == patterns[i][j] || offset[local_key].key == 0);
			offset[local_key].data_header.global_key = patterns[i][j];
		}
	}
	return true;
}

}
}
