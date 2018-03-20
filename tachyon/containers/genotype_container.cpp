#include "genotype_container.h"
#include "primitive_container.h"
#include "stride_container.h"

namespace tachyon{
namespace containers{

GenotypeContainer::GenotypeContainer(const block_type& block, const MetaContainer& meta) :
	n_entries(0),
	__meta_container(block),
	__iterators(nullptr)
{
	// Todo: if anything is uniform
	// Support
	PrimitiveContainer<U32> lengths(block.gt_support_data_container); // n_runs / objects size

	U64 offset_rle8     = 0; const char* const rle8     = block.gt_rle8_container.buffer_data_uncompressed.data();
	U64 offset_rle16    = 0; const char* const rle16    = block.gt_rle16_container.buffer_data_uncompressed.data();
	U64 offset_rle32    = 0; const char* const rle32    = block.gt_rle32_container.buffer_data_uncompressed.data();
	U64 offset_rle64    = 0; const char* const rle64    = block.gt_rle64_container.buffer_data_uncompressed.data();
	U64 offset_simple8  = 0; const char* const simple8  = block.gt_simple8_container.buffer_data_uncompressed.data();
	U64 offset_simple16 = 0; const char* const simple16 = block.gt_simple16_container.buffer_data_uncompressed.data();
	U64 offset_simple32 = 0; const char* const simple32 = block.gt_simple32_container.buffer_data_uncompressed.data();
	U64 offset_simple64 = 0; const char* const simple64 = block.gt_simple64_container.buffer_data_uncompressed.data();

	assert(block.gt_rle8_container.buffer_data_uncompressed.size()     % sizeof(BYTE) == 0);
	assert(block.gt_rle16_container.buffer_data_uncompressed.size()    % sizeof(U16)  == 0);
	assert(block.gt_rle32_container.buffer_data_uncompressed.size()    % sizeof(U32)  == 0);
	assert(block.gt_rle64_container.buffer_data_uncompressed.size()    % sizeof(U64)  == 0);
	assert(block.gt_simple8_container.buffer_data_uncompressed.size()  % sizeof(BYTE) == 0);
	assert(block.gt_simple16_container.buffer_data_uncompressed.size() % sizeof(U16)  == 0);
	assert(block.gt_simple32_container.buffer_data_uncompressed.size() % sizeof(U32)  == 0);
	assert(block.gt_simple64_container.buffer_data_uncompressed.size() % sizeof(U64)  == 0);

	this->n_entries   = meta.size();
	this->__iterators = static_cast<pointer>(::operator new[](this->size() * sizeof(value_type)));

	U64 gt_offset = 0;
	for(U32 i = 0; i < meta.size(); ++i){
		if(meta[i].hasGT()){
			if(meta[i].getGenotypeEncoding() == tachyon::core::TACHYON_GT_ENCODING::YON_GT_RLE_DIPLOID_BIALLELIC){
				if(meta[i].getGenotypeType() == tachyon::core::TACHYON_GT_PRIMITIVE_TYPE::YON_GT_BYTE){
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<BYTE>( &rle8[offset_rle8], lengths[gt_offset], this->__meta_container[i] );
					offset_rle8 += lengths[gt_offset]*sizeof(BYTE);
				} else if(meta[i].getGenotypeType() == tachyon::core::TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U16){
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U16>( &rle16[offset_rle16], lengths[gt_offset], this->__meta_container[i] );
					offset_rle16 += lengths[gt_offset]*sizeof(U16);
				} else if(meta[i].getGenotypeType() == tachyon::core::TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U32){
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U32>( &rle32[offset_rle32], lengths[gt_offset], this->__meta_container[i] );
					offset_rle32 += lengths[gt_offset]*sizeof(U32);
				} else if(meta[i].getGenotypeType() == tachyon::core::TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U64){
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U64>( &rle64[offset_rle64], lengths[gt_offset], this->__meta_container[i] );
					offset_rle64 += lengths[gt_offset]*sizeof(U64);
				} else {
					std::cerr << "unknwn type" << std::endl;
					exit(1);
				}

			} else if(meta[i].getGenotypeEncoding() == tachyon::core::TACHYON_GT_ENCODING::YON_GT_RLE_DIPLOID_NALLELIC) {
				if(meta[i].getGenotypeType() == tachyon::core::TACHYON_GT_PRIMITIVE_TYPE::YON_GT_BYTE){
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<BYTE>( &simple8[offset_simple8], lengths[gt_offset], this->__meta_container[i] );
					offset_simple8 += lengths[gt_offset]*sizeof(BYTE);
				} else if(meta[i].getGenotypeType() == tachyon::core::TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U16){
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U16>( &simple16[offset_simple16], lengths[gt_offset], this->__meta_container[i] );
					offset_simple16 += lengths[gt_offset]*sizeof(U16);
				} else if(meta[i].getGenotypeType() == tachyon::core::TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U32){
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U32>( &simple32[offset_simple32], lengths[gt_offset], this->__meta_container[i] );
					offset_simple32 += lengths[gt_offset]*sizeof(U32);
				} else if(meta[i].getGenotypeType() == tachyon::core::TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U64){
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U64>( &simple64[offset_simple64], lengths[gt_offset], this->__meta_container[i] );
					offset_simple64 += lengths[gt_offset]*sizeof(U64);
				} else {
					std::cerr << "unknwn type" << std::endl;
					exit(1);
				}
			} else {
				std::cerr << "not implemented" << std::endl;
				exit(1);
			}
			++gt_offset;
		} else {
			new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<BYTE>( );
		}
	}

	assert(offset_rle8  == block.gt_rle8_container.getSizeUncompressed());
	assert(offset_rle16 == block.gt_rle16_container.getSizeUncompressed());
	assert(offset_rle32 == block.gt_rle32_container.getSizeUncompressed());
	assert(offset_rle64 == block.gt_rle64_container.getSizeUncompressed());
	assert(offset_simple8  == block.gt_simple8_container.getSizeUncompressed());
	assert(offset_simple16 == block.gt_simple16_container.getSizeUncompressed());
	assert(offset_simple32 == block.gt_simple32_container.getSizeUncompressed());
	assert(offset_simple64 == block.gt_simple64_container.getSizeUncompressed());
}

GenotypeContainer::~GenotypeContainer(){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		(this->__iterators + i)->~GenotypeContainerInterface();

	::operator delete[](static_cast<void*>(this->__iterators));
}

}
}
