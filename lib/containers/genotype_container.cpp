#include "containers/genotype_container.h"
#include "containers/primitive_container.h"
#include "containers/stride_container.h"
#include "support/enums.h"

namespace tachyon{
namespace containers{

GenotypeContainer::GenotypeContainer(const block_type& block, const meta_container_type& meta) :
	n_entries(0),
	__meta_container(block),
	__iterators(nullptr)
{
	// Todo: if anything is uniform
	// Support
	const bool uniform_stride = block.base_containers[YON_BLK_GT_SUPPORT].header.data_header.IsUniform();
	PrimitiveContainer<uint32_t> lengths(block.base_containers[YON_BLK_GT_SUPPORT]); // n_runs / objects size

	uint64_t offset_rle8     = 0; const char* const rle8     = block.base_containers[YON_BLK_GT_INT8].buffer_data_uncompressed.data();
	uint64_t offset_rle16    = 0; const char* const rle16    = block.base_containers[YON_BLK_GT_INT16].buffer_data_uncompressed.data();
	uint64_t offset_rle32    = 0; const char* const rle32    = block.base_containers[YON_BLK_GT_INT32].buffer_data_uncompressed.data();
	uint64_t offset_rle64    = 0; const char* const rle64    = block.base_containers[YON_BLK_GT_INT64].buffer_data_uncompressed.data();

	uint64_t offset_simple8  = 0; const char* const simple8  = block.base_containers[YON_BLK_GT_S_INT8].buffer_data_uncompressed.data();
	uint64_t offset_simple16 = 0; const char* const simple16 = block.base_containers[YON_BLK_GT_S_INT16].buffer_data_uncompressed.data();
	uint64_t offset_simple32 = 0; const char* const simple32 = block.base_containers[YON_BLK_GT_S_INT32].buffer_data_uncompressed.data();
	uint64_t offset_simple64 = 0; const char* const simple64 = block.base_containers[YON_BLK_GT_S_INT64].buffer_data_uncompressed.data();

	uint64_t offset_nploid8  = 0; const char* const nploid8  = block.base_containers[YON_BLK_GT_N_INT8].buffer_data_uncompressed.data();
	uint64_t offset_nploid16 = 0; const char* const nploid16 = block.base_containers[YON_BLK_GT_N_INT16].buffer_data_uncompressed.data();
	uint64_t offset_nploid32 = 0; const char* const nploid32 = block.base_containers[YON_BLK_GT_N_INT32].buffer_data_uncompressed.data();
	uint64_t offset_nploid64 = 0; const char* const nploid64 = block.base_containers[YON_BLK_GT_N_INT64].buffer_data_uncompressed.data();

	assert(block.base_containers[YON_BLK_GT_INT8].buffer_data_uncompressed.size()    % sizeof(uint8_t) == 0);
	assert(block.base_containers[YON_BLK_GT_INT16].buffer_data_uncompressed.size()   % sizeof(uint16_t)  == 0);
	assert(block.base_containers[YON_BLK_GT_INT32].buffer_data_uncompressed.size()   % sizeof(uint32_t)  == 0);
	assert(block.base_containers[YON_BLK_GT_INT64].buffer_data_uncompressed.size()   % sizeof(uint64_t)  == 0);
	assert(block.base_containers[YON_BLK_GT_S_INT8].buffer_data_uncompressed.size()  % sizeof(uint8_t) == 0);
	assert(block.base_containers[YON_BLK_GT_S_INT16].buffer_data_uncompressed.size() % sizeof(uint16_t)  == 0);
	assert(block.base_containers[YON_BLK_GT_S_INT32].buffer_data_uncompressed.size() % sizeof(uint32_t)  == 0);
	assert(block.base_containers[YON_BLK_GT_S_INT64].buffer_data_uncompressed.size() % sizeof(uint64_t)  == 0);

	this->n_entries   = meta.size();
	this->__iterators = static_cast<pointer>(::operator new[](this->size() * sizeof(value_type)));

	uint64_t gt_offset = 0;
	uint8_t incrementor = 1;
	if(uniform_stride) incrementor = 0;

	for(uint32_t i = 0; i < meta.size(); ++i){
		if(meta[i].HasGT()){
			// Case run-length encoding diploid and biallelic and no missing
			if(meta[i].GetGenotypeEncoding() == TACHYON_GT_ENCODING::YON_GT_RLE_DIPLOID_BIALLELIC){
				if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_BYTE){
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<uint8_t>( &rle8[offset_rle8], lengths[gt_offset], this->__meta_container[i] );
					offset_rle8 += lengths[gt_offset]*sizeof(uint8_t);
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U16){
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<uint16_t>( &rle16[offset_rle16], lengths[gt_offset], this->__meta_container[i] );
					offset_rle16 += lengths[gt_offset]*sizeof(uint16_t);
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U32){
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<uint32_t>( &rle32[offset_rle32], lengths[gt_offset], this->__meta_container[i] );
					offset_rle32 += lengths[gt_offset]*sizeof(uint32_t);
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U64){
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<uint64_t>( &rle64[offset_rle64], lengths[gt_offset], this->__meta_container[i] );
					offset_rle64 += lengths[gt_offset]*sizeof(uint64_t);
				} else {
					std::cerr << utility::timestamp("ERROR","GT") << "Unknown GT encoding primitive..." << std::endl;
					exit(1);
				}

			}
			// Case run-length encoding diploid and biallelic/EOV or n-allelic
			else if(meta[i].GetGenotypeEncoding() == TACHYON_GT_ENCODING::YON_GT_RLE_DIPLOID_NALLELIC) {
				if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_BYTE){
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<uint8_t>( &simple8[offset_simple8], lengths[gt_offset], this->__meta_container[i] );
					offset_simple8 += lengths[gt_offset]*sizeof(uint8_t);
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U16){
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<uint16_t>( &simple16[offset_simple16], lengths[gt_offset], this->__meta_container[i] );
					offset_simple16 += lengths[gt_offset]*sizeof(uint16_t);
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U32){
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<uint32_t>( &simple32[offset_simple32], lengths[gt_offset], this->__meta_container[i] );
					offset_simple32 += lengths[gt_offset]*sizeof(uint32_t);
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U64){
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<uint64_t>( &simple64[offset_simple64], lengths[gt_offset], this->__meta_container[i] );
					offset_simple64 += lengths[gt_offset]*sizeof(uint64_t);
				} else {
					std::cerr << utility::timestamp("ERROR","GT") << "Unknown GT encoding primitive..." << std::endl;
					exit(1);
				}
			}
			// Case BCF-style encoding of diploids
			else if(meta[i].GetGenotypeEncoding() == TACHYON_GT_ENCODING::YON_GT_BCF_DIPLOID) {
				if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_BYTE){
					new( &this->__iterators[i] ) GenotypeContainerDiploidBCF<uint8_t>( &simple8[offset_simple8], lengths[gt_offset], this->__meta_container[i] );
					offset_simple8 += lengths[gt_offset]*sizeof(uint8_t);
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U16){
					new( &this->__iterators[i] ) GenotypeContainerDiploidBCF<uint16_t>( &simple16[offset_simple16], lengths[gt_offset], this->__meta_container[i] );
					offset_simple16 += lengths[gt_offset]*sizeof(uint16_t);
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U32){
					new( &this->__iterators[i] ) GenotypeContainerDiploidBCF<uint32_t>( &simple32[offset_simple32], lengths[gt_offset], this->__meta_container[i] );
					offset_simple32 += lengths[gt_offset]*sizeof(uint32_t);
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U64){
					new( &this->__iterators[i] ) GenotypeContainerDiploidBCF<uint64_t>( &simple64[offset_simple64], lengths[gt_offset], this->__meta_container[i] );
					offset_simple64 += lengths[gt_offset]*sizeof(uint64_t);
				}  else {
					std::cerr << utility::timestamp("ERROR","GT") << "Unknown GT encoding primitive..." << std::endl;
					exit(1);
				}
			}
			// Case RLE-encoding of nploids
			else if(meta[i].GetGenotypeEncoding() == TACHYON_GT_ENCODING::YON_GT_RLE_NPLOID) {
				if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_BYTE){
					new( &this->__iterators[i] ) GenotypeContainerNploid<uint8_t>( &nploid8[offset_nploid8], lengths[gt_offset], this->__meta_container[i] );
					offset_nploid8 += lengths[gt_offset]*(sizeof(uint8_t) + this->__meta_container[i].n_base_ploidy*sizeof(uint8_t));
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U16){
					new( &this->__iterators[i] ) GenotypeContainerNploid<uint16_t>( &nploid16[offset_nploid16], lengths[gt_offset], this->__meta_container[i] );
					offset_nploid16 += lengths[gt_offset]*(sizeof(uint16_t) + this->__meta_container[i].n_base_ploidy*sizeof(uint8_t));
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U32){
					new( &this->__iterators[i] ) GenotypeContainerNploid<uint32_t>( &nploid32[offset_nploid32], lengths[gt_offset], this->__meta_container[i] );
					offset_nploid32 += lengths[gt_offset]*(sizeof(uint32_t) + this->__meta_container[i].n_base_ploidy*sizeof(uint8_t));
				} else if(meta[i].GetGenotypeType() == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U64){
					new( &this->__iterators[i] ) GenotypeContainerNploid<uint64_t>( &nploid64[offset_nploid64], lengths[gt_offset], this->__meta_container[i] );
					offset_nploid64 += lengths[gt_offset]*(sizeof(uint64_t) + this->__meta_container[i].n_base_ploidy*sizeof(uint8_t));
				}  else {
					std::cerr << utility::timestamp("ERROR","GT") << "Unknown GT encoding primitive..." << std::endl;
					exit(1);
				}
			}
			// Case other potential encodings
			else {
				std::cerr << utility::timestamp("ERROR","GT") << "Unknown GT encoding family..." << std::endl;
				new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<uint8_t>( );
				exit(1);
			}

			// Increment offset
			gt_offset += incrementor;

		} else { // No GT available
			new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<uint8_t>( );
		}
	}

	assert(offset_rle8     == block.base_containers[YON_BLK_GT_INT8].GetSizeUncompressed());
	assert(offset_rle16    == block.base_containers[YON_BLK_GT_INT16].GetSizeUncompressed());
	assert(offset_rle32    == block.base_containers[YON_BLK_GT_INT32].GetSizeUncompressed());
	assert(offset_rle64    == block.base_containers[YON_BLK_GT_INT64].GetSizeUncompressed());
	assert(offset_simple8  == block.base_containers[YON_BLK_GT_S_INT8].GetSizeUncompressed());
	assert(offset_simple16 == block.base_containers[YON_BLK_GT_S_INT16].GetSizeUncompressed());
	assert(offset_simple32 == block.base_containers[YON_BLK_GT_S_INT32].GetSizeUncompressed());
	assert(offset_simple64 == block.base_containers[YON_BLK_GT_S_INT64].GetSizeUncompressed());
	assert(offset_nploid8  == block.base_containers[YON_BLK_GT_N_INT8].GetSizeUncompressed());
	assert(offset_nploid16 == block.base_containers[YON_BLK_GT_N_INT16].GetSizeUncompressed());
	assert(offset_nploid32 == block.base_containers[YON_BLK_GT_N_INT32].GetSizeUncompressed());
	assert(offset_nploid64 == block.base_containers[YON_BLK_GT_N_INT64].GetSizeUncompressed());
}

GenotypeContainer::~GenotypeContainer(){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		(this->__iterators + i)->~GenotypeContainerInterface();

	::operator delete[](static_cast<void*>(this->__iterators));
}

}
}
