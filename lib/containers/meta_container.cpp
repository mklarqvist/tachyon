#include "meta_container.h"
#include "primitive_container.h"
#include "stride_container.h"

namespace tachyon{
namespace containers{

MetaContainer::MetaContainer(const block_type& block) :
	n_entries(block.header.n_variants),
	__entries(static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type))))
{
	this->__ctor_setup(block);
	this->correctRelativePositions(block.header);
}

MetaContainer::~MetaContainer(void){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		(this->__entries + i)->~MetaEntry();

	::operator delete[](static_cast<void*>(this->__entries));
}

void MetaContainer::__ctor_setup(const block_type& block){
	// Build containers for hot/cold depending on what is available
	PrimitiveContainer<uint32_t>   contigs(block.base_containers[YON_BLK_CONTIG]);
	PrimitiveContainer<uint16_t>   controllers(block.base_containers[YON_BLK_CONTROLLER]);
	PrimitiveContainer<uint32_t>   positions(block.base_containers[YON_BLK_POSITION]);
	PrimitiveContainer<float> quality(block.base_containers[YON_BLK_QUALITY]);
	PrimitiveContainer<uint8_t>  refalt(block.base_containers[YON_BLK_REFALT]);
	PrimitiveContainer<int32_t>   filterID(block.base_containers[YON_BLK_ID_FILTER]);
	PrimitiveContainer<int32_t>   formatID(block.base_containers[YON_BLK_ID_FORMAT]);
	PrimitiveContainer<int32_t>   infoID(block.base_containers[YON_BLK_ID_INFO]);
	PrimitiveContainer<uint8_t>  ploidy(block.base_containers[YON_BLK_GT_PLOIDY]);

	for(uint32_t i = 0; i < this->size(); ++i){
		new( &this->__entries[i] ) value_type( );
	}

	if(contigs.size()){
		for(uint32_t i = 0; i < this->size(); ++i){
			if(contigs.isUniform()) this->__entries[i].contigID = contigs[0];
			else this->__entries[i].contigID = contigs[i];
		}
	}

	if(positions.size()){
		for(uint32_t i = 0; i < this->size(); ++i){
			if(positions.isUniform()) this->__entries[i].position = positions[0];
			else this->__entries[i].position = positions[i];
		}
	}

	if(controllers.size()){
		for(uint32_t i = 0; i < this->size(); ++i){
			if(controllers.isUniform()) this->__entries[i].controller = controllers[0];
			else this->__entries[i].controller = controllers[i];
		}
	}

	if(quality.size()){
		for(uint32_t i = 0; i < this->size(); ++i){
			if(quality.isUniform()) this->__entries[i].quality = quality[0];
			else this->__entries[i].quality = quality[i];
		}
	}

	if(filterID.size()){
		for(uint32_t i = 0; i < this->size(); ++i){
			if(filterID.isUniform()) this->__entries[i].filter_pattern_id = filterID[0];
			else this->__entries[i].filter_pattern_id = filterID[i];
		}
	}

	if(infoID.size()){
		for(uint32_t i = 0; i < this->size(); ++i){
			if(infoID.isUniform()) this->__entries[i].info_pattern_id = infoID[0];
			else this->__entries[i].info_pattern_id = infoID[i];
		}
	}

	if(formatID.size()){
		for(uint32_t i = 0; i < this->size(); ++i){
			if(formatID.isUniform()) this->__entries[i].format_pattern_id = formatID[0];
			else 	this->__entries[i].format_pattern_id = formatID[i];
		}
	}

	if(refalt.size()){ // execute only if we have data
		uint32_t refalt_position = 0;
		for(uint32_t i = 0; i < this->size(); ++i){
			if(this->__entries[i].controller.alleles_packed){
				// load from special packed
				// this is always diploid
				this->__entries[i].n_alleles = 2;
				this->__entries[i].alleles   = static_cast<value_type::allele_type*>(::operator new[](2*sizeof(value_type::allele_type)));

				// If data is <non_ref> or not
				if((refalt[refalt_position] & 15) != 5){
					//assert((refalt[refalt_position] & 15) < 5);
					const char ref = constants::REF_ALT_LOOKUP[refalt[refalt_position] & 15];
					new( &this->__entries[i].alleles[1] ) value_type::allele_type( ref );
				} else {
					const std::string s = "<NON_REF>";
					new( &this->__entries[i].alleles[1] ) value_type::allele_type( s );
				}

				// If data is <non_ref> or not
				if(((refalt[refalt_position] >> 4) & 15) != 5){
					//assert(((refalt[refalt_position] >> 4) & 15) < 5);
					const char alt = constants::REF_ALT_LOOKUP[(refalt[refalt_position] >> 4) & 15];
					new( &this->__entries[i].alleles[0] ) value_type::allele_type( alt );
				} else {
					const std::string s = "<NON_REF>";
					new( &this->__entries[i].alleles[0] ) value_type::allele_type( s );
				}
				// Do not increment in case this data is uniform
				if(refalt.isUniform() == false) ++refalt_position;
			}
			// otherwise load from literal cold
			else {
				// number of alleles is parsed from the stride container
			}
		} // end loop

		if(refalt.isUniform()) refalt_position = 1;
		assert(refalt_position == refalt.size());
	}

	if(block.base_containers[YON_BLK_ALLELES].buffer_data_uncompressed.size()){
		StrideContainer<uint32_t> strides(block.base_containers[YON_BLK_ALLELES]);
		uint32_t offset = 0;
		uint32_t stride_offset = 0;
		for(uint32_t i = 0; i < this->size(); ++i){
			if(this->__entries[i].controller.alleles_packed == false){
				this->__entries[i].n_alleles = strides[stride_offset];
				this->__entries[i].alleles   = static_cast<value_type::allele_type*>(::operator new[](strides[stride_offset]*sizeof(value_type::allele_type)));

				for(uint32_t j = 0; j < strides[stride_offset]; ++j){
					const uint16_t& l_string = *reinterpret_cast<const uint16_t* const>(&block.base_containers[YON_BLK_ALLELES].buffer_data_uncompressed[offset]);
					new( &this->__entries[i].alleles[j] ) value_type::allele_type( &block.base_containers[YON_BLK_ALLELES].buffer_data_uncompressed[offset] );
					offset += sizeof(uint16_t) + l_string;
				}
				++stride_offset;
			}
		}
		assert(offset == block.base_containers[YON_BLK_ALLELES].GetSizeUncompressed());
	}

	// Parse name
	if(block.base_containers[YON_BLK_NAMES].GetSizeUncompressed()){
		StrideContainer<uint32_t> strides(block.base_containers[YON_BLK_NAMES]);
		uint32_t offset = 0;
		assert(strides.size() == this->size());
		for(uint32_t i = 0; i < this->size(); ++i){
			this->__entries[i].name = std::string(&block.base_containers[YON_BLK_NAMES].buffer_data_uncompressed.data()[offset], strides[i]);
			offset += strides[i];
		}
	}

	// Parse ploidy.
	if(ploidy.size()){
		if(ploidy.isUniform()){
			for(uint32_t i = 0; i < this->size(); ++i){
				this->at(i).n_base_ploidy = ploidy[0];
			}
		} else {
			assert(ploidy.size() == this->size());
			for(uint32_t i = 0; i < this->size(); ++i){
				this->at(i).n_base_ploidy = ploidy[i];
				std::cerr << (int)ploidy[i] << std::endl;
			}
		}
	}
}

}
}
