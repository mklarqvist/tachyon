#include "meta_container.h"
#include "stride_container.h"
#include "primitive_container.h"

namespace tachyon{
namespace containers{

MetaContainer::MetaContainer(const block_type& block) :
	n_entries(block.header.n_variants),
	__entries(static_cast<pointer>(::operator new[](this->n_entries*sizeof(value_type))))
{
	this->__ctor_setup(block);
}

MetaContainer::~MetaContainer(void){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		(this->__entries + i)->~MetaEntry();

	::operator delete[](static_cast<void*>(this->__entries));
}

void MetaContainer::__ctor_setup(const block_type& block){
	// Build containers for hot/cold depending on what is available
	PrimitiveContainer<U32>   contigs(block.meta_contig_container);
	PrimitiveContainer<U16>   controllers(block.meta_controller_container);
	PrimitiveContainer<U32>   positions(block.meta_positions_container);
	PrimitiveContainer<float> quality(block.meta_quality_container);
	PrimitiveContainer<BYTE>  refalt(block.meta_refalt_container);
	PrimitiveContainer<S32>   filterID(block.meta_filter_map_ids);
	PrimitiveContainer<S32>   formatID(block.meta_format_map_ids);
	PrimitiveContainer<S32>   infoID(block.meta_info_map_ids);

	for(U32 i = 0; i < this->size(); ++i){
		new( &this->__entries[i] ) value_type( );
	}

	if(contigs.size()){
		for(U32 i = 0; i < this->size(); ++i){
			if(contigs.isUniform()) this->__entries[i].contigID = contigs[0];
			else this->__entries[i].contigID = contigs[i];
		}
	}

	if(positions.size()){
		for(U32 i = 0; i < this->size(); ++i){
			if(positions.isUniform()) this->__entries[i].position = positions[0];
			else this->__entries[i].position = positions[i];
		}
	}

	if(controllers.size()){
		for(U32 i = 0; i < this->size(); ++i){
			if(controllers.isUniform()) this->__entries[i].controller = controllers[0];
			else this->__entries[i].controller = controllers[i];
		}
	}

	if(quality.size()){
		for(U32 i = 0; i < this->size(); ++i){
			if(quality.isUniform()) this->__entries[i].quality = quality[0];
			else this->__entries[i].quality = quality[i];
		}
	}

	if(filterID.size()){
		for(U32 i = 0; i < this->size(); ++i){
			if(filterID.isUniform()) this->__entries[i].filter_pattern_id = filterID[0];
			else this->__entries[i].filter_pattern_id = filterID[i];
		}
	}

	if(infoID.size()){
		for(U32 i = 0; i < this->size(); ++i){
			if(infoID.isUniform()) this->__entries[i].info_pattern_id = infoID[0];
			else this->__entries[i].info_pattern_id = infoID[i];
		}
	}

	if(formatID.size()){
		for(U32 i = 0; i < this->size(); ++i){
			if(formatID.isUniform()) this->__entries[i].format_pattern_id = formatID[0];
			else 	this->__entries[i].format_pattern_id = formatID[i];
		}
	}

	U32 refalt_position = 0;
	for(U32 i = 0; i < this->size(); ++i){
		if(this->__entries[i].controller.alleles_packed){
			// load from special packed
			// this is always diploid
			this->__entries[i].n_alleles = 2;
			this->__entries[i].alleles   = static_cast<value_type::allele_type*>(::operator new[](2*sizeof(value_type::allele_type)));

			// If data is <non_ref> or not
			if((refalt[refalt_position] & 15) != 5){
				const char ref = constants::REF_ALT_LOOKUP[refalt[refalt_position] & 15];
				new( &this->__entries[i].alleles[0] ) value_type::allele_type( ref );
			} else {
				new( &this->__entries[i].alleles[0] ) value_type::allele_type( "<NON_REF>" );
			}

			// If data is <non_ref> or not
			if(((refalt[refalt_position] >> 4) & 15) != 5){
				const char alt = constants::REF_ALT_LOOKUP[(refalt[refalt_position] >> 4) & 15];
				new( &this->__entries[i].alleles[1] ) value_type::allele_type( alt );
			} else {
				new( &this->__entries[i].alleles[1] ) value_type::allele_type( "<NON_REF>" );
			}
			++refalt_position;
		}
		// otherwise load from literal cold
		else {
			// number of alleles is parsed from the stride container
		}
	}
	assert(refalt_position == refalt.size());

	if(block.meta_alleles_container.buffer_data_uncompressed.size()){
		StrideContainer<U32> strides(block.meta_alleles_container);
		U32 offset = 0;
		U32 stride_offset = 0;
		for(U32 i = 0; i < this->size(); ++i){
			if(this->__entries[i].controller.alleles_packed == false){
				this->__entries[i].n_alleles = strides[stride_offset];
				this->__entries[i].alleles   = static_cast<value_type::allele_type*>(::operator new[](strides[stride_offset]*sizeof(value_type::allele_type)));

				for(U32 j = 0; j < strides[stride_offset]; ++j){
					const U16& l_string = *reinterpret_cast<const U16* const>(&block.meta_alleles_container.buffer_data_uncompressed[offset]);
					new( &this->__entries[i].alleles[j] ) value_type::allele_type( &block.meta_alleles_container.buffer_data_uncompressed[offset] );
					offset += sizeof(U16) + l_string;
				}
				++stride_offset;
			}
		}
		assert(offset == block.meta_alleles_container.getSizeUncompressed());
	}

	// Parse name
	if(block.meta_names_container.getSizeUncompressed()){
		StrideContainer<U32> strides(block.meta_names_container);
		U32 offset = 0;
		assert(strides.size() == this->size());
		for(U32 i = 0; i < this->size(); ++i){
			this->__entries[i].name = std::string(&block.meta_names_container.buffer_data_uncompressed.data()[offset], strides[i]);
			offset += strides[i];
		}
	}
}

}
}
