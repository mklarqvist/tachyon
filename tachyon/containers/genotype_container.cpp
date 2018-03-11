#include "genotype_container.h"
#include "stride_container.h"

namespace tachyon{
namespace containers{

GenotypeContainer::GenotypeContainer(const block_type& block) :
	n_entries(0),
	__meta_container(block),
	__iterators(nullptr)
{
	this->n_entries   = this->__meta_container.size();
	this->__iterators = static_cast<pointer>(::operator new[](this->size() * sizeof(value_type)));

	if(this->n_entries == 0)
		return;

	// Aliases
	const char* const data_rle    = block.gt_rle_container.buffer_data_uncompressed.data();
	const char* const data_simple = block.gt_simple_container.buffer_data_uncompressed.data();

	if(block.gt_support_data_container.buffer_data_uncompressed.size() == 0){
		std::cerr << utility::timestamp("ERROR","GT") << "Has no genotype support data!" << std::endl;
		exit(1);
	}

	// data (0: rle, 1: simple), strides (n_objects)
	getNativeFuncDef getObjects = nullptr;
	switch(block.gt_support_data_container.getDataPrimitiveType()){
	case(YON_TYPE_8B):  getObjects = &self_type::getNative<BYTE>; break;
	case(YON_TYPE_16B): getObjects = &self_type::getNative<U16>; break;
	case(YON_TYPE_32B): getObjects = &self_type::getNative<U32>; break;
	case(YON_TYPE_64B): getObjects = &self_type::getNative<U64>; break;
	default: std::cerr << "illegal type" << std::endl; return;
	}

	// Meta data is uniform
	if(block.gt_support_data_container.header.data_header.isUniform()){
		U32 current_offset_rle    = 0;
		U32 current_offset_simple = 0;
		const U32 target = block.gt_support_data_container.header.data_header.stride;
		const U32 n_objects = (this->*getObjects)(block.gt_support_data_container.buffer_data_uncompressed, 0);
		assert(block.gt_support_data_container.header.data_header.stride > 0);

		for(U32 i = 0; i < this->n_entries; ++i){
			if(target == 1){
				if(this->__meta_container[i].getGTPrimitiveWidth() == 1)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<BYTE>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 2)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U16>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 4)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U32>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 8)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U64>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );

				current_offset_rle += n_objects * this->__meta_container[i].getGTPrimitiveWidth();
			} else if(target == 2){
				if(this->__meta_container[i].getGTPrimitiveWidth() == 1)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<BYTE>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 2)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U16>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 4)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U32>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 8)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U64>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );

				current_offset_simple += n_objects * this->__meta_container[i].getGTPrimitiveWidth();
			} else {
				std::cerr << utility::timestamp("ERROR") << "Illegal GT specification!" << std::endl;
				exit(1);
			}
		}
		return;
	}

	// Start non-uniform (standard ctor)

	if(block.gt_support_data_container.header.data_header.hasMixedStride()){
		U32 current_offset_rle    = 0;
		U32 current_offset_simple = 0;
		StrideContainer<U32> strides(block.gt_support_data_container);

		for(U32 i = 0; i < this->n_entries; ++i){
			const U32 n_objects = (this->*getObjects)(block.gt_support_data_container.buffer_data_uncompressed, i);

			if(strides[i] == 1){
				if(this->__meta_container[i].getGTPrimitiveWidth() == 1)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<BYTE>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 2)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U16>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 4)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U32>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 8)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U64>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );

				current_offset_rle += n_objects * this->__meta_container[i].getGTPrimitiveWidth();
			} else if(strides[i] == 2){
				if(this->__meta_container[i].getGTPrimitiveWidth() == 1)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<BYTE>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 2)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U16>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 4)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U32>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 8)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U64>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );

				current_offset_simple += n_objects * this->__meta_container[i].getGTPrimitiveWidth();
			} else {
				std::cerr << utility::timestamp("ERROR") << "Illegal GT specification!" << std::endl;
				exit(1);
			}
		}

		assert(current_offset_rle == block.gt_rle_container.buffer_data_uncompressed.size());
		assert(current_offset_simple == block.gt_simple_container.buffer_data_uncompressed.size());
	}
	else { // No mixed stride
		U32 current_offset_rle    = 0;
		U32 current_offset_simple = 0;
		const U32 target = block.gt_support_data_container.header.data_header.stride;

		for(U32 i = 0; i < this->n_entries; ++i){
			const U32 n_objects = (this->*getObjects)(block.gt_support_data_container.buffer_data_uncompressed, i);

			if(target == 1){
				if(this->__meta_container[i].getGTPrimitiveWidth() == 1)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<BYTE>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 2)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U16>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 4)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U32>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 8)
					new( &this->__iterators[i] ) GenotypeContainerDiploidRLE<U64>( &data_rle[current_offset_rle], n_objects, this->__meta_container[i] );

				current_offset_rle += n_objects * this->__meta_container[i].getGTPrimitiveWidth();
			} else if(target == 2){
				if(this->__meta_container[i].getGTPrimitiveWidth() == 1)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<BYTE>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 2)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U16>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 4)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U32>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );
				else if(this->__meta_container[i].getGTPrimitiveWidth() == 8)
					new( &this->__iterators[i] ) GenotypeContainerDiploidSimple<U64>( &data_simple[current_offset_simple], n_objects, this->__meta_container[i] );

				current_offset_simple += n_objects * this->__meta_container[i].getGTPrimitiveWidth();
			} else {
				std::cerr << utility::timestamp("ERROR") << "Illegal GT specification!" << std::endl;
				exit(1);
			}
		}

		assert(current_offset_rle == block.gt_rle_container.buffer_data_uncompressed.size());
		assert(current_offset_simple == block.gt_simple_container.buffer_data_uncompressed.size());
	}
}

GenotypeContainer::~GenotypeContainer(){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		(this->__iterators + i)->~GenotypeContainerInterface();

	::operator delete[](static_cast<void*>(this->__iterators));
}

}
}
