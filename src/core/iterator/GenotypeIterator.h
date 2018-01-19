#ifndef CORE_ITERATOR_GENOTYPEITERATOR_H_
#define CORE_ITERATOR_GENOTYPEITERATOR_H_

#include "../Block.h"
#include "MetaIterator.h"
#include "IntegerIterator.h"
#include "../base/GTDiploidObject.h"

namespace Tachyon{
namespace Iterator{

/**
 * We have several different GT representations
 *
 * Case diploid, bi-allelic, no NA, no missing, and run-length encoded
 * Case diploid, bi-allelic, no NA, has missing, and run-length encoded
 * Case diploid, n-allelic or has NA
 * Case diploid, n-allelic or has NA run-length encoded
 * Case m-ploid, n-allelic
 */
class GenotypeIterator{
private:
	typedef GenotypeIterator self_type;
	typedef Core::Container container_type;
	typedef Core::Block block_type;
	typedef MetaIterator meta_iterator_type;
	typedef ContainerIterator container_iterator_type;
	typedef Core::MetaEntry meta_entry_type;
	typedef IntegerIterator integer_iterator_type;
	typedef Core::GTDiploidObject gt_diploid_object_type;

	typedef const U32 (self_type::*getFunctionType)(void) const;
	typedef const U64 (self_type::*getGTFunctionType)(void) const;

public:
	GenotypeIterator(block_type& block) :
		current_position(0),
		width(0),
		status(YON_IT_START),
		__diploid_object_type(2),
		container_rle(&block.gt_rle_container),
		container_simple(&block.gt_simple_container),
		container_meta(&block.gt_support_data_container),
		iterator_gt_meta(block.gt_support_data_container),
		iterator_gt_rle(block.gt_rle_container),
		iterator_gt_simple(block.gt_simple_container),
		iterator_meta(block.getMetaIterator()),
		getNumberObjectsFunction(nullptr)
	{
		switch(this->container_meta->header.controller.type){
		case(Core::YON_GT_BYTE): this->getNumberObjectsFunction = &self_type::__getCurrentObjectLength<BYTE>; break;
		case(Core::YON_GT_U16):  this->getNumberObjectsFunction = &self_type::__getCurrentObjectLength<U16>;  break;
		case(Core::YON_GT_U32):  this->getNumberObjectsFunction = &self_type::__getCurrentObjectLength<U32>;  break;
		case(Core::YON_GT_U64):  this->getNumberObjectsFunction = &self_type::__getCurrentObjectLength<U64>;  break;
		}
		this->__updateGTPrimitive();
	}

	virtual ~GenotypeIterator(){
		delete this->iterator_meta;
	}

	inline const meta_iterator_type* getIteratorMeta(void){ return(this->iterator_meta); }

	/**<
	 * Reset and reuse
	 */
	inline void reset(void){
		this->current_position = 0;
		this->iterator_meta->reset();
		this->iterator_gt_meta.reset();
		this->status = YON_IT_START;
	}

	/**<
	 * Returns the GT object length: this is either
	 * the number of runs in an object or simple the
	 * value 1 for non-run-length encoded objects
	 * @return
	 */
	inline const U32 getCurrentObjectLength(void) const{ return((this->*getNumberObjectsFunction)()); }

	/**<
	 * Returns the target stream for housing the genotype data
	 * for the current variant
	 * @return
	 */
	inline const U32 getCurrentTargetStream(void) const{
		if(this->container_meta->header.controller.mixedStride == false) return(this->container_meta->header.stride);
		return(this->iterator_gt_meta.getStrideIterator()->getCurrentStride());
	}

	/**<
	 * Support functionality: returns the current MetaEntry from
	 * the meta entry iterator.
	 * @return Meta entry reference of current record
	 */
	inline meta_entry_type& getCurrentMeta(void) const{ return(this->iterator_meta->current()); }

	/**<
	 * Returns the current GT object. If there is more than one
	 * object at the current variant site the function
	 * getCurrentGTObjects() has to be invoked to return
	 * the pointer to the start of those objects
	 * @return
	 */
	inline const gt_diploid_object_type& getCurrentGTObject(void){
		if(this->getCurrentTargetStream() == 1) this->__diploid_object_type(this->iterator_gt_rle.current(), this->iterator_meta->current());
		else this->__diploid_object_type(this->iterator_gt_simple.current(), this->iterator_meta->current());
		return(this->__diploid_object_type);
	}

	/**<
	 * Increment the GT iterator only! Not the other structures
	 */
	void incrementGT(void){
		if(this->getCurrentTargetStream() == 1) ++this->iterator_gt_rle;
		else ++this->iterator_gt_simple;
	}

	/**<
	 * General incrementor over internal objects
	 */
	void operator++(void){
		++this->current_position;
		++(*this->iterator_meta);
		this->iterator_gt_meta.incrementSecondaryUsage(1); // internals take care of this

		this->__updateGTPrimitive();
	}

private:
	//virtual void toVCFString(std::ostream& stream) const =0;

	bool __updateGTPrimitive(void){
		// Next meta
		const meta_entry_type& meta = this->iterator_meta->current();
		const BYTE primitive        = meta.hot.controller.gt_primtive_type;
		const U32 target_stream     = this->getCurrentTargetStream();

		if(primitive == Core::YON_GT_BYTE){
			if(target_stream == 1) this->iterator_gt_rle.setType(Core::YON_TYPE_8B);
			else this->iterator_gt_simple.setType(Core::YON_TYPE_8B);
		} else if(primitive == Core::YON_GT_U16){
			if(target_stream == 1) this->iterator_gt_rle.setType(Core::YON_TYPE_16B);
			else this->iterator_gt_simple.setType(Core::YON_TYPE_16B);
		} else if(primitive == Core::YON_GT_U32){
			if(target_stream == 1) this->iterator_gt_rle.setType(Core::YON_TYPE_32B);
			else this->iterator_gt_simple.setType(Core::YON_TYPE_32B);
		} else if(primitive == Core::YON_GT_U64){
			if(target_stream == 1) this->iterator_gt_rle.setType(Core::YON_TYPE_64B);
			else this->iterator_gt_simple.setType(Core::YON_TYPE_64B);
		} else {
			std::cerr << "impossible primitive: " << (int)primitive << std::endl;
			exit(1);
		}
		return(true);
	}

	template <class T>
	inline const U32 __getCurrentObjectLength(void) const{ return(iterator_gt_meta.getDataIterator()->current<T>()); }

public:
	U32 current_position;
	U32 width;
	TACHYON_ITERATOR_STATUS status;

	gt_diploid_object_type __diploid_object_type;

	// RLE iterator
	// simple iterator
	const container_type*     container_rle;
	const container_type*     container_simple;
	const container_type*     container_meta;
	container_iterator_type   iterator_gt_meta;  // iterator over n_runs and target stream
	integer_iterator_type     iterator_gt_rle;
	integer_iterator_type     iterator_gt_simple;
	meta_iterator_type*       iterator_meta;
	// function pointer to getNumberObjects
	getFunctionType           getNumberObjectsFunction;
};


}
}


#endif /* CORE_ITERATOR_GENOTYPEITERATOR_H_ */
