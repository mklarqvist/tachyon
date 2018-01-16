#ifndef CORE_ITERATOR_GENOTYPEITERATOR_H_
#define CORE_ITERATOR_GENOTYPEITERATOR_H_

#include "../BlockEntry.h"
#include "MetaIterator.h"

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
	typedef Core::StreamContainer container_type;
	typedef Core::BlockEntry block_type;
	typedef MetaIterator meta_iterator_type;
	typedef ContainerIterator container_iterator_type;
	typedef Core::MetaEntry meta_entry_type;

	typedef const U32 (self_type::*getFunctionType)(void) const;

public:
	GenotypeIterator(block_type& block) :
		current_position(0),
		width(0),
		status(YON_IT_START),
		container_rle(&block.gt_rle_container),
		container_simple(&block.gt_simple_container),
		container_meta(&block.gt_support_data_container),
		iterator_gt_meta(block.gt_support_data_container),
		iterator_meta(block.getMetaIterator()),
		getNumberObjectsFunction(nullptr)
	{
		switch(this->container_meta->header.controller.type){
		case(Core::YON_BYTE): this->getNumberObjectsFunction = &self_type::__getCurrentObjectLength<BYTE>; break;
		case(Core::YON_U16):  this->getNumberObjectsFunction = &self_type::__getCurrentObjectLength<U16>;  break;
		case(Core::YON_U32):  this->getNumberObjectsFunction = &self_type::__getCurrentObjectLength<U32>;  break;
		case(Core::YON_U64):  this->getNumberObjectsFunction = &self_type::__getCurrentObjectLength<U64>;  break;
		}
	}

	virtual ~GenotypeIterator(){
		delete this->iterator_meta;
	}

	// std::vector<some_abstract_genotype_object> toVector(void) const;
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
	inline const U32 getCurrentTargetStream(void) const{ return(this->iterator_gt_meta.getStrideIterator()->getCurrentStride()); }

	/**<
	 * Support functionality: returns the current MetaEntry from
	 * the meta entry iterator.
	 * @return Meta entry reference of current record
	 */
	inline meta_entry_type& getCurrentMeta(void) const{ return(this->iterator_meta->current()); }

	void operator++(void){
		++this->current_position;
		++(*this->iterator_meta);
		this->iterator_gt_meta.incrementSecondaryUsage(); // internals take care of this
	}

	void operator--(void);
	void operator[](const U32 p);
	void operator+=(const U32 p);
	void operator-=(const U32 p);

private:
	//virtual void toVCFString(std::ostream& stream) const =0;

	template <class T>
	inline const U32 __getCurrentObjectLength(void) const{
		const Iterator::ContainerIteratorType<T>& it = *reinterpret_cast< const Iterator::ContainerIteratorType<T>* const >(iterator_gt_meta.getDataIterator());
		return(it.current());
	}

public:
	U32 current_position;
	U32 width;
	TACHYON_ITERATOR_STATUS status;

	// RLE iterator
	// simple iterator
	const container_type*     container_rle;
	const container_type*     container_simple;
	const container_type*     container_meta;
	container_iterator_type   iterator_gt_meta;  // iterator over n_runs and target stream
	meta_iterator_type*       iterator_meta;
	// function pointer to getNumberObjects
	getFunctionType           getNumberObjectsFunction;
};


}
}


#endif /* CORE_ITERATOR_GENOTYPEITERATOR_H_ */
