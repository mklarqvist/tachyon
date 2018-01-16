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

public:
	GenotypeIterator(block_type& block) :
		container_rle(&block.gt_rle_container),
		container_simple(&block.gt_simple_container),
		container_meta(&block.gt_support_data_container),
		iterator_gt_meta(block.gt_support_data_container),
		iterator_meta(block.getMetaIterator())
{}

	virtual ~GenotypeIterator(){}

	// std::vector<some_abstract_genotype_object> toVector(void) const;
	inline const meta_iterator_type* getIteratorMeta(void){ return(this->iterator_meta); }

private:
	//virtual void toVCFString(std::ostream& stream) const =0;

	/**<
	 * Returns the target stream for holding the current
	 * @return
	 */
	const U32& getCurrentGTTargetStream(void);

	/**<
	 * Returns the GT object length: this is either
	 * the number of runs in an object or simple the
	 * value 1 for non-run-length encoded objects
	 * @return
	 */
	const U32& getCurrentGTLength(void);

public:
	// RLE iterator
	// simple iterator
	const container_type*     container_rle;
	const container_type*     container_simple;
	const container_type*     container_meta;
	container_iterator_type   iterator_gt_meta;  // iterator over n_runs and target stream
	meta_iterator_type* iterator_meta;
};


}
}


#endif /* CORE_ITERATOR_GENOTYPEITERATOR_H_ */
