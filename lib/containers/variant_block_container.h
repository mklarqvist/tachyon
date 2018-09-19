#ifndef CONTAINERS_VARIANT_BLOCK_CONTAINER_H_
#define CONTAINERS_VARIANT_BLOCK_CONTAINER_H_

#include "variant_block.h"
#include "containers/genotype_container.h"
#include "containers/interval_container.h"
#include "core/variant_record.h"

#include "algorithm/compression/compression_manager.h"
#include "algorithm/encryption/encryption_decorator.h"

namespace tachyon {
namespace containers {

/**<
 * Decouples the low-level object `VariantBlock` that holds all internal data and
 * strides. This container should preferably be the primary object the end-user
 * interacts with.
 */
class VariantBlockContainer {
private:
	typedef VariantBlockContainer self_type;
    typedef DataContainerHeader   data_header_type;
    typedef VariantBlock          block_type;
    typedef VariantBlockHeader    block_header_type;
    typedef VariantBlockFooter    block_footer_type;
    typedef VariantHeader         global_header_type;
	typedef containers::VariantBlock             block_entry_type;
	typedef containers::GenotypeContainer        gt_container_type;
	typedef containers::IntervalContainer        interval_container_type;
	typedef DataBlockSettings                    block_settings_type;

	typedef algorithm::CompressionManager  compression_manager_type;
	typedef EncryptionDecorator            encryption_manager_type;

	typedef std::unordered_map<int, int>    map_type;

public:
	VariantBlockContainer(void);
	VariantBlockContainer(const global_header_type& header);
	VariantBlockContainer(const self_type& other);
	VariantBlockContainer(self_type&& other) noexcept;
	VariantBlockContainer& operator=(const self_type& other);
	VariantBlockContainer& operator=(self_type&& other) noexcept;
	~VariantBlockContainer(void);

	// Adds global header pointer to this object.
	self_type& operator<<(const global_header_type& header){
		this->header_ = &header;
		return(*this);
	}

	inline void reset(void){ this->block_.clear(); }

	/**< @brief Reads one or more separate digital objects from disk
	 * Primary function for reading partial data from disk. Data
	 * read in this way is not checked for integrity until later.
	 * @param stream   Input stream
	 * @param settings Settings record describing reading parameters
	 * @return         Returns FALSE if there was a problem, TRUE otherwise
	 */
	bool ReadBlock(std::ifstream& stream, block_settings_type& settings);

	// Accessors
	inline block_type& GetBlock(void){ return(this->block_); }
	inline const block_type& GetBlock(void) const{ return(this->block_); }

	// Checkers
	inline bool AnyEncrypted(void) const{ return(this->block_.header.controller.any_encrypted); }
	inline bool HasGenotypes(void) const{ return(this->block_.header.controller.has_gt); }
	inline bool HasPermutedGenotypes(void) const{ return(this->block_.header.controller.has_gt_permuted); }

	inline void AllocateGenotypeMemory(void){
		delete [] this->gt_exp;
		this->gt_exp = new yon_gt_rcd[this->header_->GetNumberSamples()];
	}
	inline yon_gt_rcd* GetAllocatedGenotypeMemory(void){ return(this->gt_exp); }
	inline yon_gt_rcd* GetAllocatedGenotypeMemory(void) const{ return(this->gt_exp); }

private:
	block_type                block_;
	compression_manager_type  compression_manager;
	encryption_manager_type   encryption_manager;
	const global_header_type* header_;
	// External memory allocation for linear use of lazy-evaluated
	// expansion of genotype records. This is critical when the sample
	// numbers are becoming large as allocating/deallocating hundreds
	// of thousands of pointers for every variant is very time consuming.
	yon_gt_rcd* gt_exp;
};

}
}


#endif /* CONTAINERS_VARIANT_BLOCK_CONTAINER_H_ */
