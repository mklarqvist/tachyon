#include "variant_block_container.h"

namespace tachyon{
namespace containers{

VariantBlockContainer::VariantBlockContainer() :
	header_(nullptr),
	gt_exp(nullptr)
{

}

VariantBlockContainer::VariantBlockContainer(const global_header_type& header) :
	header_(&header),
	gt_exp(nullptr)
{

}

VariantBlockContainer::VariantBlockContainer(const self_type& other) :
	block_(other.block_),
	compression_manager(other.compression_manager),
	encryption_manager(other.encryption_manager),
	header_(other.header_),
	gt_exp(nullptr)
{
	if(other.gt_exp != nullptr){
		this->gt_exp = new yon_gt_rcd*[this->header_->GetNumberSamples()];
		for(uint32_t i = 0; i < this->header_->GetNumberSamples(); ++i)
			this->gt_exp[i] = other.gt_exp[i];
	}
}

VariantBlockContainer::VariantBlockContainer(self_type&& other) noexcept :
	block_(std::move(other.block_)),
	compression_manager(std::move(other.compression_manager)),
	encryption_manager(std::move(other.encryption_manager)),
	header_(nullptr),
	gt_exp(nullptr)
{
	std::swap(this->gt_exp,  other.gt_exp);
	std::swap(this->header_, other.header_);
}

VariantBlockContainer& VariantBlockContainer::operator=(const self_type& other){
	delete [] this->gt_exp;
	return *this = VariantBlockContainer(other);
}

VariantBlockContainer& VariantBlockContainer::operator=(self_type&& other) noexcept{
	if(this == &other){
		// precautions against self-moves
		return *this;
	}

	this->block_  = std::move(other.block_);
	this->compression_manager = std::move(other.compression_manager);
	this->encryption_manager = std::move(other.encryption_manager);
	this->header_ = nullptr;
	delete [] this->gt_exp; this->gt_exp  = nullptr;
	std::swap(this->gt_exp, other.gt_exp);
	std::swap(this->header_, other.header_);
	return *this;
}

VariantBlockContainer::~VariantBlockContainer(void)
{
	// Do not delete header pointer as this object
	// never owns that data.
	delete [] this->gt_exp;
}

bool VariantBlockContainer::ReadBlock(std::ifstream& stream, block_settings_type& settings){
	return(this->block_.read(stream, settings, *this->header_));
}

}
}
