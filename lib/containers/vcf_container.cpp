#include "vcf_container.h"

namespace tachyon {
namespace containers {

VcfContainer::VcfContainer(void) :
	n_carry_over_(0),
	n_entries_(0),
	n_capacity_(500),
	entries_(new pointer[500])
{
	for (size_type i = 0; i < this->capacity(); ++i)
		this->entries_[i] = nullptr;
}

VcfContainer::VcfContainer(const size_type& start_capacity) :
	n_carry_over_(0),
	n_entries_(0),
	n_capacity_(start_capacity),
	entries_(new pointer[start_capacity])
{
	for (size_type i = 0; i < this->capacity(); ++i)
		this->entries_[i] = nullptr;
}

VcfContainer::~VcfContainer() {
	if (this->entries_ != nullptr) {
		for (std::size_t i = 0; i < this->n_entries_; ++i)
			bcf_destroy(this->entries_[i]);

		::operator delete[](static_cast<void*>(this->entries_));
	}
}

VcfContainer::VcfContainer(self_type&& other) noexcept :
	n_carry_over_(0),
	n_entries_(other.n_entries_),
	n_capacity_(other.n_capacity_),
	entries_(other.entries_)
{
	other.entries_ = new pointer[other.n_capacity_];
	for (size_type i = 0; i < other.capacity(); ++i)
		other.entries_[i] = nullptr;

	if (other.n_carry_over_) {
		other.entries_[0] = this->at(this->size()-1);
		assert(this->at(this->size()-1) != nullptr);
		this->entries_[this->size()-1] = nullptr;
		other.n_carry_over_ = 0;
		other.n_entries_    = 1;
		--this->n_entries_;
	} else {
		other.n_entries_ = 0;
		other.n_carry_over_ = 0;
	}
}

VcfContainer& VcfContainer::operator=(self_type&& other) noexcept
{
	if (this->entries_ != nullptr) {
		for (std::size_t i = 0; i < this->n_entries_; ++i)
			bcf_destroy(this->entries_[i]);

		::operator delete[](static_cast<void*>(this->entries_));
	}

	this->n_carry_over_ = 0;
	this->n_capacity_   = other.n_capacity_;
	this->entries_      = other.entries_;
	this->n_entries_    = other.n_entries_;

	other.entries_ = new pointer[other.n_capacity_];
	for (size_type i = 0; i < other.capacity(); ++i)
		other.entries_[i] = nullptr;

	if (other.n_carry_over_) {
		other.entries_[0] = this->at(this->size()-1);
		assert(this->at(this->size()-1) != nullptr);
		this->entries_[this->size()-1] = nullptr;
		other.n_carry_over_ = 0;
		other.n_entries_    = 1;
		--this->n_entries_;
	} else {
		other.n_entries_ = 0;
		other.n_carry_over_ = 0;
	}
	return(*this);
}

void VcfContainer::resize(const size_t new_size) {
	if (new_size < this->capacity()) {
		for (size_t i = new_size; i < this->n_entries_; ++i)
			bcf_destroy(this->entries_[i]);

		if (this->n_entries_ >= new_size) this->n_entries_ = new_size;
		return;
	}

	pointer* temp = new pointer[new_size];
	for (size_t i = 0; i < new_size; ++i) temp[i] = nullptr;
	for (size_t i = 0; i < this->size(); ++i) temp[i] = this->entries_[i];

	delete [] this->entries_;
	this->entries_    = temp;
	this->n_capacity_ = new_size;
}

bool VcfContainer::GetVariants(const int32_t n_variants,
                               const int64_t n_bases,
                               std::unique_ptr<io::VcfReader>& reader,
							   const int unpack_level)
{
	// Make sure enough records are available for
	// overloading. If this is not the case we
	// resize the container to the target size
	// plus some padding.
	if (this->size() + n_variants + this->n_carry_over_ >= this->capacity())
		this->resize(this->size() + n_variants + this->n_carry_over_ + 64);


	// Look at the back of the vector and check if
	// that entry is a valid pointer to a htslib
	// bcf1_t structure. If it is not then initialize
	// a new bcf1_t structure. Next, read a new
	// vcf records from the file stream and overload
	// the bcf1_t record.
	VcfContainer::pointer bcf1_ = this->end();
	if (bcf1_ == nullptr) bcf1_  = bcf_init();
	if (reader->next(bcf1_, unpack_level) == false)
		return false;

	// Push the overloaded bcf1_t structure to this
	// container in the back.
	*this += bcf1_;

	// Track the current position and contig.
	int64_t first_pos    = bcf1_->pos;
	int32_t first_contig = bcf1_->rid;
	if (this->size() != 1) {
		first_pos    = this->entries_[0]->pos;
		first_contig = this->entries_[0]->rid;
	}

	if (bcf1_->pos - first_pos > n_bases || first_contig != bcf1_->rid) {
		this->n_carry_over_ = 1;
		return(this->size() - 1);
	}

	for (int32_t i = 1; i < n_variants; ++i) {
		bcf1_ = this->end();
		if (bcf1_ == nullptr) bcf1_  = bcf_init();

		if (reader->next(bcf1_, unpack_level) == false)
			return(this->size());

		*this += bcf1_;

		if (bcf1_->pos - first_pos > n_bases || first_contig != bcf1_->rid) {
			this->n_carry_over_ = 1;
			return(this->size() - 1);
		}
	}

	return(this->size());
}

GenotypeSummary VcfContainer::GetGenotypeSummary(const uint32_t position, const uint64_t& n_samples) const{
	GenotypeSummary g;

	// If there are no FORMAT fields there cannot exist any
	// GT data.
	if (this->at(position)->n_fmt == 0)
		return(g);

	// Iterate through the allowed primitive types for genotypes to collect summary
	// statistics for genotypes at this loci. Information collected includes the
	// base ploidy, if there's any mixed phasing, the number of missing genotypes, and
	// the number of samples that has a special end-of-vector encoding.
	// Only the signed primitives int8_t, int16_t, and int32_t are valid for genotypes.
	switch(this->at(position)->d.fmt[0].type) {
	case(BCF_BT_INT8):  g.Evaluate<int8_t> (n_samples, this->at(position)->d.fmt[0]); break;
	case(BCF_BT_INT16): g.Evaluate<int16_t>(n_samples, this->at(position)->d.fmt[0]); break;
	case(BCF_BT_INT32): g.Evaluate<int32_t>(n_samples, this->at(position)->d.fmt[0]); break;
	case(BCF_BT_NULL):
	case(BCF_BT_FLOAT):
	case(BCF_BT_CHAR):
	default:
		std::cerr << "Illegal genotype primtive type: " << this->at(position)->d.fmt[0].type << std::endl;
	}

	return(g);
}

void VcfContainer::clear(void) {
	uint32_t start_pos = 0;
	if (this->n_carry_over_) {
		assert(this->size() != 0);
		std::swap(this->entries_[this->size() - 1], this->entries_[0]);
		start_pos = 1;
		this->n_carry_over_ = 0;
	}

	for (uint32_t i = start_pos; i < this->size(); ++i)
		bcf_clear(this->entries_[i]);

	this->n_entries_ = start_pos;
}

}
}
