#ifndef CONTAINERS_VCF_CONTAINER_H_
#define CONTAINERS_VCF_CONTAINER_H_

#include "components/generic_iterator.h"
#include "io/vcf_utils.h"

namespace tachyon{
namespace containers{

class VcfContainer{
public:
	typedef VcfContainer       self_type;
	typedef bcf1_t             value_type;
	typedef value_type&        reference;
	typedef const value_type&  const_reference;
	typedef value_type*        pointer;
	typedef const value_type*  const_pointer;
	typedef std::ptrdiff_t     difference_type;
	typedef std::size_t        size_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
	VcfContainer(void);
	VcfContainer(const size_type& start_capacity);
	VcfContainer(const VcfContainer& other) = delete; // Disallow copy ctor

	VcfContainer& operator=(self_type&& other) noexcept
	{
		if(this->entries_ != nullptr){
			for(std::size_t i = 0; i < this->n_entries_; ++i)
				bcf_destroy(this->entries_[i]);

			::operator delete[](static_cast<void*>(this->entries_));
		}

		this->n_carry_over_ = 0;
		this->n_capacity_   = other.n_capacity_;
		this->entries_      = other.entries_;
		this->n_entries_    = other.n_entries_;

		other.entries_ = new pointer[other.n_capacity_];
		for(size_type i = 0; i < other.capacity(); ++i)
			other.entries_[i] = nullptr;

		if(other.n_carry_over_){
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

	~VcfContainer();

	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline size_type sizeWithoutCarryOver(void) const{ return(this->n_entries_ - this->n_carry_over_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }
	inline pointer front(void){ return(this->entries_[0]); }
	inline const_pointer front(void) const{ return(this->entries_[0]); }
	inline pointer back(void){ return(this->entries_[this->size() == 0 ? 0 : this->size() - 1 - this->n_carry_over_]); }
	inline const_pointer back(void) const{ return(this->entries_[this->size() == 0 ? 0 : this->size() - 1 - this->n_carry_over_]); }

	inline void operator+=(const pointer entry){ this->entries_[this->n_entries_++] = entry; }
	inline pointer operator[](const uint32_t position){ return(this->entries_[position]); }
	inline const_pointer operator[](const uint32_t position) const{ return(this->entries_[position]); }
	inline pointer at(const uint32_t position){ return(this->entries_[position]); }
	inline const_pointer at(const uint32_t position) const{ return(this->entries_[position]); }

	inline pointer end(void){ return(this->entries_[this->n_entries_]); }
	inline const_pointer end(void) const{ return(this->entries_[this->n_entries_]); }

	void resize(const size_t new_size);
	bool GetVariants(const int32_t n_variants, const int64_t n_bases, std::unique_ptr<io::VcfReader>& reader);

	// Calculate genotype summary statistics from a lazy evaluated bcf1_t struct.
	// Warning: this function does NOT check if the FORMAT field GT exists either
	// in the header or in the structure itself. The assumption is that it does
	// exist and according to the Bcf specification has to be the first FORMAT
	// field set.
	io::VcfGenotypeSummary GetGenotypeSummary(const uint32_t position, const uint64_t& n_samples) const;
	void clear(void);

public:
	uint32_t  n_carry_over_;
	size_type n_entries_;
	size_type n_capacity_;
	pointer*  entries_;
};

}
}



#endif /* CONTAINERS_VCF_CONTAINER_H_ */
