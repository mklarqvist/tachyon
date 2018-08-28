#ifndef CORE_VARIANT_IMPORTER_CONTAINER_STATS_H_
#define CORE_VARIANT_IMPORTER_CONTAINER_STATS_H_

#include "containers/data_container.h"

namespace tachyon{
namespace support{

struct VariantImporterStatsObject{
	typedef VariantImporterStatsObject self_type;
	typedef containers::DataContainer data_container_type;

	VariantImporterStatsObject(void) : cost_uncompressed(0), cost_compressed(0){}
	~VariantImporterStatsObject(){}

	self_type& operator+=(const self_type& other){
		this->cost_compressed += other.cost_compressed;
		this->cost_uncompressed += other.cost_uncompressed;
		return(*this);
	}

	self_type& operator+=(const data_container_type& container){
		this->cost_uncompressed += container.GetObjectSizeUncompressed();
		this->cost_compressed   += container.GetObjectSize();
		return(*this);
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.cost_compressed << '\t' << entry.cost_uncompressed << '\t' <<
				(entry.cost_compressed ? (double)entry.cost_uncompressed/entry.cost_compressed : 0);
		return(stream);
	}

	uint64_t cost_uncompressed;
	uint64_t cost_compressed;
};

class VariantImporterContainerStats{
public:
	typedef VariantImporterContainerStats self_type;
	typedef VariantImporterStatsObject    value_type;
    typedef std::size_t                   size_type;
    typedef value_type&                   reference;
    typedef const value_type&             const_reference;
    typedef value_type*                   pointer;
    typedef const value_type*             const_pointer;

public:
    VariantImporterContainerStats() :
    	n_entries_(0),
		n_capacity_(0),
		entries_(nullptr)
    {
    }

    VariantImporterContainerStats(const size_type start_capacity) :
    	n_entries_(0),
		n_capacity_(start_capacity),
		entries_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
    {

    }

    ~VariantImporterContainerStats(){
    	for(std::size_t i = 0; i < this->size(); ++i)
			(this->entries_ + i)->~VariantImporterStatsObject();

		::operator delete[](static_cast<void*>(this->entries_));
    }

    VariantImporterContainerStats(const self_type& other) :
    	n_entries_(other.n_entries_),
		n_capacity_(other.n_capacity_),
		entries_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
    {
    	for(uint32_t i = 0; i < this->size(); ++i)
    		new( &this->entries_[i] ) value_type( other.at(i) );
    }

    VariantImporterContainerStats(self_type&& other) noexcept :
    	n_entries_(other.n_entries_),
		n_capacity_(other.n_capacity_),
		entries_(nullptr)
    {
    	std::swap(this->entries_, other.entries_);
    }

    VariantImporterContainerStats& operator=(self_type&& other) noexcept {
    	this->n_entries_ = other.n_entries_;
    	this->n_capacity_ = other.n_capacity_;
    	// Clear local entries if any.
    	for(std::size_t i = 0; i < this->size(); ++i)
			(this->entries_ + i)->~VariantImporterStatsObject();

		::operator delete[](static_cast<void*>(this->entries_));
		this->entries_ = nullptr;

    	std::swap(this->entries_, other.entries_);
    	return(*this);
	}

    VariantImporterContainerStats& operator=(const self_type& other) {
		this->n_entries_ = other.n_entries_;
		this->n_capacity_ = other.n_capacity_;
		// Clear local entries if any.
		for(std::size_t i = 0; i < this->size(); ++i)
			(this->entries_ + i)->~VariantImporterStatsObject();

		::operator delete[](static_cast<void*>(this->entries_));
		this->entries_ = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));

		for(uint32_t i = 0; i < this->size(); ++i)
			new( &this->entries_[i] ) value_type( other.at(i) );

		return(*this);
	}

	// Element access
	inline reference at(const size_type& position){ return(this->entries_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->entries_[position]); }
	inline reference operator[](const size_type& position){ return(this->entries_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->entries_[position]); }
	inline pointer data(void){ return(this->entries_); }
	inline const_pointer data(void) const{ return(this->entries_); }
	inline reference front(void){ return(this->entries_[0]); }
	inline const_reference front(void) const{ return(this->entries_[0]); }
	inline reference back(void){ return(this->entries_[this->n_entries_ - 1]); }
	inline const_reference back(void) const{ return(this->entries_[this->n_entries_ - 1]); }

	// Capacity
	inline bool empty(void) const{ return(this->n_entries_ == 0); }
	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	self_type& operator+=(const self_type& other){
		while(other.size() > this->size()) this->resize();

		for(uint32_t i = 0; i < other.size(); ++i)
			this->at(i) += other[i];

		return(*this);
	}

	//
	inline self_type& operator+=(const const_reference entry){
		if(this->size() + 1 == this->n_capacity_)
			this->resize();


		this->entries_[this->n_entries_++] = entry;
		return(*this);
	}
	inline self_type& add(const const_reference entry){ return(*this += entry); }

	void resize(const uint32_t new_size){
		if(new_size == this->capacity()) return;
		if(new_size < this->capacity()){
			if(new_size < this->size()){
				for(uint32_t i = new_size; i < this->size(); ++i)
					(this->entries_ + i)->~VariantImporterStatsObject();
				this->n_entries_ = new_size;
			}
			return;
		}

		pointer temp = this->entries_;

		this->n_capacity_ = new_size;
		this->entries_    = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));

		// Lift over values from old addresses
		for(uint32_t i = 0; i < this->size(); ++i)
			new( &this->entries_[i] ) value_type( temp[i] );

		if(temp != nullptr){
			for(std::size_t i = 0; i < this->size(); ++i)
				(temp + i)->~VariantImporterStatsObject();

			::operator delete[](static_cast<void*>(temp));
		}
	}

	inline void resize(void){ return(this->resize(this->capacity()*2)); }

	void Allocate(const uint32_t entries){
		// Delete previous entries.
		if(this->entries_ != nullptr){
			for(std::size_t i = 0; i < this->size(); ++i)
				(this->entries_ + i)->~VariantImporterStatsObject();
			::operator delete[](static_cast<void*>(this->entries_));

			this->entries_ = nullptr;
		}

		if(entries > this->capacity()) this->resize(entries);

		this->n_entries_ = entries;
		for(uint32_t i = 0; i < this->size(); ++i)
			new( &this->entries_[i] ) value_type( );
	}

private:
	size_type n_entries_;
	size_type n_capacity_;
	pointer   entries_;
};

}
}

#endif /* CORE_VARIANT_IMPORTER_CONTAINER_STATS_H_ */
