#ifndef CONTAINERS_PRIMITIVE_GROUP_CONTAINER_STRING_H_
#define CONTAINERS_PRIMITIVE_GROUP_CONTAINER_STRING_H_

#include "components/generic_iterator.h"
#include "primitive_group_container.h"

namespace tachyon{
namespace containers{

template <>
class PrimitiveGroupContainer<std::string> : public PrimitiveGroupContainerInterface{
private:
    typedef PrimitiveGroupContainer self_type;
    typedef PrimitiveContainer<std::string> value_type;
    typedef std::size_t             size_type;
    typedef value_type&             reference;
    typedef const value_type&       const_reference;
    typedef value_type*             pointer;
    typedef const value_type*       const_pointer;
    typedef std::ptrdiff_t          difference_type;
    typedef DataContainer           data_container_type;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
    PrimitiveGroupContainer();
    PrimitiveGroupContainer(const data_container_type& container, const uint32_t& offset, const uint32_t& n_entries, const uint32_t strides_each);
    ~PrimitiveGroupContainer(void);

    inline PrimitiveGroupContainerInterface* Clone(){ return(new self_type(*this)); }
	PrimitiveGroupContainerInterface& Move(PrimitiveGroupContainerInterface& src){
		self_type* o = reinterpret_cast<self_type*>(&src);
		this->n_capacity_ = o->n_capacity_;
		this->n_objects_  = o->n_objects_;

		for(std::size_t i = 0; i < this->size(); ++i)
			((this->containers_ + i)->~PrimitiveContainer)();

		::operator delete[](static_cast<void*>(this->containers_));

		this->containers_ = nullptr;

		std::swap(this->containers_, o->containers_);

		return(*this);
	}

	uint32_t BalanceVector(){
		if(this->size() == 0) return 0;

		uint32_t max_stride = 0;
		for(int i = 0; i < this->size(); ++i){
			if(this->at(i).data_.size() > max_stride)
				max_stride = this->at(i).data_.size();
		}

		return max_stride;
	}

	DataContainer ToDataContainer(void){
		if(this->size() == 0)
			return DataContainer();

		// Find the longest string in this Format record.
		const uint32_t stride    = this->BalanceVector();
		const uint32_t l_entries = this->size() * l_entries;

		DataContainer d;
		d.data_uncompressed.resize(l_entries + 128);
		d.strides_uncompressed.resize(l_entries + 128);
		d.header.data_header.stride = stride;

		for(uint32_t i = 0; i < this->size(); ++i){
			uint32_t delta = stride - this->at(i).data_.size();
			this->at(i).UpdateDataContainer(d, false);
			// If the target string is short then the largest string
			// available we pad the data added with null-terminations.
			// This will force the importing algorithm to skip premature
			// terminations while still maintaining the required matrix
			// properties of a Format record.
			for(int j = 0; j < delta; ++j)
				d.data_uncompressed += (char)'\0';

		}

		d.AddStride(stride);

		return(d);
	}

	DataContainer& UpdateDataContainer(DataContainer& container){
		if(this->size() == 0)
			return container;

		uint32_t l_entries = 0;
		uint32_t ref_size = this->at(0).data_.size();
		for(uint32_t i = 0; i < this->size(); ++i){
			assert(this->at(i).data_.size() == ref_size);
			l_entries += this->at(i).data_.size();
		}

		for(uint32_t i = 0; i < this->size(); ++i)
			this->at(i).UpdateDataContainer(container, false);

		container.AddStride(ref_size);

		return(container);
	}

	void resize(void);
	void resize(const size_t new_size);

	// Element access
	inline reference at(const size_type& position){ return(this->containers_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->containers_[position]); }
	inline reference operator[](const size_type& position){ return(this->containers_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->containers_[position]); }
	inline pointer data(void){ return(this->containers_); }
	inline const_pointer data(void) const{ return(this->containers_); }
	inline reference front(void){ return(this->containers_[0]); }
	inline const_reference front(void) const{ return(this->containers_[0]); }
	inline reference back(void){ return(this->containers_[this->n_objects_ - 1]); }
	inline const_reference back(void) const{ return(this->containers_[this->n_objects_ - 1]); }

	// Capacity
	inline bool empty(void) const{ return(this->n_objects_ == 0); }
	inline const size_type& size(void) const{ return(this->n_objects_); }

	// Iterator
	inline iterator begin(){ return iterator(&this->containers_[0]); }
	inline iterator end(){ return iterator(&this->containers_[this->n_objects_]); }
	inline const_iterator begin() const{ return const_iterator(&this->containers_[0]); }
	inline const_iterator end() const{ return const_iterator(&this->containers_[this->n_objects_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->containers_[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->containers_[this->n_objects_]); }

	io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint64_t position) const{
		this->at(position).ToVcfString(buffer);
		return(buffer);
	}

	bcf1_t* UpdateHtslibVcfRecordFormatInt32(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const{
		return(rec);
	}

	bcf1_t* UpdateHtslibVcfRecordFormatFloat(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const{
		return(rec);
	}

	bcf1_t* UpdateHtslibVcfRecordFormatString(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const{
		const char** dst = new const char*[this->size()];

		for(uint32_t i = 0; i < this->size(); ++i){
			dst[i] = this->at(i).data_.data();
			//std::cerr << dst[i] << std::endl;
		}

		bcf_update_format_string(hdr, rec, tag.data(), dst, this->size());

		delete [] dst;
		return(rec);
	}

private:
    pointer   containers_;
};

}
}



#endif /* CONTAINERS_PRIMITIVE_GROUP_CONTAINER_STRING_H_ */
