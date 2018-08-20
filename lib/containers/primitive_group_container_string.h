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
