#ifndef CONTAINERS_PRIMITIVE_GROUP_CONTAINER_H_
#define CONTAINERS_PRIMITIVE_GROUP_CONTAINER_H_

#include "components/generic_iterator.h"
#include "primitive_container.h"

namespace tachyon{
namespace containers{

class PrimitiveGroupContainerInterface {
public:
	typedef PrimitiveGroupContainerInterface self_type;
	typedef std::size_t size_type;

public:
	PrimitiveGroupContainerInterface() : n_capacity_(0), n_objects_(0){ }
	PrimitiveGroupContainerInterface(const size_type n_objects) : n_capacity_(n_objects), n_objects_(n_objects){ }
	virtual ~PrimitiveGroupContainerInterface(){}

	/*
	virtual PrimitiveGroupContainerInterface* Clone() =0;
	virtual PrimitiveGroupContainerInterface& Move(PrimitiveGroupContainerInterface& other) =0;

	virtual void resize(void) =0;
	virtual void resize(const size_t new_size) =0;
	*/

	/**<
	 * Convert the data in a given PrimitiveContainer to valid Vcf
	 * formatting. Requires a valid array offset to the target
	 * container of interest.
	 * @param buffer Destination buffer.
	 * @return       Returns a reference to the destination buffer.
	 */
	virtual io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint64_t position) const =0;

	/**<
	 * Update a given htslib bcf1_t record destination Format field
	 * with the information in a given PrimitiveGroupContainer. Takes as
	 * arguments pointers to the dst bcf1_t record and global Vcf
	 * header and a reference to the dst tag name. The dst tag name
	 * has to be present in the bcf_hdr_t structure.
	 *
	 * The different UpdateHtslibVcfRecordFormat* functions is required
	 * to point to the different primitive types supported in the htslib
	 * vcf specification (int32, float, and special case of strings).
	 *
	 * @param rec Pointer to destination htslib vcf record.
	 * @param hdr Pointer to global Vcf header.
	 * @param tag Reference of destination tag string.
	 * @return    Returns a pointer to the input bcf1_t record.
	 */
	virtual bcf1_t* UpdateHtslibVcfRecordFormatInt32(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const =0;
	virtual bcf1_t* UpdateHtslibVcfRecordFormatFloat(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const =0;
	virtual bcf1_t* UpdateHtslibVcfRecordFormatString(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const =0;

	// Capacity
	inline bool empty(void) const{ return(this->n_objects_ == 0); }
	inline const size_type& size(void) const{ return(this->n_objects_); }

public:
	size_type n_capacity_;
    size_type n_objects_;
};

template <class return_type>
class PrimitiveGroupContainer : public PrimitiveGroupContainerInterface {
public:
    typedef PrimitiveGroupContainer self_type;
    typedef PrimitiveContainer<return_type> value_type;
    typedef std::size_t             size_type;
    typedef value_type&             reference;
    typedef const value_type&       const_reference;
    typedef value_type*             pointer;
    typedef const value_type*       const_pointer;
    typedef std::ptrdiff_t          difference_type;
    typedef DataContainer           data_container_type;

    typedef yonRawIterator<return_type>       iterator;
	typedef yonRawIterator<const return_type> const_iterator;

public:
    PrimitiveGroupContainer();
    PrimitiveGroupContainer(const data_container_type& container, const uint32_t& offset, const uint32_t n_objects, const uint32_t strides_each);
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

	// Iterator
	inline iterator begin(){ return iterator(&this->containers_[0]); }
	inline iterator end(){ return iterator(&this->containers_[this->n_objects_]); }
	inline const_iterator begin() const{ return const_iterator(&this->containers_[0]); }
	inline const_iterator end() const{ return const_iterator(&this->containers_[this->n_objects_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->containers_[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->containers_[this->n_objects_]); }

	inline io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint64_t position) const{
		return(utility::ToVcfString(buffer, this->at(position).data(), this->at(position).size()));
	}

	bcf1_t* UpdateHtslibVcfRecordFormatInt32(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const;
	bcf1_t* UpdateHtslibVcfRecordFormatFloat(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const;
	bcf1_t* UpdateHtslibVcfRecordFormatString(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const;

private:
	template <class actual_primitive_type>
	void Setup(const data_container_type& container, const uint32_t& offset, const uint32_t n_objects, const uint32_t strides_each);

private:
    pointer containers_;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_type>
PrimitiveGroupContainer<return_type>::PrimitiveGroupContainer() : containers_(nullptr){}

template <class return_type>
PrimitiveGroupContainer<return_type>::PrimitiveGroupContainer(const data_container_type& container,
                                                              const uint32_t& offset,
                                                              const uint32_t n_objects,
                                                              const uint32_t strides_each) :
	PrimitiveGroupContainerInterface(n_objects),
	containers_(static_cast<pointer>(::operator new[](this->n_objects_*sizeof(value_type))))
{

	if(container.header.data_header.IsSigned()){
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->Setup<int8_t>(container, offset, n_objects, strides_each));  break;
		case(YON_TYPE_16B):    (this->Setup<int16_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_32B):    (this->Setup<int32_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_64B):    (this->Setup<int64_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container, offset, n_objects, strides_each));   break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container, offset, n_objects, strides_each));  break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default: std::cerr << utility::timestamp("ERROR","PGC") << "Disallowed: " << container.header.data_header.GetPrimitiveType() << std::endl; return;
		}
	} else {
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->Setup<uint8_t>(container, offset, n_objects, strides_each));  break;
		case(YON_TYPE_16B):    (this->Setup<uint16_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_32B):    (this->Setup<uint32_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_64B):    (this->Setup<uint64_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container, offset, n_objects, strides_each));    break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container, offset, n_objects, strides_each));   break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default: std::cerr << utility::timestamp("ERROR","PGC") << "Disallowed: " << container.header.data_header.GetPrimitiveType() << std::endl; return;
		}
	}
}

template <class return_type>
PrimitiveGroupContainer<return_type>::~PrimitiveGroupContainer(){
	for(std::size_t i = 0; i < this->n_objects_; ++i)
		((this->containers_ + i)->~value_type)();

	::operator delete[](static_cast<void*>(this->containers_));
}

template <class return_type>
template <class actual_primitive_type>
void PrimitiveGroupContainer<return_type>::Setup(const data_container_type& container, const uint32_t& offset, const uint32_t n_objects, const uint32_t strides_each){
	uint32_t current_offset = offset;
	for(uint32_t i = 0; i < this->n_objects_; ++i){
		new( &this->containers_[i] ) value_type( container, current_offset, strides_each );
		current_offset += strides_each * sizeof(actual_primitive_type);
	}
}

template <class return_type>
bcf1_t* PrimitiveGroupContainer<return_type>::UpdateHtslibVcfRecordFormatInt32(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const{
	uint32_t n_records = 0;
	for(uint32_t i = 0; i < this->size(); ++i)
		n_records += this->at(i).size();

	int32_t* dst = new int32_t[n_records];
	uint32_t n_offset = 0;
	for(uint32_t i = 0; i < this->size(); ++i){
		utility::FormatDataHtslib(this->at(i).data(), &dst[n_offset], this->at(i).size());
		n_offset += this->at(i).size();
	}

	bcf_update_format_int32(hdr, rec, tag.data(), dst, n_records);

	delete [] dst;
	return(rec);
}

template <class return_type>
bcf1_t* PrimitiveGroupContainer<return_type>::UpdateHtslibVcfRecordFormatFloat(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const{
	uint32_t n_records = 0;
	for(uint32_t i = 0; i < this->size(); ++i)
		n_records += this->at(i).size();

	float* dst = new float[n_records];
	uint32_t n_offset = 0;
	for(uint32_t i = 0; i < this->size(); ++i){
		utility::FormatDataHtslib(this->at(i).data(), &dst[n_offset], this->at(i).size());
		n_offset += this->at(i).size();
	}

	bcf_update_format_float(hdr, rec, tag.data(), dst, n_records);

	delete [] dst;
	return(rec);
}

template <class return_type>
bcf1_t* PrimitiveGroupContainer<return_type>::UpdateHtslibVcfRecordFormatString(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const{
	std::cerr << utility::timestamp("ERROR","PGC") << "Illegal conversion from non-character primitive group container to string" << std::endl;
	exit(1);
	return(rec);
}

}
}



#endif /* CONTAINERS_PRIMITIVE_GROUP_CONTAINER_H_ */
