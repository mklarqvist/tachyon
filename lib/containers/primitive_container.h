#ifndef CONTAINER_PRIMITIVECONTAINER_H_
#define CONTAINER_PRIMITIVECONTAINER_H_

#include <typeinfo>

#include "components/generic_iterator.h"
#include "variant_block.h"
#include "math/summary_statistics.h"
#include "utility/support_vcf.h"

namespace tachyon{
namespace containers{

class PrimitiveContainerInterface {
public:
	typedef PrimitiveContainerInterface self_type;
	typedef std::size_t size_type;

public:
	PrimitiveContainerInterface(void) : is_uniform_(false), n_entries_(0){}
	PrimitiveContainerInterface(const bool uniform, const size_t size) : is_uniform_(uniform), n_entries_(size){}
	virtual ~PrimitiveContainerInterface(){}

	 // Capacity
	inline bool empty(void) const{ return(this->n_entries_ == 0); }
	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline bool isUniform(void) const{ return(this->is_uniform_); }

	virtual io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer) const =0;

protected:
	bool    is_uniform_;
    size_t  n_entries_;
};

template <class return_type>
class PrimitiveContainer : public PrimitiveContainerInterface {
private:
    typedef std::size_t       size_type;
    typedef return_type       value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef std::ptrdiff_t    difference_type;
    typedef DataContainer     container_type;

    typedef yonRawIterator<value_type>       iterator;
    typedef yonRawIterator<const value_type> const_iterator;

public:
    PrimitiveContainer();
    PrimitiveContainer(const return_type value);
    PrimitiveContainer(const container_type& container);
    PrimitiveContainer(const container_type& container, const U32& offset, const U32 n_entries);
    ~PrimitiveContainer(void);

    // Element access
    inline reference at(const size_type& position){ return(this->__entries[position]); }
    inline const_reference at(const size_type& position) const{ return(this->__entries[position]); }
    inline reference operator[](const size_type& position){ return(this->__entries[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->__entries[position]); }
    inline pointer data(void){ return(this->__entries); }
    inline const_pointer data(void) const{ return(this->__entries); }
    inline reference front(void){ return(this->__entries[0]); }
    inline const_reference front(void) const{ return(this->__entries[0]); }
    inline reference back(void){ return(this->__entries[this->n_entries_ - 1]); }
    inline const_reference back(void) const{ return(this->__entries[this->n_entries_ - 1]); }


    // Iterator
    inline iterator begin(){ return iterator(&this->__entries[0]); }
    inline iterator end(){ return iterator(&this->__entries[this->n_entries_]); }
    inline const_iterator begin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator end() const{ return const_iterator(&this->__entries[this->n_entries_]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->__entries[this->n_entries_]); }

    io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer) const{
    	utility::to_vcf_string(buffer, this->data(), this->size());
    	return(buffer);
    }

private:
    template <class native_primitive>
    void __setup(const container_type& container, const U32& offset);

    template <class native_primitive>
    void __setupSigned(const container_type& container, const U32& offset);

private:
    pointer __entries;
};

template <>
class PrimitiveContainer<std::string> : public PrimitiveContainerInterface{
public:
    typedef PrimitiveContainer   self_type;
    typedef std::string          value_type;
    typedef value_type&          reference;
    typedef const value_type&    const_reference;
    typedef value_type*          pointer;
    typedef const value_type*    const_pointer;
    typedef std::ptrdiff_t       difference_type;
    typedef std::size_t          size_type;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
    PrimitiveContainer(){}
	PrimitiveContainer(const char* data, const size_t l_data) :
		PrimitiveContainerInterface(false, l_data),
		data_(data, l_data)
    {}
	~PrimitiveContainer(void){}

	// Element access
	inline pointer data(void){ return(&this->data_); }
	inline const_pointer data(void) const{ return(&this->data_); }
	const bool empty(void) const{ return(this->data_.size() == 0); }

	 // Iterator
	inline iterator begin(){ return iterator(&this->data_); }
	inline iterator end(){ return iterator(&this->data_ + this->n_entries_); }
	inline const_iterator begin() const{ return const_iterator(&this->data_); }
	inline const_iterator end() const{ return const_iterator(&this->data_ + this->n_entries_); }
	inline const_iterator cbegin() const{ return const_iterator(&this->data_); }
	inline const_iterator cend() const{ return const_iterator(&this->data_ + this->n_entries_); }

	io::BasicBuffer& to_vcf_string(io::BasicBuffer& buffer) const{
		if(this->data_.size() == 0){
			buffer += '.';
			return(buffer);
		}
		buffer += this->data_;
		return(buffer);
	}

public:
	std::string data_;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer() :
	__entries(nullptr)
{

}


template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer(const return_type value) :
	__entries(new return_type[1])
{
	this->n_entries_ = 1;
	this->__entries[0] = value;
}

template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer(const container_type& container) :
	__entries(nullptr)
{
	if(container.header.data_header.GetPrimitiveWidth() == -1)
		return;

	assert(container.buffer_data_uncompressed.size() % container.header.data_header.GetPrimitiveWidth() == 0);

	this->n_entries_ = container.buffer_data_uncompressed.size() / container.header.data_header.GetPrimitiveWidth();
	this->__entries = new value_type[this->n_entries_];

	if(this->n_entries_ == 0)
		return;

	this->is_uniform_ = container.header.data_header.IsUniform();

	if(container.header.data_header.IsSigned()){
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->__setupSigned<SBYTE>(container, 0));  break;
		case(YON_TYPE_16B):    (this->__setupSigned<S16>(container, 0));    break;
		case(YON_TYPE_32B):    (this->__setupSigned<S32>(container, 0));    break;
		case(YON_TYPE_64B):    (this->__setupSigned<S64>(container, 0));    break;
		case(YON_TYPE_FLOAT):  (this->__setup<float>(container, 0));        break;
		case(YON_TYPE_DOUBLE): (this->__setup<double>(container, 0));       break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default:
			std::cerr << "Disallowed: " << container.header.data_header.GetPrimitiveType() << std::endl;
			return;
		}
	} else {
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->__setup<BYTE>(container, 0));   break;
		case(YON_TYPE_16B):    (this->__setup<U16>(container, 0));    break;
		case(YON_TYPE_32B):    (this->__setup<U32>(container, 0));    break;
		case(YON_TYPE_64B):    (this->__setup<U64>(container, 0));    break;
		case(YON_TYPE_FLOAT):  (this->__setup<float>(container, 0));  break;
		case(YON_TYPE_DOUBLE): (this->__setup<double>(container, 0)); break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default:
			std::cerr << "Disallowed: " << container.header.data_header.GetPrimitiveType() << std::endl;
			return;
		}
	}
}

template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer(const container_type& container,
                                                              const U32& offset,
                                                              const U32  n_entries) :
    PrimitiveContainerInterface(false, n_entries),
	__entries(new value_type[n_entries])
{
	if(container.header.data_header.IsSigned()){
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->__setupSigned<SBYTE>(container, offset));  break;
		case(YON_TYPE_16B):    (this->__setupSigned<S16>(container, offset));    break;
		case(YON_TYPE_32B):    (this->__setupSigned<S32>(container, offset));    break;
		case(YON_TYPE_64B):    (this->__setupSigned<S64>(container, offset));    break;
		case(YON_TYPE_FLOAT):  (this->__setup<float>(container, offset));        break;
		case(YON_TYPE_DOUBLE): (this->__setup<double>(container, offset));       break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default: std::cerr << "Disallowed" << std::endl; return;
		}
	} else {
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->__setup<BYTE>(container, offset));   break;
		case(YON_TYPE_16B):    (this->__setup<U16>(container, offset));    break;
		case(YON_TYPE_32B):    (this->__setup<U32>(container, offset));    break;
		case(YON_TYPE_64B):    (this->__setup<U64>(container, offset));    break;
		case(YON_TYPE_FLOAT):  (this->__setup<float>(container, offset));  break;
		case(YON_TYPE_DOUBLE): (this->__setup<double>(container, offset)); break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default: std::cerr << "Disallowed" << std::endl; return;
		}
	}
}

template <class return_type>
PrimitiveContainer<return_type>::~PrimitiveContainer(void){
	delete [] this->__entries;
}

template <class return_type>
template <class native_primitive>
void PrimitiveContainer<return_type>::__setup(const container_type& container, const U32& offset){
	const native_primitive* const data = reinterpret_cast<const native_primitive* const>(&container.buffer_data_uncompressed.buffer[offset]);

	for(U32 i = 0; i < this->size(); ++i)
		this->__entries[i] = data[i];
}

template <class return_type>
template <class native_primitive>
void PrimitiveContainer<return_type>::__setupSigned(const container_type& container, const U32& offset){
	const native_primitive* const data = reinterpret_cast<const native_primitive* const>(&container.buffer_data_uncompressed.buffer[offset]);

	if(sizeof(native_primitive) == sizeof(return_type)){
		return(this->__setup<native_primitive>(container, offset));
	}
	else {
		for(U32 i = 0; i < this->size(); ++i){
			// If the data is missing in the native format
			if(data[i] == std::numeric_limits<native_primitive>::min()){
				this->__entries[i] = std::numeric_limits<return_type>::min();
			}
			// If the data is EOV in the native format
			else if(data[i] == std::numeric_limits<native_primitive>::min() + 1){
				//std::cerr << "is eov: " << std::bitset<sizeof(native_primitive)*8>(data[i]) << "\t" << std::bitset<sizeof(return_type)*8>(std::numeric_limits<return_type>::min() + 1) << std::endl;
				this->__entries[i] = std::numeric_limits<return_type>::min() + 1;
			}
			else
				this->__entries[i] = data[i];
		}
	}
}

}
}



#endif /* PrimitiveContainer_H_ */
