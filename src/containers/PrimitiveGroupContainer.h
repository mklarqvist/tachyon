#ifndef CONTAINERS_PRIMITIVEGROUPCONTAINER_H_
#define CONTAINERS_PRIMITIVEGROUPCONTAINER_H_

namespace Tachyon{
namespace Core{

template <class T>
class PrimitiveGroupContainer{
private:
    typedef PrimitiveGroupContainer self_type;
    typedef PrimitiveContainer<T> value_type;
    typedef std::size_t           size_type;
    typedef value_type&           reference;
    typedef const value_type&     const_reference;
    typedef value_type*           pointer;
    typedef const value_type*     const_pointer;
    typedef std::ptrdiff_t        difference_type;
    typedef DataContainer         data_container_type;

public:
    // PrimitiveGroupContainer();
    PrimitiveGroupContainer(const data_container_type& container, const U32& offset, const U32 n_objects, const U32 strides_each);
    ~PrimitiveGroupContainer(void);

private:
    size_type __n_objects;
    pointer   __containers;
};

template <class T>
PrimitiveGroupContainer<T>::PrimitiveGroupContainer(const data_container_type& container, const U32& offset, const U32 n_objects, const U32 strides_each) :
	__n_objects(n_objects),
	__containers(static_cast<pointer>(::operator new[](this->__n_objects*sizeof(value_type))))
{
	U32 current_offset = offset;
	for(U32 i = 0; i < this->__n_objects; ++i){
		new( &this->__containers[i] ) value_type( container, current_offset, strides_each );
		current_offset += strides_each * sizeof(T);
	}
}

template <class T>
PrimitiveGroupContainer<T>::~PrimitiveGroupContainer(){
	for(std::size_t i = 0; i < this->__n_objects; ++i)
		((this->__containers + i)->~PrimitiveContainer)();

	::operator delete[](static_cast<void*>(this->__containers));
}

}
}



#endif /* CONTAINERS_PRIMITIVEGROUPCONTAINER_H_ */
