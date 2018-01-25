#ifndef ITERATOR_ITERATORINTEGERREFERENCE_H_
#define ITERATOR_ITERATORINTEGERREFERENCE_H_

#include "../io/BasicBuffer.h"

namespace tachyon{
namespace iterator{

template <class return_primitive_type = U32>
class IteratorIntegerReference{
    typedef IteratorIntegerReference self_type;
    typedef io::BasicBuffer          buffer_type;
    typedef return_primitive_type    T;

public:
    IteratorIntegerReference(char* const data, const U32 length);
    IteratorIntegerReference(const self_type& other);
    IteratorIntegerReference(self_type&& other) noexcept;
    IteratorIntegerReference& operator=(const self_type& other);
    IteratorIntegerReference& operator=(self_type&& other) noexcept;
    virtual ~IteratorIntegerReference();

    // Iterators
    inline void advance(const U32& steps){ (*this) += steps; }
    inline void retreat(const U32& steps){ (*this) -= steps; }

    // Element access
    virtual T operator[](const U32& position) const =0;
    virtual T operator*(void) const =0;
    virtual T at(const U32& position) const =0;
    virtual T front(void) const =0;
    virtual T back(void) const =0;
    virtual std::vector<return_primitive_type> data(void) const =0;

    // Capacity
    inline const size_t& size(void) const{ return(this->n_entries); }
    inline const bool empty(void) const{ return(this->n_entries == 0); }

    // Basic mathematical functions
    virtual float mathMean(void) const =0;
    virtual float mathAverage(void) const =0;
    virtual U64 mathMax(void) const =0;
    virtual U64 mathMin(void) const =0;

    // Convertions to string
    //virtual void toString(std::ostream& stream = std::cout) =0;
    //virtual void toString(buffer_type& buffer) =0;

private:
    //friend std::ostream& operator<<(std::ostream& os, const self_type& iterator);

protected:
    size_t current_position;
    size_t n_entries;
    char*  __data;
};

template <class actual_primitive_type, class return_primitive_type = U32>
class IteratorIntegerReferenceImpl : public IteratorIntegerReference<return_primitive_type>{
private:
    typedef IteratorIntegerReferenceImpl self_type;
    typedef IteratorIntegerReference<return_primitive_type> parent_type;
    typedef io::BasicBuffer          buffer_type;
    typedef return_primitive_type    T;

public:
    IteratorIntegerReferenceImpl(char* const data, const U32 length);
    ~IteratorIntegerReferenceImpl();

    // Element access
    inline T operator[](const U32& position) const{ return(*reinterpret_cast<const actual_primitive_type* const>(&this->__data[position*sizeof(actual_primitive_type)])); }
    inline T operator*() const{ return(*reinterpret_cast<const actual_primitive_type* const>(&this->__data[this->current_position*sizeof(actual_primitive_type)])); }
    inline T at(const U32& position) const{ return(*reinterpret_cast<const actual_primitive_type* const>(&this->__data[position*sizeof(actual_primitive_type)])); }
    inline T front(void) const{ return(*reinterpret_cast<const actual_primitive_type* const>(&this->__data[0])); }
    inline T back(void) const{
        if(this->n_entries == 0)
            return(this->front());
        return(*reinterpret_cast<const actual_primitive_type* const>(&this->__data[(this->n_entries - 1)*sizeof(actual_primitive_type)]));
    }

    std::vector<return_primitive_type> data(void) const{
    	std::vector<return_primitive_type> ret;
    	ret.resize(this->n_entries);

    	for(U32 i = 0; i < this->n_entries; ++i)
    		ret[i] = this->at(i);

    	return(ret);
    }

    // Basic mathematical functions
    float mathMean(void) const{
        if(this->n_entries == 0) return 0;

        float total = 0;
        for(U32 i = 0; i < this->n_entries; ++i)
            total += this->at(i);

        return(total / this->n_entries);
    }

    inline float mathAverage(void) const{ return(this->mathMean()); }

    U64 mathMax(void) const{
        if(this->n_entries == 0) return 0;
        U64 current_max  = 0;

        for(U32 i = 0; i < this->n_entries; ++i)
            if(this->at(i) > current_max) current_max = this->at(i);

        return(current_max);
    }

    U64 mathMin(void) const{
        if(this->n_entries == 0) return 0;
        U64 current_min = std::numeric_limits<U64>::max();

        for(U32 i = 0; i < this->n_entries; ++i)
            if(this->at(i) < current_min) current_min = this->at(i);

        return(current_min);
    }

    // Convertions to string
    //void toString(std::ostream& stream = std::cout);
    //void toString(buffer_type& buffer);
};

template <class return_primitive_type>
IteratorIntegerReference<return_primitive_type>::IteratorIntegerReference(char* const data, const U32 length) :
	current_position(0),
	n_entries(length),
	__data(data)
{

}

template <class return_primitive_type>
IteratorIntegerReference<return_primitive_type>::~IteratorIntegerReference(){}

template <class actual_primitive_type, class return_primitive_type>
IteratorIntegerReferenceImpl<actual_primitive_type, return_primitive_type>::IteratorIntegerReferenceImpl(char* const data, const U32 length) :
	parent_type(data, length)
{
}

template <class actual_primitive_type, class return_primitive_type>
IteratorIntegerReferenceImpl<actual_primitive_type, return_primitive_type>::~IteratorIntegerReferenceImpl(){}


}
}



#endif /* ITERATOR_ITERATORINTEGERREFERENCE_H_ */
