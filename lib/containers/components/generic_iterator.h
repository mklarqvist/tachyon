#ifndef CONTAINERS_COMPONENTS_GENERIC_ITERATOR_H_
#define CONTAINERS_COMPONENTS_GENERIC_ITERATOR_H_

#include <cstddef>
#include <iterator>

namespace tachyon{

//-------------------------------------------------------------------
// Raw iterator with random access
//-------------------------------------------------------------------
template<typename DataType>
class yonRawIterator : public std::iterator<std::random_access_iterator_tag,
                                            DataType,
                                            ptrdiff_t,
                                            DataType*,
                                            DataType&>
{
public:
    typedef yonRawIterator       self_type;
    typedef DataType             value_type;
    typedef value_type&          reference;
    typedef const value_type&    const_reference;
    typedef value_type*          pointer;
    typedef const value_type*    const_pointer;
    typedef std::ptrdiff_t       difference_type;
    typedef std::size_t          size_type;

public:
	yonRawIterator(DataType* ptr = nullptr){m_ptr = ptr;}
	yonRawIterator(const self_type& rawIterator) = default;
    ~yonRawIterator(){}

    self_type& operator=(const self_type& rawIterator) = default;
    self_type& operator=(pointer ptr){m_ptr = ptr;return (*this);}

    operator bool() const {
        if(m_ptr) return true;
        else return false;
    }

    bool operator==(const self_type& rawIterator)const{return (m_ptr == rawIterator.getConstPtr());}
    bool operator!=(const self_type& rawIterator)const{return (m_ptr != rawIterator.getConstPtr());}

    self_type& operator+=(const ptrdiff_t& movement){m_ptr += movement;return (*this);}
    self_type& operator-=(const ptrdiff_t& movement){m_ptr -= movement;return (*this);}
    self_type& operator++(){++m_ptr;return (*this);}
    self_type& operator--(){--m_ptr;return (*this);}
    self_type  operator++(int){auto temp(*this);++m_ptr;return temp;}
    self_type  operator--(int){auto temp(*this);--m_ptr;return temp;}
    self_type  operator+(const ptrdiff_t& movement){auto oldPtr = m_ptr;m_ptr+=movement;auto temp(*this);m_ptr = oldPtr;return temp;}
    self_type  operator-(const ptrdiff_t& movement){auto oldPtr = m_ptr;m_ptr-=movement;auto temp(*this);m_ptr = oldPtr;return temp;}

    ptrdiff_t operator-(const self_type& rawIterator){return std::distance(rawIterator.getPtr(),this->getPtr());}

    reference operator*(){return *m_ptr;}
    const_reference operator*()const{return *m_ptr;}
    pointer operator->(){return m_ptr;}
    pointer getPtr()const{return m_ptr;}
    const_pointer getConstPtr()const{return m_ptr;}

protected:
    pointer m_ptr;
};

}



#endif /* CONTAINERS_COMPONENTS_GENERIC_ITERATOR_H_ */
