#ifndef CONTAINERS_INFO_CONTAINER_STRING_H_
#define CONTAINERS_INFO_CONTAINER_STRING_H_

#include "meta_container.h"
#include "stride_container.h"
#include "info_container.h"

namespace tachyon{
namespace containers{

/**<
 * InfoContainer for strings
 */
template <>
class InfoContainer<std::string> : public InfoContainerInterface{
private:
    typedef InfoContainer        self_type;
    typedef PrimitiveContainer<std::string> value_type;
    typedef value_type&          reference;
    typedef const value_type&    const_reference;
    typedef value_type*          pointer;
    typedef const value_type*    const_pointer;
    typedef std::ptrdiff_t       difference_type;
    typedef std::size_t          size_type;
    typedef io::BasicBuffer      buffer_type;
    typedef DataContainer        data_container_type;
    typedef MetaContainer        meta_container_type;
    typedef StrideContainer<uint32_t> stride_container_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
    InfoContainer();
    InfoContainer(const data_container_type& container);
    InfoContainer(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches);
    ~InfoContainer(void);

    // Element access
    inline reference at(const size_type& position){ return(this->containers_[position]); }
    inline const_reference at(const size_type& position) const{ return(this->containers_[position]); }
    inline reference operator[](const size_type& position){ return(this->containers_[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->containers_[position]); }
    inline pointer data(void){ return(this->containers_); }
    inline const_pointer data(void) const{ return(this->containers_); }
    inline reference front(void){ return(this->containers_[0]); }
    inline const_reference front(void) const{ return(this->containers_[0]); }
    inline reference back(void){ return(this->containers_[this->n_entries - 1]); }
    inline const_reference back(void) const{ return(this->containers_[this->n_entries - 1]); }

    // Capacity
    inline bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }

    // Iterator
    inline iterator begin(){ return iterator(&this->containers_[0]); }
    inline iterator end()  { return iterator(&this->containers_[this->n_entries]); }
    inline const_iterator begin()  const{ return const_iterator(&this->containers_[0]); }
    inline const_iterator end()    const{ return const_iterator(&this->containers_[this->n_entries]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->containers_[0]); }
    inline const_iterator cend()   const{ return const_iterator(&this->containers_[this->n_entries]); }

    inline std::ostream& ToVcfString(std::ostream& stream, const uint32_t position) const{ return(stream << this->at(position).data_); }
    inline io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint32_t position) const{ buffer += this->at(position).data_; return(buffer); }
    inline io::BasicBuffer& ToJsonString(io::BasicBuffer& buffer, const uint32_t position) const{
    	if(this->at(position).size() == 0){
    		buffer += "null";
    		return(buffer);
    	}
    	buffer += '"'; buffer += this->at(position).data_; buffer += '"';
    	return(buffer);
    }

    bool emptyPosition(const uint32_t& position) const{ return(this->at(position).empty()); }

private:
    // For mixed strides
    void Setup(const data_container_type& container);
    void SetupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches);
    // For fixed strides
	void Setup(const data_container_type& container, const uint32_t stride_size);
	void SetupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const uint32_t stride_size);

private:
    pointer containers_;
};

}
}



#endif /* CONTAINERS_INFO_CONTAINER_STRING_H_ */
