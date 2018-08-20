#ifndef CONTAINERS_FORMAT_CONTAINER_STRING_H_
#define CONTAINERS_FORMAT_CONTAINER_STRING_H_

#include "format_container.h"
#include "primitive_group_container_string.h"

namespace tachyon{
namespace containers{

/**<
 * Primary class for FORMAT data: special case for char
 * interpreted as std::string
 *
 * String containers have the special property that there
 * can only be one element. This is an intrinsic limitation
 * of the specification
 */
template <>
class FormatContainer<std::string>  : public FormatContainerInterface{
public:
    typedef FormatContainer       self_type;
    typedef PrimitiveGroupContainer<std::string> value_type;
    typedef value_type&           reference;
    typedef const value_type&     const_reference;
    typedef value_type*           pointer;
    typedef const value_type*     const_pointer;
    typedef std::ptrdiff_t        difference_type;
    typedef std::size_t           size_type;
    typedef io::BasicBuffer       buffer_type;
    typedef DataContainer         data_container_type;
    typedef MetaContainer         meta_container_type;
    typedef StrideContainer<uint32_t>  stride_container_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
    FormatContainer();
    FormatContainer(const data_container_type& container, const uint64_t n_samples);
    FormatContainer(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const uint64_t n_samples); // use when balancing
    ~FormatContainer(void);

    // Element access
    inline reference at(const size_type& position){ return(this->__containers[position]); }
    inline const_reference at(const size_type& position) const{ return(this->__containers[position]); }
    inline reference operator[](const size_type& position){ return(this->__containers[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->__containers[position]); }
    inline pointer data(void){ return(this->__containers); }
    inline const_pointer data(void) const{ return(this->__containers); }
    inline reference front(void){ return(this->__containers[0]); }
    inline const_reference front(void) const{ return(this->__containers[0]); }
    inline reference back(void){ return(this->__containers[this->n_entries - 1]); }
    inline const_reference back(void) const{ return(this->__containers[this->n_entries - 1]); }

    // Capacity
    inline bool empty(void) const{ return(this->n_entries == 0); }
    inline const size_type& size(void) const{ return(this->n_entries); }

    // Iterator
    inline iterator begin(){ return iterator(&this->__containers[0]); }
    inline iterator end()  { return iterator(&this->__containers[this->n_entries]); }
    inline const_iterator begin()  const{ return const_iterator(&this->__containers[0]); }
    inline const_iterator end()    const{ return const_iterator(&this->__containers[this->n_entries]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__containers[0]); }
    inline const_iterator cend()   const{ return const_iterator(&this->__containers[this->n_entries]); }

    // Type-specific
	inline std::ostream& ToVcfString(std::ostream& stream,
	                                 const uint32_t position,
	                                 const uint64_t sample) const
	{
		//utility::to_vcf_string(stream, this->at(position).at(sample).data_);
		//this->at(position).at(sample).data_;
		assert(2 == 1);
		return(stream);
	}

	bcf1_t* UpdateHtslibVcfRecord(const uint32_t position, bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const{
		return(this->at(position).UpdateHtslibVcfRecordFormatString(rec, hdr, tag));
	}

	inline io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint32_t position, const uint64_t sample) const{ buffer += this->at(position).at(sample).data_; return(buffer); }
	inline io::BasicBuffer& ToJsonString(io::BasicBuffer& buffer, const uint32_t position, const uint64_t sample) const{ buffer += this->at(position).at(sample).data_; return(buffer); }
	inline bool emptyPosition(const uint32_t& position) const{ return(this->at(position).empty()); }
	inline bool emptyPosition(const uint32_t& position, const uint64_t& sample) const{ return(this->at(position).at(sample).empty()); }

private:
    /**<
     *
     * @param container
     * @param n_samples
     */
    void __setup(const data_container_type& container, const uint64_t& n_samples);

    /**<
     *
     * @param data_container
     * @param meta_container
     * @param pattern_matches
     * @param n_samples
     */
	void __setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const uint64_t& n_samples);

    /**<
     *
     * @param data_container
     * @param meta_container
     * @param pattern_matches
     * @param n_samples
     * @param stride_size
     */
    void __setupBalanced(const data_container_type& data_container, const meta_container_type& meta_container, const std::vector<bool>& pattern_matches, const uint64_t& n_samples, const uint32_t stride_size);

    /**<
     *
     * @param container
     * @param n_samples
     * @param stride_size
     */
	void __setup(const data_container_type& container, const uint64_t& n_samples, const uint32_t stride_size);

private:
    pointer __containers;
};

}
}



#endif /* CONTAINERS_FORMAT_CONTAINER_STRING_H_ */
