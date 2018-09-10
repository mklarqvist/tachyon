#ifndef CONTAINER_INFOCONTAINER_H_
#define CONTAINER_INFOCONTAINER_H_

#include "data_container.h"
#include "primitive_container.h"
#include "stride_container.h"
#include "meta_container.h"
#include "utility/support_vcf.h"
#include "info_container_interface.h"

namespace tachyon{
namespace containers{

template <class return_type>
class InfoContainer : public InfoContainerInterface{
public:
    typedef InfoContainer                   self_type;
    typedef PrimitiveContainer<return_type> value_type;
    typedef value_type&                     reference;
    typedef const value_type&               const_reference;
    typedef value_type*                     pointer;
    typedef const value_type*               const_pointer;
    typedef std::ptrdiff_t                  difference_type;
    typedef std::size_t                     size_type;
    typedef io::BasicBuffer                 buffer_type;
    typedef DataContainer                   data_container_type;
    typedef MetaContainer                   meta_container_type;
    typedef StrideContainer<uint32_t>            stride_container_type;

    typedef yonRawIterator<value_type>       iterator;
   	typedef yonRawIterator<const value_type> const_iterator;

public:
    InfoContainer();
    InfoContainer(const bool is_flag);
    InfoContainer(const data_container_type& container);
    InfoContainer(const data_container_type& data_container,
                  const meta_container_type& meta_container,
                  const std::vector<bool>& pattern_matches);
    ~InfoContainer(void);

    DataContainer ToDataContainer(void){
    	if(this->size() == 0)
    		return(DataContainer());

		uint32_t n_entries = 0;
		for(uint32_t i = 0; i < this->size(); ++i)
			n_entries += this->at(i).size();

		if(n_entries == 0)
			return(DataContainer());

		DataContainer d;
		d.data_uncompressed.resize(n_entries + 128);
		d.strides_uncompressed.resize(n_entries + 128);

		for(uint32_t i = 0; i < this->size(); ++i)
			this->at(i).UpdateDataContainer(d);

		return(d);
	}

	DataContainer& UpdateDataContainer(DataContainer& container){
		if(this->size() == 0)
			return(container);

		uint32_t n_entries = 0;
		for(uint32_t i = 0; i < this->size(); ++i)
			n_entries += this->at(i).size();

		if(n_entries == 0)
			return(container);

		if(container.data_uncompressed.size() + n_entries > container.data_uncompressed.capacity())
			container.data_uncompressed.resize((container.data_uncompressed.size()+n_entries)*2);

		for(uint32_t i = 0; i < this->size(); ++i)
			this->at(i).UpdateDataContainer(container);

		return(container);
	}

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

    // Iterator
    inline iterator begin(){ return iterator(&this->containers_[0]); }
    inline iterator end()  { return iterator(&this->containers_[this->n_entries]); }
    inline const_iterator begin()  const{ return const_iterator(&this->containers_[0]); }
    inline const_iterator end()    const{ return const_iterator(&this->containers_[this->n_entries]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->containers_[0]); }
    inline const_iterator cend()   const{ return const_iterator(&this->containers_[this->n_entries]); }

    // Type-specific
    inline std::ostream& ToVcfString(std::ostream& stream, const uint32_t position) const{
    	//utility::to_vcf_string(stream, this->at(position).data(), this->at(position).size());
    	return(stream);
    }

    inline io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint32_t position) const{
    	utility::ToVcfString(buffer, this->at(position).data(), this->at(position).size());
    	return(buffer);
    }

    inline io::BasicBuffer& ToJsonString(io::BasicBuffer& buffer, const uint32_t position) const{
    	//utility::to_json_string(buffer, this->at(position));
		return(buffer);
    }

    inline bool emptyPosition(const uint32_t& position) const{ return(this->at(position).empty()); }

private:
    // For mixed strides
    template <class actual_primitive>
    void Setup(const data_container_type& container);

    template <class actual_primitive>
	void SetupBalanced(const data_container_type& data_container,
	                   const meta_container_type& meta_container,
	                   const std::vector<bool>& pattern_matches);

	void SetupBalancedFlag(const data_container_type& data_container,
	                       const meta_container_type& meta_container,
	                       const std::vector<bool>& pattern_matches);


    // For fixed strides
	template <class actual_primitive>
	void Setup(const data_container_type& container, const uint32_t stride_size);

	template <class actual_primitive>
	void SetupBalanced(const data_container_type& data_container,
	                   const meta_container_type& meta_container,
	                   const std::vector<bool>& pattern_matches, const uint32_t stride_size);

private:
    pointer containers_;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_type>
InfoContainer<return_type>::InfoContainer(void) :
	containers_(nullptr)
{

}

template <class return_type>
InfoContainer<return_type>::InfoContainer(const bool is_flag) :
	containers_(static_cast<pointer>(::operator new[](1*sizeof(value_type))))
{
	this->n_entries  = 1;
	this->n_capacity = 1;
	// Set the primitive container value to 0. This
	// is required for the yon1_t structures to point
	// to something that is not simply a nullpointer.
	// It Has no other practical uses.
	new( &this->containers_[0] ) value_type( 0 );
}

template <class return_type>
InfoContainer<return_type>::InfoContainer(const data_container_type& data_container,
                                          const meta_container_type& meta_container,
                                            const std::vector<bool>& pattern_matches) :
	containers_(nullptr)
{
	if(data_container.data_uncompressed.size() == 0 && data_container.header.data_header.GetPrimitiveType() != YON_TYPE_BOOLEAN){
		return;
	}

	if(data_container.header.data_header.HasMixedStride()){
		if(data_container.header.data_header.IsSigned()){
			switch(data_container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->SetupBalanced<int8_t>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_16B):    (this->SetupBalanced<int16_t>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_32B):    (this->SetupBalanced<int32_t>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_64B):    (this->SetupBalanced<int64_t>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_FLOAT):  (this->SetupBalanced<float>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_DOUBLE): (this->SetupBalanced<double>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_BOOLEAN):(this->SetupBalancedFlag(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)data_container.header.data_header.controller.type << std::endl; return;
			}
		} else {
			switch(data_container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->SetupBalanced<uint8_t>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_16B):    (this->SetupBalanced<uint16_t>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_32B):    (this->SetupBalanced<uint32_t>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_64B):    (this->SetupBalanced<uint64_t>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_FLOAT):  (this->SetupBalanced<float>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_DOUBLE): (this->SetupBalanced<double>(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_BOOLEAN):(this->SetupBalancedFlag(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)data_container.header.data_header.controller.type << std::endl; return;
			}
		}
	} else {
		if(data_container.header.data_header.IsSigned()){
			switch(data_container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->SetupBalanced<int8_t>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_16B):    (this->SetupBalanced<int16_t>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_32B):    (this->SetupBalanced<int32_t>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_64B):    (this->SetupBalanced<int64_t>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_FLOAT):  (this->SetupBalanced<float>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_DOUBLE): (this->SetupBalanced<double>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_BOOLEAN):(this->SetupBalancedFlag(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)data_container.header.data_header.controller.type << std::endl; return;
			}
		} else {
			switch(data_container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->SetupBalanced<uint8_t>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_16B):    (this->SetupBalanced<uint16_t>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_32B):    (this->SetupBalanced<uint32_t>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_64B):    (this->SetupBalanced<uint64_t>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_FLOAT):  (this->SetupBalanced<float>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_DOUBLE): (this->SetupBalanced<double>(data_container, meta_container, pattern_matches, data_container.header.data_header.stride));  break;
			case(YON_TYPE_BOOLEAN):(this->SetupBalancedFlag(data_container, meta_container, pattern_matches));  break;
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)data_container.header.data_header.controller.type << std::endl; return;
			}
		}
	}
}

template <class return_type>
InfoContainer<return_type>::InfoContainer(const data_container_type& container) :
	containers_(nullptr)
{
	if(container.data_uncompressed.size() == 0 && container.header.data_header.GetPrimitiveType() != YON_TYPE_BOOLEAN)
		return;


	if(container.header.data_header.HasMixedStride()){
		if(container.header.data_header.IsSigned()){
			switch(container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->Setup<int8_t>(container));  break;
			case(YON_TYPE_16B):    (this->Setup<int16_t>(container));    break;
			case(YON_TYPE_32B):    (this->Setup<int32_t>(container));    break;
			case(YON_TYPE_64B):    (this->Setup<int64_t>(container));    break;
			case(YON_TYPE_FLOAT):  (this->Setup<float>(container));  break;
			case(YON_TYPE_DOUBLE): (this->Setup<double>(container)); break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return;
			}
		} else {
			switch(container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->Setup<uint8_t>(container));   break;
			case(YON_TYPE_16B):    (this->Setup<uint16_t>(container));    break;
			case(YON_TYPE_32B):    (this->Setup<uint32_t>(container));    break;
			case(YON_TYPE_64B):    (this->Setup<uint64_t>(container));    break;
			case(YON_TYPE_FLOAT):  (this->Setup<float>(container));  break;
			case(YON_TYPE_DOUBLE): (this->Setup<double>(container)); break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return;
			}
		}
	} else {
		if(container.header.data_header.IsSigned()){
			switch(container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->Setup<int8_t>(container, container.header.data_header.stride));  break;
			case(YON_TYPE_16B):    (this->Setup<int16_t>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_32B):    (this->Setup<int32_t>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_64B):    (this->Setup<int64_t>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_FLOAT):  (this->Setup<float>(container, container.header.data_header.stride));  break;
			case(YON_TYPE_DOUBLE): (this->Setup<double>(container, container.header.data_header.stride)); break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return;
			}
		} else {
			switch(container.header.data_header.GetPrimitiveType()){
			case(YON_TYPE_8B):     (this->Setup<uint8_t>(container, container.header.data_header.stride));   break;
			case(YON_TYPE_16B):    (this->Setup<uint16_t>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_32B):    (this->Setup<uint32_t>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_64B):    (this->Setup<uint64_t>(container, container.header.data_header.stride));    break;
			case(YON_TYPE_FLOAT):  (this->Setup<float>(container, container.header.data_header.stride));  break;
			case(YON_TYPE_DOUBLE): (this->Setup<double>(container, container.header.data_header.stride)); break;
			case(YON_TYPE_BOOLEAN):
			case(YON_TYPE_CHAR):
			case(YON_TYPE_STRUCT):
			case(YON_TYPE_UNKNOWN):
			default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return;

			}
		}
	}
}

template <class return_type>
InfoContainer<return_type>::~InfoContainer(){
	for(std::size_t i = 0; i < this->size(); ++i)
		((this->containers_ + i)->~PrimitiveContainer)();

	::operator delete[](static_cast<void*>(this->containers_));
}

template <class return_type>
template <class actual_primitive>
void InfoContainer<return_type>::Setup(const data_container_type& container){
	if(container.strides_uncompressed.size() == 0)
		return;

	// Worst case number of entries
	this->n_capacity = container.data_uncompressed.size() / sizeof(actual_primitive);
	if(this->capacity() == 0)
		return;

	this->containers_ = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
	stride_container_type strides(container);

	uint32_t current_offset = 0;
	uint32_t i = 0;
	while(true){
		new( &this->containers_[i] ) value_type( container, current_offset, strides[i] );
		current_offset += strides[i] * sizeof(actual_primitive);
		++this->n_entries;
		++i;

		// Break condition
		if(current_offset == container.data_uncompressed.size())
			break;

		// Assertion of critical error
		assert(current_offset < container.data_uncompressed.size());
	}
	assert(current_offset == container.data_uncompressed.size());
}

template <class return_type>
template <class actual_primitive>
void InfoContainer<return_type>::SetupBalanced(const data_container_type& data_container,
                                                 const meta_container_type& meta_container,
                                                   const std::vector<bool>& pattern_matches)
{
	this->n_entries = meta_container.size();
	if(this->size() == 0)
		return;

	this->containers_ = static_cast<pointer>(::operator new[](this->size()*sizeof(value_type)));
	stride_container_type strides(data_container);

	uint32_t current_offset = 0;
	uint32_t stride_offset  = 0;

	for(uint32_t i = 0; i < this->size(); ++i){
		// There are no INFO fields
		if(meta_container[i].GetInfoPatternId() == -1){
			new( &this->containers_[i] ) value_type( );
		}
		// If pattern matches
		else if(pattern_matches[meta_container[i].GetInfoPatternId()]){
			new( &this->containers_[i] ) value_type( data_container, current_offset, strides[stride_offset] );
			current_offset += strides[stride_offset] * sizeof(actual_primitive);
			++stride_offset;
		}
		// Otherwise place an empty
		else {
			new( &this->containers_[i] ) value_type( );
		}
	}

	assert(current_offset == data_container.data_uncompressed.size());
	assert(stride_offset == strides.size());
}

template <class return_type>
void InfoContainer<return_type>::SetupBalancedFlag(const data_container_type& data_container,
                                                 const meta_container_type& meta_container,
                                                   const std::vector<bool>& pattern_matches)
{
	this->n_entries = meta_container.size();
	std::cerr << "in flag ctor info: " << this->size() << std::endl;
	if(this->size() == 0)
		return;

	this->containers_ = static_cast<pointer>(::operator new[](this->size()*sizeof(value_type)));

	for(uint32_t i = 0; i < this->size(); ++i){
		// There are no INFO fields
		if(meta_container[i].GetInfoPatternId() == -1){
			new( &this->containers_[i] ) value_type( false );
		}
		// If pattern matches
		else if(pattern_matches[meta_container[i].GetInfoPatternId()]){
			std::cerr << "match add true: " << i << std::endl;
			new( &this->containers_[i] ) value_type( true );
		}
		// Otherwise place an empty
		else {
			new( &this->containers_[i] ) value_type( false );
		}
	}

}

template <class return_type>
template <class actual_primitive>
void InfoContainer<return_type>::Setup(const data_container_type& container, const uint32_t stride_size){
	this->n_entries = container.data_uncompressed.size() / sizeof(actual_primitive);

	if(this->size() == 0)
		return;

	this->containers_ = static_cast<pointer>(::operator new[](this->size()*sizeof(value_type)));

	uint32_t current_offset = 0;
	for(uint32_t i = 0; i < this->size(); ++i){
		//const actual_primitive* const data = reinterpret_cast<const actual_primitive* const>(&container.buffer_data_uncompressed[current_offset]);
		new( &this->containers_[i] ) value_type( container, current_offset, stride_size );
		current_offset += stride_size * sizeof(actual_primitive);
	}
	assert(current_offset == container.data_uncompressed.size());
}

template <class return_type>
template <class actual_primitive>
void InfoContainer<return_type>::SetupBalanced(const data_container_type& data_container,
                                                 const meta_container_type& meta_container,
                                                   const std::vector<bool>& pattern_matches,
                                                                 const uint32_t  stride_size)
{
	this->n_entries = meta_container.size();
	if(this->size() == 0)
		return;

	this->containers_ = static_cast<pointer>(::operator new[](this->size()*sizeof(value_type)));

	uint32_t current_offset = 0;
	// Case 1: if data is uniform
	if(data_container.header.data_header.IsUniform()){
		for(uint32_t i = 0; i < this->size(); ++i){
			// There are no INFO fields
			if(meta_container[i].GetInfoPatternId() == -1){
				new( &this->containers_[i] ) value_type( );
			}
			// If pattern matches
			else if(pattern_matches[meta_container[i].GetInfoPatternId()]){
				new( &this->containers_[i] ) value_type( data_container, 0, stride_size );
			} else {
				new( &this->containers_[i] ) value_type( );
			}
		}
		current_offset += stride_size * sizeof(actual_primitive);
	}
	// Case 2: if data is not uniform
	else {
		for(uint32_t i = 0; i < this->size(); ++i){
			// If pattern matches
			if(pattern_matches[meta_container[i].GetInfoPatternId()]){
				new( &this->containers_[i] ) value_type( data_container, current_offset, stride_size );
				current_offset += stride_size * sizeof(actual_primitive);
			}
			// Otherwise place an empty
			else {
				new( &this->containers_[i] ) value_type( );
			}
		}
	}
	assert(current_offset == data_container.data_uncompressed.size());
}

}
}


#endif /* InfoContainer_H_ */
