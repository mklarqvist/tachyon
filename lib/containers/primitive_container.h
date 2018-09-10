#ifndef CONTAINER_PRIMITIVECONTAINER_H_
#define CONTAINER_PRIMITIVECONTAINER_H_

#include <typeinfo>

#include "components/generic_iterator.h"
#include "math/summary_statistics.h"
#include "utility/support_vcf.h"

namespace tachyon{
namespace containers{

class PrimitiveContainerInterface {
public:
	typedef PrimitiveContainerInterface self_type;
	typedef std::size_t size_type;

public:
	PrimitiveContainerInterface(void) :
		is_uniform_(false),
		n_entries_(0),
		n_capacity_(0)
	{

	}

	PrimitiveContainerInterface(const bool uniform, const size_t size) :
		is_uniform_(uniform),
		n_entries_(size),
		n_capacity_(size)
	{

	}

	virtual ~PrimitiveContainerInterface(){}

	/**<
	 * Virtual clone (copy) constructor for the situation when
	 * iteratively copying from an array of PrimitiveContainerInterface
	 * objects.
	 * @return Returns a pointer to the newly cloned object.
	 */
	virtual PrimitiveContainerInterface* Clone() =0;

	/**<
	 * Virtual move constructor for the situation when iteratively
	 * moving from an array of PrimitiveContainerInterface objects.
	 * The src interface reference will be reinterpreted to the
	 * correct derived class in the virtual function invocation in
	 * the derived class.
	 *
	 * Example:
	 *     PrimitiveContainerInterface** x;
	 *     PrimitiveContainerInterface** y;
	 *     //Move x[0] to y[0];
	 *     y[0]->Move(x[0]);
	 *
	 * @param src Reference to a PrimitiveContainerInterface.
	 * @return    Returns a reference to the moved object.
	 */
	virtual PrimitiveContainerInterface& Move(PrimitiveContainerInterface& src) =0;

	virtual void resize(void) =0;
	virtual void resize(const size_t new_size) =0;

	/**<
	 * Overloaded += functions for adding new data to a
	 * destination PrimitiveContainer. These functions have
	 * to convert a src primitive type to the appropriate
	 * dst primitive type: this may involve truncating bits
	 * if the dst word width is smaller than the src word
	 * width (for example uint32_t -> uint8_t). Similarly,
	 * adding signed primitives to the container has to be
	 * checked for missing and sentinel node symbols and
	 * converted to the equivalent symbol in the return
	 * primitive type space.
	 * @param entry Input (src) primitive type or string.
	 */
	virtual void operator+=(const int8_t entry)   =0;
	virtual void operator+=(const int16_t entry)  =0;
	virtual void operator+=(const int32_t entry)  =0;
	virtual void operator+=(const int64_t entry)  =0;
	virtual void operator+=(const uint8_t entry)  =0;
	virtual void operator+=(const uint16_t entry) =0;
	virtual void operator+=(const uint32_t entry) =0;
	virtual void operator+=(const uint64_t entry) =0;
	virtual void operator+=(const char entry)     =0;
	virtual void operator+=(const float entry)    =0;
	virtual void operator+=(const double entry)   =0;
	virtual void operator+=(const std::string& entry) =0;

	/**<
	 * Convert a PrimitiveContainer into a DataContainer. This is primarily
	 * done for writing out a Tachyon archive.
	 * @return Returns a DataContainer with the contextual representation of the data into this container.
	 */
	virtual DataContainer ToDataContainer(void) =0;

	/**<
	 * Add data from a PrimitiveContainer into an already existing DataContainer.
	 * This is primarily done when iterative updating a DataContainer as for example
	 * when iterating over PrimitiveGroupContainer classes.
	 * @param container Destination DataContainer.
	 * @return          Returns a reference to the input DataContainer with the contextual representation of the data in this container added to it.
	 */
	virtual DataContainer& UpdateDataContainer(DataContainer& container) =0;

	// Capacity
	inline bool empty(void) const{ return(this->n_entries_ == 0); }
	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }
	inline bool IsUniform(void) const{ return(this->is_uniform_); }

	/**<
	 * Convert the data in a given PrimitiveContainer to valid Vcf
	 * formatting.
	 * @param buffer Destination buffer.
	 * @return       Returns a reference to the destination buffer.
	 */
	virtual io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer) const =0;

	/**<
	 * Update a given htslib bcf1_t record destination Info field
	 * with the information in a given PrimitiveContainer. Takes as
	 * arguments pointers to the dst bcf1_t record and global Vcf
	 * header and a reference to the dst tag name. The dst tag name
	 * has to be present in the bcf_hdr_t structure.
	 * @param rec Pointer to destination htslib vcf record.
	 * @param hdr Pointer to global Vcf header.
	 * @param tag Reference of destination tag string.
	 * @return    Returns a pointer to the input bcf1_t record.
	 */
	virtual bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
	                                          bcf_hdr_t* hdr,
	                                          const std::string& tag) const =0;

protected:
	bool    is_uniform_;
    size_t  n_entries_;
    size_t  n_capacity_;
};

template <class return_type>
class PrimitiveContainer : public PrimitiveContainerInterface {
public:
	typedef PrimitiveContainer self_type;
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
    PrimitiveContainer(const container_type& container,
                       const uint32_t& offset,
                       const uint32_t n_entries);
    ~PrimitiveContainer(void);
    PrimitiveContainer(const self_type& other);
    PrimitiveContainer(self_type&& other) noexcept;
    PrimitiveContainer& operator=(const self_type& other);
    PrimitiveContainer& operator=(self_type&& other) noexcept;

    inline PrimitiveContainerInterface* Clone(){ return(new self_type(*this)); }
    PrimitiveContainerInterface& Move(PrimitiveContainerInterface& src){
    	self_type* o   = reinterpret_cast<self_type*>(&src);
    	this->is_uniform_ = o->is_uniform_;
    	this->n_capacity_ = o->n_capacity_;
    	this->n_entries_  = o->n_entries_;
    	this->entries_    = nullptr;
    	std::swap(this->entries_, o->entries_);
    	return(*this);
    }

    void resize(void);
	void resize(const size_t new_size);

	void operator+=(const int8_t entry);
	void operator+=(const int16_t entry);
	void operator+=(const int32_t entry);
	void operator+=(const int64_t entry);
	void operator+=(const uint8_t entry);
	void operator+=(const uint16_t entry);
	void operator+=(const uint32_t entry);
	void operator+=(const uint64_t entry);
	void operator+=(const char entry);
	void operator+=(const float entry);
	void operator+=(const double entry);
	void operator+=(const std::string& entry);

	DataContainer ToDataContainer(void){
		if(this->size() == 0)
			return(DataContainer());

		DataContainer d;
		d.data_uncompressed.resize(this->size() + 128);
		d.strides_uncompressed.resize(this->size() + 128);

		for(uint32_t i = 0; i < this->size(); ++i)
			d.Add(this->at(i));

		d.AddStride(this->size());
		++d;

		return(d);
	}

	DataContainer& UpdateDataContainer(DataContainer& container){
		if(this->size() == 0)
			return(container);

		if(container.data_uncompressed.size() + this->size() > container.data_uncompressed.capacity())
			container.data_uncompressed.resize((container.data_uncompressed.size()+this->size())*2);

		for(uint32_t i = 0; i < this->size(); ++i)
			container.Add(this->at(i));

		container.AddStride(this->size());
		++container;

		return(container);
	}

    // Element access
    inline reference at(const size_type& position){ return(this->entries_[position]); }
    inline const_reference at(const size_type& position) const{ return(this->entries_[position]); }
    inline reference operator[](const size_type& position){ return(this->entries_[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->entries_[position]); }
    inline pointer data(void){ return(this->entries_); }
    inline const_pointer data(void) const{ return(this->entries_); }
    inline reference front(void){ return(this->entries_[0]); }
    inline const_reference front(void) const{ return(this->entries_[0]); }
    inline reference back(void){ return(this->entries_[this->n_entries_ - 1]); }
    inline const_reference back(void) const{ return(this->entries_[this->n_entries_ - 1]); }

    // Iterator
    inline iterator begin(){ return iterator(&this->entries_[0]); }
    inline iterator end(){ return iterator(&this->entries_[this->n_entries_]); }
    inline const_iterator begin() const{ return const_iterator(&this->entries_[0]); }
    inline const_iterator end() const{ return const_iterator(&this->entries_[this->n_entries_]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->entries_[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->entries_[this->n_entries_]); }

    inline io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer) const{
    	return(utility::ToVcfString(buffer, this->data(), this->size()));
    }

    inline bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
                                             bcf_hdr_t* hdr,
	                                         const std::string& tag) const
    {
    	return(utility::UpdateHtslibVcfRecordInfo(rec, hdr, tag, this->data(), this->size()));
    }

private:
    template <class native_primitive>
    void Setup(const container_type& container, const uint32_t& offset);

    template <class native_primitive>
    void SetupSigned(const container_type& container, const uint32_t& offset);

private:
    pointer entries_;
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
    PrimitiveContainer(){ this->n_capacity_ = 1; }
	PrimitiveContainer(const char* data, const size_t l_data) :
		PrimitiveContainerInterface(false, 1),
		data_(data, l_data)
    {

    }
	~PrimitiveContainer(void){}

	inline PrimitiveContainerInterface* Clone(){ return(new self_type(*this)); }
	PrimitiveContainerInterface& Move(PrimitiveContainerInterface& src){
		self_type* o   = reinterpret_cast<self_type*>(&src);
		this->is_uniform_ = o->is_uniform_;
		this->n_capacity_ = o->n_capacity_;
		this->n_entries_  = o->n_entries_;
		this->data_ = std::move(o->data_);
		return(*this);
	}

	inline void resize(void){}
	inline void resize(const size_t new_size){}

	inline void operator+=(const int8_t entry){ this->data_ += std::to_string(entry); }
	inline void operator+=(const int16_t entry){ this->data_ += std::to_string(entry); }
	inline void operator+=(const int32_t entry){ this->data_ += std::to_string(entry); }
	inline void operator+=(const int64_t entry){ this->data_ += std::to_string(entry); }
	inline void operator+=(const uint8_t entry){ this->data_ += std::to_string(entry); };
	inline void operator+=(const uint16_t entry){ this->data_ += std::to_string(entry); }
	inline void operator+=(const uint32_t entry){ this->data_ += std::to_string(entry); }
	inline void operator+=(const uint64_t entry){ this->data_ += std::to_string(entry); }
	inline void operator+=(const char entry){ this->data_ += std::to_string(entry); }
	inline void operator+=(const float entry){ this->data_ += std::to_string(entry); }
	inline void operator+=(const double entry){ this->data_ += std::to_string(entry); }
	inline void operator+=(const std::string& entry){ this->data_ += entry; }

	DataContainer ToDataContainer(void){
		DataContainer d;
		d.data_uncompressed.resize(this->size() + 128);
		d.strides_uncompressed.resize(this->size() + 128);

		for(uint32_t i = 0; i < this->size(); ++i){
			d.AddString(this->data_);
			d.AddStride(this->data_.size());
			++d;
		}

		return(d);
	}

	DataContainer& UpdateDataContainer(DataContainer& container){
		if(container.data_uncompressed.size() + this->size() > container.data_uncompressed.capacity())
			container.data_uncompressed.resize((container.data_uncompressed.size()+this->size())*2);

		for(uint32_t i = 0; i < this->size(); ++i){
			container.AddString(this->data_);
			container.AddStride(this->data_.size());
			++container;
		}

		return(container);
	}

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

	io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer) const{
		if(this->data_.size() == 0){
			buffer += '.';
			return(buffer);
		}
		buffer += this->data_;
		return(buffer);
	}

	inline bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
	                                         bcf_hdr_t* hdr,
	                                         const std::string& tag) const
	{
		return(utility::UpdateHtslibVcfRecordInfo(rec, hdr, tag, this->data_));
	}

public:
	std::string data_;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer() :
	entries_(nullptr)
{

}


template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer(const return_type value) :
	entries_(new return_type[1])
{
	this->n_entries_ = 1;
	this->entries_[0] = value;
}

template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer(const container_type& container) :
	entries_(nullptr)
{
	if(container.header.data_header.GetPrimitiveWidth() == -1)
		return;

	assert(container.data_uncompressed.size() % container.header.data_header.GetPrimitiveWidth() == 0);

	this->n_entries_ = container.data_uncompressed.size() / container.header.data_header.GetPrimitiveWidth();
	this->entries_ = new value_type[this->n_entries_];

	if(this->n_entries_ == 0)
		return;

	this->is_uniform_ = container.header.data_header.IsUniform();

	if(container.header.data_header.IsSigned()){
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->SetupSigned<int8_t>(container, 0));  break;
		case(YON_TYPE_16B):    (this->SetupSigned<int16_t>(container, 0)); break;
		case(YON_TYPE_32B):    (this->SetupSigned<int32_t>(container, 0)); break;
		case(YON_TYPE_64B):    (this->SetupSigned<int64_t>(container, 0)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container, 0));         break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container, 0));        break;
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
		case(YON_TYPE_8B):     (this->Setup<uint8_t>(container, 0));  break;
		case(YON_TYPE_16B):    (this->Setup<uint16_t>(container, 0)); break;
		case(YON_TYPE_32B):    (this->Setup<uint32_t>(container, 0)); break;
		case(YON_TYPE_64B):    (this->Setup<uint64_t>(container, 0)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container, 0));    break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container, 0));   break;
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
                                                    const uint32_t& offset,
                                                    const uint32_t  n_entries) :
    PrimitiveContainerInterface(false, n_entries),
	entries_(new value_type[n_entries])
{
	if(container.header.data_header.IsSigned()){
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->SetupSigned<int8_t>(container, offset));  break;
		case(YON_TYPE_16B):    (this->SetupSigned<int16_t>(container, offset)); break;
		case(YON_TYPE_32B):    (this->SetupSigned<int32_t>(container, offset)); break;
		case(YON_TYPE_64B):    (this->SetupSigned<int64_t>(container, offset)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container, offset));         break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container, offset));        break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default: std::cerr << "Disallowed" << std::endl; return;
		}
	} else {
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->Setup<uint8_t>(container, offset));  break;
		case(YON_TYPE_16B):    (this->Setup<uint16_t>(container, offset)); break;
		case(YON_TYPE_32B):    (this->Setup<uint32_t>(container, offset)); break;
		case(YON_TYPE_64B):    (this->Setup<uint64_t>(container, offset)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container, offset));    break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container, offset));   break;
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
	delete [] this->entries_;
}

template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer(const self_type& other) :
PrimitiveContainerInterface(other),
	entries_(new value_type[this->n_capacity_])
{
	memcpy(this->entries_, other.entries_, sizeof(value_type)*this->size());
}

template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer(self_type&& other) noexcept  :
	PrimitiveContainerInterface(other),
	entries_(nullptr)
{
	std::swap(this->entries_, other.entries_);
}

template <class return_type>
PrimitiveContainer<return_type>& PrimitiveContainer<return_type>::operator=(const self_type& other){
	delete [] this->entries_;
	this->is_uniform_ = other.is_uniform_;
	this->n_entries_  = other.n_entries_;
	this->n_capacity_ = other.n_capacity_;
	this->entries_    = new value_type[this->capacity()];
	memcpy(this->entries_, other.entries_, sizeof(value_type)*this->size());
	return(*this);
}

template <class return_type>
PrimitiveContainer<return_type>& PrimitiveContainer<return_type>::operator=(self_type&& other) noexcept{
	if(this == &other){
		// precautions against self-moves
		return *this;
	}

	delete [] this->entries_; this->entries_ = nullptr;
	this->is_uniform_ = other.is_uniform_;
	this->n_entries_  = other.n_entries_;
	this->n_capacity_ = other.n_capacity_;
	std::swap(this->entries_, other.entries_);
	return(*this);
}

template <class return_type>
template <class native_primitive>
void PrimitiveContainer<return_type>::Setup(const container_type& container, const uint32_t& offset){
	const native_primitive* const data = reinterpret_cast<const native_primitive* const>(&container.data_uncompressed.data()[offset]);

	for(uint32_t i = 0; i < this->size(); ++i)
		this->entries_[i] = data[i];
}

template <class return_type>
template <class native_primitive>
void PrimitiveContainer<return_type>::SetupSigned(const container_type& container, const uint32_t& offset){
	const native_primitive* const data = reinterpret_cast<const native_primitive* const>(&container.data_uncompressed.data()[offset]);

	if(sizeof(native_primitive) == sizeof(return_type)){
		return(this->Setup<native_primitive>(container, offset));
	}
	else {
		for(uint32_t i = 0; i < this->size(); ++i){
			// If the data is missing in the native format.
			if(data[i] == std::numeric_limits<native_primitive>::min()){
				this->entries_[i] = std::numeric_limits<return_type>::min();
			}
			// If the data is EOV in the native format.
			else if(data[i] == std::numeric_limits<native_primitive>::min() + 1){
				this->entries_[i] = std::numeric_limits<return_type>::min() + 1;
			}
			else
				this->entries_[i] = data[i];
		}
	}
}

template <class return_type>
void PrimitiveContainer<return_type>::resize(void){
	pointer temp       = this->entries_;
	this->entries_     = new value_type[this->n_capacity_*2];
	this->n_capacity_ *= 2;
	memcpy(this->entries_, temp, sizeof(value_type)*this->n_entries_);
	delete [] temp;
}

template <class return_type>
void PrimitiveContainer<return_type>::resize(const size_t new_size){
	// if new size < current capacity
	if(new_size < this->n_capacity_){
		// if new size < current number of entries
		if(new_size < this->n_entries_){
			this->n_entries_ = new_size;
			return;
		}
		return;
	}

	pointer temp      = this->entries_;
	this->entries_    = new value_type[new_size];
	this->n_capacity_ = new_size;
	memcpy(this->entries_, temp, sizeof(value_type)*this->n_entries_);
	delete [] temp;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const int8_t entry){
	if(this->n_entries_ + 1 == this->n_capacity_)
		this->resize();

	if(entry == INT8_MIN) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min();
	else if(entry == INT8_MIN+1) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min()+1;
	else this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const int16_t entry){
	if(this->n_entries_ + 1 == this->n_capacity_)
		this->resize();

	if(entry == INT16_MIN) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min();
	else if(entry == INT16_MIN+1) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min()+1;
	else this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const int32_t entry){
	if(this->n_entries_ + 1 == this->n_capacity_)
		this->resize();

	if(entry == INT32_MIN) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min();
	else if(entry == INT32_MIN+1) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min()+1;
	else this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const int64_t entry){
	if(this->n_entries_ + 1 == this->n_capacity_)
		this->resize();

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const uint8_t entry){
	if(this->n_entries_ + 1 == this->n_capacity_)
		this->resize();

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const uint16_t entry){
	if(this->n_entries_ + 1 == this->n_capacity_)
		this->resize();

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const uint32_t entry){
	if(this->n_entries_ + 1 == this->n_capacity_)
		this->resize();

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const uint64_t entry){
	if(this->n_entries_ + 1 == this->n_capacity_)
		this->resize();

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const char entry){
	if(this->n_entries_ + 1 == this->n_capacity_)
		this->resize();

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const float entry){
	if(this->n_entries_ + 1 == this->n_capacity_)
		this->resize();

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const double entry){
	if(this->n_entries_ + 1 == this->n_capacity_)
		this->resize();

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const std::string& entry){}

}
}



#endif /* PrimitiveContainer_H_ */
