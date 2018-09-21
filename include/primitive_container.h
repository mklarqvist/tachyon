#ifndef CONTAINER_PRIMITIVECONTAINER_H_
#define CONTAINER_PRIMITIVECONTAINER_H_

#include <typeinfo>

#include "data_container.h"
#include "containers/components/generic_iterator.h"
#include "support_vcf.h"

namespace tachyon {

class PrimitiveContainerInterface {
public:
	typedef PrimitiveContainerInterface self_type;
	typedef std::size_t size_type;

public:
	PrimitiveContainerInterface(void);
	PrimitiveContainerInterface(const bool uniform, const size_t size);
	virtual ~PrimitiveContainerInterface();

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
	virtual PrimitiveContainerInterface* Move(void) =0;

	/**<
	 * Expand records up to some number. This is used primarily when
	 * balancing vector of vectors (matrices) for writing as both yon
	 * and bcf requires matrices to be balanced. Data is always padded
	 * with the appropriate EOV value.
	 * @param to Target number of items to expand to.
	 */
	virtual void ExpandEmpty(const uint32_t to) =0;

	virtual void resize(void) =0;
	virtual void resize(const size_t new_size) =0;

	inline void operator++(void){ ++this->n_entries_; }
	void operator--(void){ if(this->n_entries_ != 0) --this->n_entries_; }

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
	 *
	 * Conversions have to take place when the dst container
	 * is signed and the primitive types are different. For
	 * example, this container is int32_t and adding a
	 * int8_t. In this case, the special values for EOV and
	 * missing have to be expanded to the new, larger, bit
	 * space.
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
	virtual yon1_dc_t ToDataContainer(void) =0;

	/**<
	 * Add data from a PrimitiveContainer into an already existing DataContainer.
	 * This is primarily done when iterative updating a DataContainer as for example
	 * when iterating over PrimitiveGroupContainer classes.
	 * @param container Destination DataContainer.
	 * @return          Returns a reference to the input DataContainer with the contextual representation of the data in this container added to it.
	 */
	virtual yon1_dc_t& UpdateDataContainer(yon1_dc_t& container, const bool update_stride) =0;

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

	/**<
	 * Return the word width of the templated return type in the implementation
	 * classes. This is useful for knowing the byte width of objects stored in the
	 * specialized children without having knowledge of their explicit typing.
	 * @return Returns the byte width of the return type.
	 */
	virtual uint8_t GetWordWidth(void) =0;

protected:
	bool    is_uniform_;
    size_t  n_entries_;
    size_t  n_capacity_;
};

template <class return_type>
class PrimitiveContainer : public PrimitiveContainerInterface {
public:
	typedef PrimitiveContainer self_type;
    typedef std::size_t        size_type;
    typedef return_type        value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;
    typedef std::ptrdiff_t     difference_type;
    typedef yon1_dc_t          container_type;

    typedef yonRawIterator<value_type>       iterator;
    typedef yonRawIterator<const value_type> const_iterator;

public:
    PrimitiveContainer();
    PrimitiveContainer(const return_type value);
    PrimitiveContainer(const container_type& container);
    PrimitiveContainer(const container_type& container,
                       const uint32_t& offset,
                       const uint32_t n_entries);
    PrimitiveContainer(const self_type& other);
    PrimitiveContainer(return_type* values, const uint32_t n_entries);
    PrimitiveContainer(self_type&& other) noexcept;
    PrimitiveContainer& operator=(const self_type& other);
    PrimitiveContainer& operator=(self_type&& other) noexcept;
    ~PrimitiveContainer(void);

    inline PrimitiveContainerInterface* Clone(){ return(new self_type(*this)); }

    PrimitiveContainerInterface* Move(void);

    void resize(void);
	void resize(const size_t new_size);

	// Add data to this object.
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

	yon1_dc_t ToDataContainer(void);
	yon1_dc_t& UpdateDataContainer(yon1_dc_t& container, const bool update_stride = true);
	void ExpandEmpty(const uint32_t to);

	inline uint8_t GetWordWidth(void){ return(sizeof(return_type)); }

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

// Partial specialization for std::string
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
    PrimitiveContainer(void);
	PrimitiveContainer(const char* data, const size_t l_data);
	PrimitiveContainer(const PrimitiveContainer& other);
	PrimitiveContainer& operator=(const PrimitiveContainer& other);
	~PrimitiveContainer(void);

	inline PrimitiveContainerInterface* Clone(){ return(new self_type(*this)); }
	PrimitiveContainerInterface* Move(void);
	void ExpandEmpty(const uint32_t to);

	inline uint8_t GetWordWidth(void){ return(sizeof(char)); }

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

	yon1_dc_t ToDataContainer(void);
	yon1_dc_t& UpdateDataContainer(yon1_dc_t& container, const bool update_stride = true);

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

	io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer) const;

	inline bcf1_t* UpdateHtslibVcfRecordInfo(bcf1_t* rec,
	                                         bcf_hdr_t* hdr,
	                                         const std::string& tag) const
	{
		return(utility::UpdateHtslibVcfRecordInfo(rec, hdr, tag, this->data_));
	}

public:
	std::string data_;
};

class PrimitiveGroupContainerInterface {
public:
	typedef PrimitiveGroupContainerInterface self_type;
	typedef std::size_t size_type;

public:
	PrimitiveGroupContainerInterface();
	PrimitiveGroupContainerInterface(const size_type n_objects);
	virtual ~PrimitiveGroupContainerInterface();

	virtual void resize(void) =0;
	virtual void resize(const size_t new_size) =0;

	/**<
	 * Virtual clone (copy) constructor for the situation when
	 * iteratively copying from an array of PrimitiveGroupContainerInterface
	 * objects.
	 * @return Returns a pointer to the newly cloned object.
	 */
	virtual PrimitiveGroupContainerInterface* Clone() =0;

	/**<
	 * Virtual move constructor for the situation when iteratively
	 * moving from an array of PrimitiveGroupContainerInterface objects.
	 * The src interface reference will be reinterpreted to the
	 * correct derived class in the virtual function invocation in
	 * the derived class.
	 *
	 * Example:
	 *     PrimitiveGroupContainerInterface** x;
	 *     PrimitiveGroupContainerInterface** y;
	 *     //Move x[0] to y[0];
	 *     y[0]->Move(x[0]);
	 *
	 * @param src Reference to a PrimitiveGroupContainerInterface.
	 * @return    Returns a reference to the moved object.
	 */
	virtual PrimitiveGroupContainerInterface* Move(void) =0;

	/**<
	 * Convert a PrimitiveGroupContainer into a DataContainer. This is primarily
	 * done for writing out a Tachyon archive.
	 * @return Returns a DataContainer with the contextual representation of the data into this container.
	 */
	virtual yon1_dc_t ToDataContainer(void) =0;

	/**<
	 * Add data from a PrimitiveGroupContainer into an already existing DataContainer.
	 * @param container Destination DataContainer.
	 * @return          Returns a reference to the input DataContainer with the contextual representation of the data in this container added to it.
	 */
	virtual yon1_dc_t& UpdateDataContainer(yon1_dc_t& container) =0;


	/**<
	 * Convert the data in a given PrimitiveGroupContainer to valid Vcf
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
    typedef yon1_dc_t           data_container_type;

    typedef yonRawIterator<return_type>       iterator;
	typedef yonRawIterator<const return_type> const_iterator;

public:
    PrimitiveGroupContainer();
    PrimitiveGroupContainer(const self_type& other);
    PrimitiveGroupContainer(const data_container_type& container,
                            const uint32_t& offset,
	                        const uint32_t n_objects,
	                        const uint32_t strides_each);
    ~PrimitiveGroupContainer(void);

    inline PrimitiveGroupContainerInterface* Clone(){ return(new self_type(*this)); }

    PrimitiveGroupContainerInterface* Move(void);

    yon1_dc_t ToDataContainer(void);
    yon1_dc_t& UpdateDataContainer(yon1_dc_t& container);

    void resize(void);
	void resize(const size_t new_size);

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
	void Setup(const data_container_type& container,
	           const uint32_t& offset,
	           const uint32_t n_objects,
	           const uint32_t strides_each);

private:
    pointer containers_;
};

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
    typedef yon1_dc_t               data_container_type;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
    PrimitiveGroupContainer();
    PrimitiveGroupContainer(const self_type& other);
    PrimitiveGroupContainer(const data_container_type& container,
                            const uint32_t& offset,
	                        const uint32_t& n_entries,
	                        const uint32_t strides_each);
    ~PrimitiveGroupContainer(void);

    inline PrimitiveGroupContainerInterface* Clone(){ return(new self_type(*this)); }

	PrimitiveGroupContainerInterface* Move(void);
	uint32_t BalanceVector(void);
	yon1_dc_t ToDataContainer(void);
	yon1_dc_t& UpdateDataContainer(yon1_dc_t& container);

	void resize(void);
	void resize(const size_t new_size);

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

	inline io::BasicBuffer& ToVcfString(io::BasicBuffer& buffer, const uint64_t position) const{
		return(this->at(position).ToVcfString(buffer));
	}

	inline bcf1_t* UpdateHtslibVcfRecordFormatInt32(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const{
		return(rec);
	}

	inline bcf1_t* UpdateHtslibVcfRecordFormatFloat(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const{
		return(rec);
	}

	bcf1_t* UpdateHtslibVcfRecordFormatString(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const;

private:
    pointer containers_;
};


// IMPLEMENTATION -------------------------------------------------------------
// templated


template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer() :
	PrimitiveContainerInterface(false, 0),
	entries_(nullptr)
{

}


template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer(const return_type value) :
	PrimitiveContainerInterface(false, 1),
	entries_(new return_type[1])
{
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
		case(YON_TYPE_8B):     (this->SetupSigned<int8_t>(container,  0)); break;
		case(YON_TYPE_16B):    (this->SetupSigned<int16_t>(container, 0)); break;
		case(YON_TYPE_32B):    (this->SetupSigned<int32_t>(container, 0)); break;
		case(YON_TYPE_64B):    (this->SetupSigned<int64_t>(container, 0)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container,  0));        break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container, 0));        break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default:
			std::cerr << utility::timestamp("ERROR") << "Disallowed: " << container.header.data_header.GetPrimitiveType() << std::endl;
			return;
		}
	} else {
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->Setup<uint8_t>(container,  0)); break;
		case(YON_TYPE_16B):    (this->Setup<uint16_t>(container, 0)); break;
		case(YON_TYPE_32B):    (this->Setup<uint32_t>(container, 0)); break;
		case(YON_TYPE_64B):    (this->Setup<uint64_t>(container, 0)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container,    0)); break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container,   0)); break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default:
			std::cerr << utility::timestamp("ERROR") << "Disallowed: " << container.header.data_header.GetPrimitiveType() << std::endl;
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
		case(YON_TYPE_8B):     (this->SetupSigned<int8_t>(container,  offset)); break;
		case(YON_TYPE_16B):    (this->SetupSigned<int16_t>(container, offset)); break;
		case(YON_TYPE_32B):    (this->SetupSigned<int32_t>(container, offset)); break;
		case(YON_TYPE_64B):    (this->SetupSigned<int64_t>(container, offset)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container,  offset));        break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container, offset));        break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default: std::cerr << utility::timestamp("ERROR") << "Disallowed: " << container.header.data_header.GetPrimitiveType() << std::endl; return;
		}
	} else {
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->Setup<uint8_t>(container,  offset)); break;
		case(YON_TYPE_16B):    (this->Setup<uint16_t>(container, offset)); break;
		case(YON_TYPE_32B):    (this->Setup<uint32_t>(container, offset)); break;
		case(YON_TYPE_64B):    (this->Setup<uint64_t>(container, offset)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container,  offset));   break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container, offset));   break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default: std::cerr << utility::timestamp("ERROR") << "Disallowed: " << container.header.data_header.GetPrimitiveType() << std::endl; return;
		}
	}
}

template <class return_type>
PrimitiveContainer<return_type>::PrimitiveContainer(return_type* values, const uint32_t n_entries) :
    PrimitiveContainerInterface(false, n_entries),
	entries_(new value_type[n_entries])
{
	for(int i = 0; i < n_entries; ++i)
		this->entries_[i] = values[i];
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
PrimitiveContainerInterface* PrimitiveContainer<return_type>::Move(void){
	self_type* o   = new self_type();
	o->is_uniform_ = this->is_uniform_;
	o->n_capacity_ = this->n_capacity_;
	o->n_entries_  = this->n_entries_;
	std::swap(this->entries_, o->entries_);
	this->n_capacity_ = 0;
	this->n_entries_  = 0;
	return(o);
}

template <class return_type>
yon1_dc_t PrimitiveContainer<return_type>::ToDataContainer(void){
	if(this->size() == 0)
		return(yon1_dc_t());

	yon1_dc_t d;
	d.data_uncompressed.resize(this->size()*sizeof(return_type) + 128);
	d.strides_uncompressed.resize(this->size()*sizeof(return_type) + 128);
	d.header.data_header.stride = this->size();

	for(uint32_t i = 0; i < this->size(); ++i) d.Add(this->at(i));
	d.AddStride(this->size());
	++d;

	return(d);
}

template <class return_type>
yon1_dc_t& PrimitiveContainer<return_type>::UpdateDataContainer(yon1_dc_t& container, const bool update_stride){
	if(this->size() == 0)
		return(container);

	if(container.data_uncompressed.size() + this->size()*sizeof(return_type) > container.data_uncompressed.capacity())
		container.data_uncompressed.resize((container.data_uncompressed.size()+this->size()*sizeof(return_type))*2);

	for(uint32_t i = 0; i < this->size(); ++i)
		container.Add(this->at(i));

	if(update_stride) container.AddStride(this->size());
	++container;

	return(container);
}

template <class return_type>
void PrimitiveContainer<return_type>::ExpandEmpty(const uint32_t to){
	// Cannot resize into mixed stride with EOV values
	// if there are no negative values available.
	assert(std::numeric_limits<return_type>::min() != 0);

	if(to > this->capacity())
		this->resize(to + 10);

	for(int i = 0; i < n_entries_; to)
		entries_[i] = std::numeric_limits<return_type>::min() + 1;
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
	if(this->n_entries_ == this->n_capacity_)
		this->resize();

	if(entry < std::numeric_limits<return_type>::min()){
		std::cerr << utility::timestamp("WARNING") << "Truncating value: " << (int)entry << "->" << (int)return_type(entry) << std::endl;
	}

	if(entry == INT8_MIN) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min();
	else if(entry == INT8_MIN+1) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min()+1;
	else this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const int16_t entry){
	if(this->n_entries_ == this->n_capacity_)
		this->resize();

	if(entry < std::numeric_limits<return_type>::min()){
		std::cerr << utility::timestamp("WARNING") << "Truncating value: " << (int)entry << "->" << (int)return_type(entry) << std::endl;
	}

	if(entry == INT16_MIN) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min();
	else if(entry == INT16_MIN+1) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min()+1;
	else this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const int32_t entry){
	if(this->n_entries_ == this->n_capacity_)
		this->resize();

	if(entry < std::numeric_limits<return_type>::min()){
		std::cerr << utility::timestamp("WARNING") << "Truncating value: " << (int)entry << "->" << (int)return_type(entry) << std::endl;
	}

	if(entry == INT32_MIN) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min();
	else if(entry == INT32_MIN+1) this->entries_[this->n_entries_++] = std::numeric_limits<return_type>::min()+1;
	else this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const int64_t entry){
	if(this->n_entries_ == this->n_capacity_)
		this->resize();

	if(entry < std::numeric_limits<return_type>::min()){
		std::cerr << utility::timestamp("WARNING") << "Truncating value: " << (int64_t)entry << "->" << (int64_t)return_type(entry) << std::endl;
	}

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const uint8_t entry){
	if(this->n_entries_ == this->n_capacity_)
		this->resize();

	if(entry > std::numeric_limits<return_type>::max()){
		std::cerr << utility::timestamp("WARNING") << "Truncating value: " << (uint32_t)entry << "->" << (uint32_t)return_type(entry) << std::endl;
	}

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const uint16_t entry){
	if(this->n_entries_ == this->n_capacity_)
		this->resize();

	if(entry > std::numeric_limits<return_type>::max()){
		std::cerr << utility::timestamp("WARNING") << "Truncating value: " << (uint32_t)entry << "->" << (uint32_t)return_type(entry) << std::endl;
	}

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const uint32_t entry){
	if(this->n_entries_ == this->n_capacity_)
		this->resize();

	if(entry > std::numeric_limits<return_type>::max()){
		std::cerr << utility::timestamp("WARNING") << "Truncating value: " << (uint32_t)entry << "->" << (uint32_t)return_type(entry) << std::endl;
	}


	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const uint64_t entry){
	if(this->n_entries_ == this->n_capacity_)
		this->resize();

	if(entry > std::numeric_limits<return_type>::max()){
		std::cerr << utility::timestamp("WARNING") << "Truncating value: " << (uint64_t)entry << "->" << (uint64_t)return_type(entry) << std::endl;
	}

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const char entry){
	if(this->n_entries_ == this->n_capacity_)
		this->resize();

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const float entry){
	if(this->n_entries_ == this->n_capacity_)
		this->resize();

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const double entry){
	if(this->n_entries_ == this->n_capacity_)
		this->resize();

	this->entries_[this->n_entries_++] = entry;
}

template <class return_type>
void PrimitiveContainer<return_type>::operator+=(const std::string& entry){
	std::cerr << utility::timestamp("WARNING") << "Cannot add a string to this container" << std::endl;
}


// IMPLEMENTATION GROUP --------------------------------------------------------
// templated


template <class return_type>
PrimitiveGroupContainer<return_type>::PrimitiveGroupContainer() : containers_(nullptr){}

template <class return_type>
PrimitiveGroupContainer<return_type>::PrimitiveGroupContainer(const self_type& other) :
	PrimitiveGroupContainerInterface(other),
	containers_(static_cast<pointer>(::operator new[](this->n_capacity_*sizeof(value_type))))
{
	for(int i = 0; i < this->size(); ++i)
		new( &this->containers_[i] ) value_type( other.at(i) );
}

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
		case(YON_TYPE_8B):     (this->Setup<int8_t>(container,  offset, n_objects, strides_each)); break;
		case(YON_TYPE_16B):    (this->Setup<int16_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_32B):    (this->Setup<int32_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_64B):    (this->Setup<int64_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container,   offset, n_objects, strides_each)); break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container,  offset, n_objects, strides_each)); break;
		case(YON_TYPE_BOOLEAN):
		case(YON_TYPE_CHAR):
		case(YON_TYPE_STRUCT):
		case(YON_TYPE_UNKNOWN):
		default: std::cerr << utility::timestamp("ERROR","PGC") << "Disallowed: " << container.header.data_header.GetPrimitiveType() << std::endl; return;
		}
	} else {
		switch(container.header.data_header.GetPrimitiveType()){
		case(YON_TYPE_8B):     (this->Setup<uint8_t>(container,  offset, n_objects, strides_each)); break;
		case(YON_TYPE_16B):    (this->Setup<uint16_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_32B):    (this->Setup<uint32_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_64B):    (this->Setup<uint64_t>(container, offset, n_objects, strides_each)); break;
		case(YON_TYPE_FLOAT):  (this->Setup<float>(container,    offset, n_objects, strides_each)); break;
		case(YON_TYPE_DOUBLE): (this->Setup<double>(container,   offset, n_objects, strides_each)); break;
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

template <class return_type>
PrimitiveGroupContainerInterface* PrimitiveGroupContainer<return_type>::Move(void){
	self_type* o = new self_type();
	o->n_capacity_ = this->n_capacity_;
	o->n_objects_  = this->n_objects_;
	this->n_capacity_ = 0;
	this->n_objects_  = 0;
	std::swap(o->containers_, this->containers_);
	return(o);
}

template <class return_type>
yon1_dc_t PrimitiveGroupContainer<return_type>::ToDataContainer(void){
	if(this->size() == 0)
		return yon1_dc_t();

	uint32_t n_entries = 0;
	uint32_t ref_size = this->at(0).size();
	for(uint32_t i = 0; i < this->size(); ++i){
		assert(this->at(i).size() == ref_size);
		n_entries += this->at(i).size();
	}

	yon1_dc_t d;
	d.data_uncompressed.resize(n_entries*sizeof(return_type) + 128);
	d.strides_uncompressed.resize(n_entries*sizeof(return_type) + 128);
	d.header.data_header.stride = ref_size;

	for(uint32_t i = 0; i < this->size(); ++i)
		this->at(i).UpdateDataContainer(d, false);

	d.AddStride(ref_size);

	return(d);
}

template <class return_type>
yon1_dc_t& PrimitiveGroupContainer<return_type>::UpdateDataContainer(yon1_dc_t& container){
	if(this->size() == 0)
		return(container);

	uint32_t n_entries = 0;
	uint32_t ref_size = this->at(0).size();
	for(uint32_t i = 0; i < this->size(); ++i){
		assert(this->at(i).size() == ref_size);
		n_entries += this->at(i).size();
	}

	if(container.data_uncompressed.size() + n_entries > container.data_uncompressed.capacity())
		container.data_uncompressed.resize((container.data_uncompressed.size()+n_entries)*2);

	for(uint32_t i = 0; i < this->size(); ++i)
		this->at(i).UpdateDataContainer(container, false);

	container.AddStride(ref_size);

	return(container);
}

template <class return_type>
void PrimitiveGroupContainer<return_type>::resize(void){
	pointer temp       = this->containers_;
	this->n_capacity_ *= 2;
	this->containers_  = static_cast<pointer>(::operator new[](this->n_capacity_*sizeof(value_type)));

	for(uint32_t i = 0; i < this->size(); ++i)
		new( &this->containers_[i] ) value_type( temp[i] );

	// Delete old data.
	for(std::size_t i = 0; i < this->size(); ++i)
		((temp + i)->~value_type)();

	::operator delete[](static_cast<void*>(temp));
}

template <class return_type>
void PrimitiveGroupContainer<return_type>::resize(const size_t new_size){
	// if new size < current capacity
	if(new_size < this->n_capacity_){
		// if new size < current number of entries
		if(new_size < this->n_objects_){
			this->n_objects_ = new_size;
			return;
		}
		return;
	}

	pointer temp       = this->containers_;
	this->n_capacity_  = new_size;
	this->containers_  = static_cast<pointer>(::operator new[](this->n_capacity_*sizeof(value_type)));

	for(uint32_t i = 0; i < this->size(); ++i)
		new( &this->containers_[i] ) value_type( temp[i] );

	// Delete old data.
	for(std::size_t i = 0; i < this->size(); ++i)
		((temp + i)->~value_type)();

	::operator delete[](static_cast<void*>(temp));
}

}

#endif /* PrimitiveContainer_H_ */
