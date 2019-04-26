#include "primitive_container.h"

namespace tachyon {

PrimitiveContainerInterface::PrimitiveContainerInterface(void) :
	is_uniform_(false),
	n_entries_(0),
	n_capacity_(0)
{

}

PrimitiveContainerInterface::PrimitiveContainerInterface(const bool uniform, const size_t size) :
	is_uniform_(uniform),
	n_entries_(size),
	n_capacity_(size)
{

}

PrimitiveContainerInterface::~PrimitiveContainerInterface() {}

PrimitiveContainer<std::string>::PrimitiveContainer(void) :
	PrimitiveContainerInterface(false, 1)
{
	this->n_capacity_ = 1;
	this->n_entries_  = 0;
}

PrimitiveContainer<std::string>::PrimitiveContainer(const char* data, const size_t l_data) :
	PrimitiveContainerInterface(false, 1),
	data_(data, l_data)
{

}

PrimitiveContainer<std::string>::PrimitiveContainer(const PrimitiveContainer& other) :
	PrimitiveContainerInterface(other.is_uniform_, 1),
	data_(other.data_.data(), other.data_.size())
{
}

PrimitiveContainer<std::string>& PrimitiveContainer<std::string>::operator=(const PrimitiveContainer& other) {
	this->is_uniform_ = other.is_uniform_;
	this->n_entries_ = other.n_entries_;
	data_ = other.data_;
	return(*this);
}

PrimitiveContainer<std::string>::~PrimitiveContainer(void) {}

PrimitiveContainerInterface* PrimitiveContainer<std::string>::Move(void) {
	self_type* o      = new self_type();
	o->is_uniform_ = this->is_uniform_;
	o->n_capacity_ = this->n_capacity_;
	o->n_entries_  = this->n_entries_;
	this->data_ = std::move(o->data_);
	this->data_.clear();
	this->n_entries_  = 0;
	this->n_capacity_ = 1;
	return(o);
}

void PrimitiveContainer<std::string>::ExpandEmpty(const uint32_t to) {
	const size_t sz = data_.size();
	for (int i = sz; i < to; ++i) data_.push_back('\0');
}

yon1_dc_t PrimitiveContainer<std::string>::ToDataContainer(void) {
	if (this->size() == 0)
		return yon1_dc_t();

	yon1_dc_t d;
	d.data_uncompressed.resize(this->size() + 128);
	d.strides_uncompressed.resize(this->size() + 128);
	d.header.data_header.stride = this->size();

	for (uint32_t i = 0; i < this->size(); ++i) {
		d.AddString(this->data_.c_str(), this->data_.size());
		d.AddStride(this->data_.size());
		++d;
	}

	return(d);
}

yon1_dc_t& PrimitiveContainer<std::string>::UpdateDataContainer(yon1_dc_t& container, const bool update_stride) {
	if (this->size() == 0)
		return container;

	if (container.data_uncompressed.size() + this->size() > container.data_uncompressed.capacity())
		container.data_uncompressed.resize((container.data_uncompressed.size()+this->size())*2);

	for (uint32_t i = 0; i < this->size(); ++i) {
		// Add C-style strings to buffer.
		container.AddString(this->data_.c_str(), this->data_.size());
		if (update_stride) container.AddStride(this->data_.size());
		++container;
	}

	return(container);
}

yon_buffer_t& PrimitiveContainer<std::string>::ToVcfString(yon_buffer_t& buffer) const {
	if (this->data_.size() == 0) {
		buffer += '.';
		return(buffer);
	}
	buffer += this->data_;
	return(buffer);
}

PrimitiveGroupContainerInterface::PrimitiveGroupContainerInterface() : n_capacity_(0), n_objects_(0) { }
PrimitiveGroupContainerInterface::PrimitiveGroupContainerInterface(const size_type n_objects) : n_capacity_(n_objects), n_objects_(n_objects) { }
PrimitiveGroupContainerInterface::~PrimitiveGroupContainerInterface() {}

PrimitiveGroupContainer<std::string>::PrimitiveGroupContainer() : containers_(nullptr) {}

PrimitiveGroupContainer<std::string>::PrimitiveGroupContainer(const self_type& other) :
	PrimitiveGroupContainerInterface(other),
	containers_(static_cast<pointer>(::operator new[](this->n_capacity_*sizeof(value_type))))
{
	for (int i = 0; i < this->size(); ++i)
		new( &this->containers_[i] ) value_type( other.at(i) );
}

PrimitiveGroupContainer<std::string>::PrimitiveGroupContainer(const data_container_type& container,
	                                                          const uint32_t& offset,
	                                                          const uint32_t& n_entries,
	                                                          const uint32_t strides_each) :
	PrimitiveGroupContainerInterface(n_entries), // limitation
	containers_(static_cast<pointer>(::operator new[](this->size()*sizeof(value_type))))
{
	uint32_t current_offset = offset;
	for (size_type i = 0; i < this->size(); ++i) {
		// check length
		size_type j = 0;
		for (; j < strides_each; ++j) {
			// Find premature end-of-string marker and stop.
			if (container.data_uncompressed[current_offset + j] == '\0') {
				break;
			}
		}
		new( &this->containers_[i] ) value_type( &container.data_uncompressed[current_offset], j );
		current_offset += strides_each;
	}
}

PrimitiveGroupContainer<std::string>::~PrimitiveGroupContainer() {
	for (std::size_t i = 0; i < this->size(); ++i)
		((this->containers_ + i)->~PrimitiveContainer)();

	::operator delete[](static_cast<void*>(this->containers_));
}

void PrimitiveGroupContainer<std::string>::resize(void) {
	pointer temp       = this->containers_;
	this->n_capacity_ *= 2;
	this->containers_  = static_cast<pointer>(::operator new[](this->n_capacity_*sizeof(value_type)));

	for (uint32_t i = 0; i < this->size(); ++i)
		new( &this->containers_[i] ) value_type( temp[i] );

	// Delete old data.
	for (std::size_t i = 0; i < this->size(); ++i)
		((temp + i)->~PrimitiveContainer)();

	::operator delete[](static_cast<void*>(temp));
}

void PrimitiveGroupContainer<std::string>::resize(const size_t new_size) {
	// if new size < current capacity
	if (new_size < this->n_capacity_) {
		// if new size < current number of entries
		if (new_size < this->n_objects_) {
			this->n_objects_ = new_size;
			return;
		}
		return;
	}

	pointer temp       = this->containers_;
	this->n_capacity_  = new_size;
	this->containers_  = static_cast<pointer>(::operator new[](this->n_capacity_*sizeof(value_type)));

	for (uint32_t i = 0; i < this->size(); ++i)
		new( &this->containers_[i] ) value_type( temp[i] );

	// Delete old data.
	for (std::size_t i = 0; i < this->size(); ++i)
		((temp + i)->~PrimitiveContainer)();

	::operator delete[](static_cast<void*>(temp));
}

PrimitiveGroupContainerInterface* PrimitiveGroupContainer<std::string>::Move(void) {
	self_type* o = new self_type();
	o->n_capacity_ = this->n_capacity_;
	o->n_objects_  = this->n_objects_;
	this->n_capacity_ = 0;
	this->n_objects_  = 0;
	std::swap(o->containers_, this->containers_);
	return(o);
}

uint32_t PrimitiveGroupContainer<std::string>::BalanceVector() {
	if (this->size() == 0) return 0;

	uint32_t max_stride = 0;
	for (int i = 0; i < this->size(); ++i) {
		if (this->at(i).data_.size() > max_stride)
			max_stride = this->at(i).data_.size();
	}

	return max_stride;
}

yon1_dc_t PrimitiveGroupContainer<std::string>::ToDataContainer(void) {
	if (this->size() == 0)
		return yon1_dc_t();

	// Find the longest string in this Format record.
	const uint32_t stride    = this->BalanceVector();
	const uint32_t l_entries = this->size() * l_entries;

	yon1_dc_t d;
	d.data_uncompressed.resize(l_entries + 128);
	d.strides_uncompressed.resize(l_entries + 128);
	d.header.data_header.stride = stride;

	for (uint32_t i = 0; i < this->size(); ++i) {
		// If the target string is short then the largest string
		// available we pad the data added with null-terminations.
		// This will force the importing algorithm to skip premature
		// terminations while still maintaining the required matrix
		// properties of a Format record.
		this->at(i).ExpandEmpty(stride);
		// Add actual data.
		this->at(i).UpdateDataContainer(d, false);
	}

	d.AddStride(stride);

	return(d);
}

yon1_dc_t& PrimitiveGroupContainer<std::string>::UpdateDataContainer(yon1_dc_t& container) {
	if (this->size() == 0)
		return container;

	uint32_t l_entries = 0;
	uint32_t ref_size = this->at(0).data_.size();
	for (uint32_t i = 0; i < this->size(); ++i) {
		assert(this->at(i).data_.size() == ref_size);
		l_entries += this->at(i).data_.size();
	}

	for (uint32_t i = 0; i < this->size(); ++i)
		this->at(i).UpdateDataContainer(container, false);

	container.AddStride(ref_size);

	return(container);
}

bcf1_t* PrimitiveGroupContainer<std::string>::UpdateHtslibVcfRecordFormatString(bcf1_t* rec, bcf_hdr_t* hdr, const std::string& tag) const {
	const char** dst = new const char*[this->size()];

	for (uint32_t i = 0; i < this->size(); ++i) {
		dst[i] = this->at(i).data_.data();
		//std::cerr << dst[i] << std::endl;
	}

	bcf_update_format_string(hdr, rec, tag.data(), dst, this->size());

	delete [] dst;
	return(rec);
}

}
