#include "basic_buffer.h"

namespace tachyon{
namespace io{

BasicBuffer::BasicBuffer() :
	owns_data_(true),
	n_chars_(0),
	width_(0),
	iterator_position_(0),
	buffer_(nullptr)
{}

BasicBuffer::BasicBuffer(const uint64_t size) :
	owns_data_(true),
	n_chars_(0),
	width_(size),
	iterator_position_(0),
	buffer_(new char[size])
{}

BasicBuffer::BasicBuffer(char* target, const size_t length) :
	owns_data_(false),
	n_chars_(length),
	width_(length),
	iterator_position_(0),
	buffer_(target)
{}

BasicBuffer::BasicBuffer(const uint64_t size, char* target) :
	owns_data_(false),
	n_chars_(0),
	width_(size),
	iterator_position_(0),
	buffer_(target)
{}

BasicBuffer::BasicBuffer(const self_type& other) :
	owns_data_(other.owns_data_),
	n_chars_(other.n_chars_),
	width_(other.width_),
	iterator_position_(other.iterator_position_),
	buffer_(new char[other.width_])
{
	memcpy(this->buffer_, other.buffer_, other.size());
}

BasicBuffer::~BasicBuffer(){
	if(this->owns_data_)
		delete [] this->buffer_;
}

BasicBuffer& BasicBuffer::operator=(const self_type& other){
	if(this->owns_data_) delete [] this->buffer_;
	this->owns_data_  = other.owns_data_;
	this->n_chars_    = other.n_chars_;
	this->width_      = other.width_;
	this->iterator_position_ = other.iterator_position_;
	this->buffer_     = new char[other.width_];
	memcpy(this->buffer_, other.buffer_, other.size());
	return(*this);
}

BasicBuffer& BasicBuffer::operator=(self_type&& other) noexcept{
	if(this->owns_data_) delete [] this->buffer_;
	this->buffer_ = nullptr;
	this->owns_data_  = other.owns_data_;
	this->n_chars_    = other.n_chars_;
	this->width_      = other.width_;
	this->iterator_position_ = other.iterator_position_;
	std::swap(this->buffer_, other.buffer_);
	other.n_chars_ = 0;
	other.width_   = 0;
	other.iterator_position_ = 0;
	return(*this);
}

void BasicBuffer::resize(const uint64_t new_size){
	if(this->n_chars_ == 0 && new_size == 0) return;
	if(new_size < this->capacity()){
		if(this->n_chars_ > new_size)
			this->n_chars_ = new_size;
		return;
	}

	char* temp = new char[new_size];
	assert(this->size() < new_size);
	memcpy(temp, this->buffer_, this->size());
	delete [] this->buffer_;
	this->buffer_ = temp;
	this->width_  = new_size;
}

void BasicBuffer::resize(const self_type& other){
	if(other.size() >= this->capacity()){
		this->resize(other.capacity());
	}
}

void BasicBuffer::Add(const char* data, const uint32_t length){
	if(this->size() + length >= this->capacity())
		this->resize((this->size() + length) * 2);

	memcpy(&this->buffer_[this->n_chars_], &data[0], length);
	this->n_chars_ += length;
}

void BasicBuffer::AddReadble(const int8_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%d", value);
	this->n_chars_ += ret;
}

void BasicBuffer::AddReadble(const int16_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%d", value);
	this->n_chars_ += ret;
}

void BasicBuffer::AddReadble(const int32_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%d", value);
	this->n_chars_ += ret;
}

void BasicBuffer::AddReadble(const uint8_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%u", value);
	this->n_chars_ += ret;
}

void BasicBuffer::AddReadble(const uint16_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%u", value);
	this->n_chars_ += ret;
}

void BasicBuffer::AddReadble(const uint32_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%u", value);
	this->n_chars_ += ret;
}

void BasicBuffer::AddReadble(const uint64_t& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%" PRIu64, value);
	this->n_chars_ += ret;
}

void BasicBuffer::AddReadble(const float& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%g", value);
	this->n_chars_ += ret;
}

void BasicBuffer::AddReadble(const double& value){
	if(this->n_chars_ + 100 >= this->width_)
		this->resize(std::max(this->width_ + 100, this->width_*2));
	const int ret = sprintf(&this->buffer_[this->n_chars_], "%g", value);
	this->n_chars_ += ret;
}

void BasicBuffer::AddReadble(const std::string& value){
	if(this->n_chars_ + value.size() >= this->width_)
		this->resize(std::max(this->n_chars_ + value.size() + 100, this->width_*2));
	*this += value;
}

BasicBuffer& BasicBuffer::operator+=(const self_type& other){
	if(this->size() + other.size() >= this->capacity())
		this->resize((this->size() + other.size()) * 2);

	memcpy(&this->buffer_[this->n_chars_], other.buffer_, other.n_chars_);
	this->n_chars_ += other.n_chars_;

	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const char& value){
	if(this->n_chars_ + sizeof(char) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	this->buffer_[this->n_chars_] = value;
	++this->n_chars_;
	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const int8_t& value){
	if(this->n_chars_ + sizeof(int8_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	this->buffer_[this->n_chars_] = value;
	++this->n_chars_;
	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const uint8_t& value){
	if(this->n_chars_ + sizeof(uint8_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	uint8_t* p = reinterpret_cast<uint8_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(uint8_t);
	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const float& value){
	if(this->n_chars_ + sizeof(float) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	float* p = reinterpret_cast<float*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(float);
	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const uint16_t& value){
	if(this->n_chars_ + sizeof(uint16_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	uint16_t* p = reinterpret_cast<uint16_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(uint16_t);
	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const int16_t& value){
	if(this->n_chars_ + sizeof(int16_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	short* p = reinterpret_cast<short*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(short);
	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const uint32_t& value){
	if(this->n_chars_ + sizeof(uint32_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	uint32_t* p = reinterpret_cast<uint32_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(uint32_t);
	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const int32_t& value){
	if(this->n_chars_ + sizeof(int32_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	int32_t* p = reinterpret_cast<int32_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(int32_t);
	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const double& value){
	if(this->n_chars_ + sizeof(double) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	double* p = reinterpret_cast<double*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(double);
	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const uint64_t& value){
	if(this->n_chars_ + sizeof(uint64_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	uint64_t* p = reinterpret_cast<uint64_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(uint64_t);
	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const int64_t& value){
	if(this->n_chars_ + sizeof(int64_t) >= this->width_)
		this->resize(std::max(this->width_ + 1000, this->width_*2));

	int64_t* p = reinterpret_cast<int64_t*>(&this->buffer_[this->n_chars_]);
	*p = value;
	this->n_chars_ += sizeof(int64_t);
	return *this;
}

BasicBuffer& BasicBuffer::operator+=(const std::string& value){
	if(this->n_chars_ + value.size() + sizeof(uint8_t) >= this->width_){
		uint64_t resize_to = std::max(this->n_chars_ + value.size() + sizeof(uint8_t) + 1000, this->width_ * 2);
		this->resize(resize_to);
	}

	for(uint32_t i = 0; i < value.size(); ++i){
		this->buffer_[this->n_chars_] = value[i];
		++this->n_chars_;
	}

	return *this;
}

void SerializeString(const std::string& string, io::BasicBuffer& buffer){
	uint32_t size_helper = string.size();
	buffer += size_helper;
	buffer += string;
}

void DeserializeString(std::string& string, io::BasicBuffer& buffer){
	uint32_t size_helper;
	buffer >> size_helper;
	string.resize(size_helper);
	buffer.read(&string[0], size_helper);
}

}
}
