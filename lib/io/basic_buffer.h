#ifndef BASICBUFFER_H_
#define BASICBUFFER_H_

#include <cstddef>
#include <iostream>
#include <inttypes.h>
#include <stdint.h>
#include <cassert>

#include "support/helpers.h"

namespace tachyon {
namespace io{

struct BasicBuffer{
private:
    typedef BasicBuffer       self_type;
    typedef char              value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef std::ptrdiff_t    difference_type;
    typedef std::size_t       size_type;

public:
	BasicBuffer() : owns_data_(true), n_chars(0), width(0), iterator_position_(0), buffer(nullptr){}
	BasicBuffer(const uint64_t size) : owns_data_(true), n_chars(0), width(size), iterator_position_(0), buffer(new char[size]){}
	BasicBuffer(char* target, const size_t length) : owns_data_(false), n_chars(length), width(length), iterator_position_(0), buffer(target){}
	BasicBuffer(const uint64_t size, char* target) : owns_data_(false), n_chars(0), width(size), iterator_position_(0), buffer(target){}

	BasicBuffer(const self_type& other) :
		owns_data_(other.owns_data_),
		n_chars(other.n_chars),
		width(other.width),
		iterator_position_(other.iterator_position_),
		buffer(new char[other.width])
	{
		memcpy(this->buffer, other.buffer, other.size());
	}

	virtual ~BasicBuffer(){
		if(this->owns_data_)
			delete [] this->buffer;
	}

	self_type& operator=(const self_type& other){
		if(this->owns_data_) delete [] this->buffer;
		this->owns_data_ = other.owns_data_;
		this->n_chars    = other.n_chars;
		this->width      = other.width;
		this->iterator_position_ = other.iterator_position_;
		this->buffer     = new char[other.width];
		memcpy(this->buffer, other.buffer, other.size());
		return(*this);
	}

	self_type& operator=(self_type&& other) noexcept{
		if(this->owns_data_) delete [] this->buffer;
		this->buffer = nullptr;
		std::swap(this->buffer, other.buffer);
		this->owns_data_ = other.owns_data_;
		this->n_chars    = other.n_chars;
		this->width      = other.width;
		this->iterator_position_ = other.iterator_position_;
		return(*this);
	}

	inline reference back(void){ return(this->buffer[this->n_chars-1]); }
	inline reference front(void){ return(this->buffer[0]); }
	inline const_reference back(void) const{ return(this->buffer[this->n_chars-1]); }
	inline const_reference front(void) const{ return(this->buffer[0]); }

	class iterator{
	private:
		typedef iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		reference operator*() const{ return *ptr_; }
		pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	class const_iterator{
	private:
		typedef const_iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		const_iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		const_reference operator*() const{ return *ptr_; }
		const_pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	// Iterator
	inline iterator begin(){ return iterator(&this->buffer[0]); }
	inline iterator end()  { return iterator(&this->buffer[this->n_chars]); }
	inline const_iterator begin()  const{ return const_iterator(&this->buffer[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->buffer[this->n_chars]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->buffer[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->buffer[this->n_chars]); }

	inline void set(const size_t size){
		this->n_chars = 0;
		this->width = size;
		if(this->buffer != nullptr)
			delete [] this->buffer;

		this->buffer = new char[size];
	}

	inline void set(const size_t size, char* target){
		this->n_chars = 0;
		this->width = size;
		this->buffer = target;
	}

	inline virtual void set(char* target){
		this->n_chars = 0;
		this->width = 0;
		this->buffer = target;
	}

	inline void reset(){ this->n_chars = 0; this->iterator_position_ = 0; }
	inline void resetIterator(){ this->iterator_position_ = 0; }
	inline void move(const uint64_t to){ this->n_chars = to; }
	inline const uint64_t& size(void) const{ return this->n_chars; }
	inline const uint64_t& capacity(void) const{ return this->width; }

	void resize(const uint64_t new_size){
		if(this->n_chars == 0 && new_size == 0) return;
		char* temp = new char[new_size];
		//std::cerr << this->size() << "<" << new_size << std::endl;
		assert(this->size() < new_size);
		memcpy(temp, this->buffer, this->size());
		delete [] this->buffer;
		this->buffer = temp;
		this->width = new_size;
	}

	void resize(const self_type& other){
		if(other.size() >= this->capacity()){
			this->resize(other.capacity());
		}
	}

	void Add(const char* data, const uint32_t length){
		if(this->size() + length >= this->capacity())
			this->resize((this->size() + length) * 2);

		memcpy(&this->buffer[this->n_chars], &data[0], length);
		this->n_chars += length;
	}

	void AddReadble(const int8_t& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(std::max(this->width + 100, this->width*2));
		const int ret = sprintf(&this->buffer[this->n_chars], "%d", value);
		this->n_chars += ret;
	}

	void AddReadble(const int16_t& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(std::max(this->width + 100, this->width*2));
		const int ret = sprintf(&this->buffer[this->n_chars], "%d", value);
		this->n_chars += ret;
	}

	void AddReadble(const int32_t& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(std::max(this->width + 100, this->width*2));
		const int ret = sprintf(&this->buffer[this->n_chars], "%d", value);
		this->n_chars += ret;
	}

	void AddReadble(const uint8_t& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(std::max(this->width + 100, this->width*2));
		const int ret = sprintf(&this->buffer[this->n_chars], "%u", value);
		this->n_chars += ret;
	}

	void AddReadble(const uint16_t& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(std::max(this->width + 100, this->width*2));
		const int ret = sprintf(&this->buffer[this->n_chars], "%u", value);
		this->n_chars += ret;
	}

	void AddReadble(const uint32_t& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(std::max(this->width + 100, this->width*2));
		const int ret = sprintf(&this->buffer[this->n_chars], "%u", value);
		this->n_chars += ret;
	}

	void AddReadble(const uint64_t& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(std::max(this->width + 100, this->width*2));
		const int ret = sprintf(&this->buffer[this->n_chars], "%" PRIu64, value);
		this->n_chars += ret;
	}

	void AddReadble(const float& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(std::max(this->width + 100, this->width*2));
		const int ret = sprintf(&this->buffer[this->n_chars], "%g", value);
		this->n_chars += ret;
	}

	void AddReadble(const double& value){
		if(this->n_chars + 100 >= this->width)
			this->resize(std::max(this->width + 100, this->width*2));
		const int ret = sprintf(&this->buffer[this->n_chars], "%g", value);
		this->n_chars += ret;
	}

	void AddReadble(const std::string& value){
		if(this->n_chars + value.size() >= this->width)
			this->resize(std::max(this->n_chars + value.size() + 100, this->width*2));
		*this += value;
	}

	inline self_type& operator+=(const self_type& other){
		if(this->size() + other.size() >= this->capacity())
			this->resize((this->size() + other.size()) * 2);

		memcpy(&this->buffer[this->n_chars], other.buffer, other.n_chars);
		this->n_chars += other.n_chars;

		return *this;
	}

	inline self_type& operator+=(const char& value){
		if(this->n_chars + sizeof(char) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		this->buffer[this->n_chars] = value;
		++this->n_chars;
		return *this;
	}

	inline self_type& operator+=(const int8_t& value){
		if(this->n_chars + sizeof(int8_t) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		this->buffer[this->n_chars] = value;
		++this->n_chars;
		return *this;
	}

	inline self_type& operator+=(const uint8_t& value){
		if(this->n_chars + sizeof(uint8_t) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		uint8_t* p = reinterpret_cast<uint8_t*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(uint8_t);
		return *this;
	}

	inline self_type& operator+=(const float& value){
		if(this->n_chars + sizeof(float) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		float* p = reinterpret_cast<float*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(float);
		return *this;
	}

	inline self_type& operator+=(const uint16_t& value){
		if(this->n_chars + sizeof(uint16_t) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		uint16_t* p = reinterpret_cast<uint16_t*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(uint16_t);
		return *this;
	}

	inline self_type& operator+=(const int16_t& value){
		if(this->n_chars + sizeof(int16_t) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		short* p = reinterpret_cast<short*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(short);
		return *this;
	}

	inline self_type& operator+=(const uint32_t& value){
		if(this->n_chars + sizeof(uint32_t) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		uint32_t* p = reinterpret_cast<uint32_t*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(uint32_t);
		return *this;
	}

	inline self_type& operator+=(const int32_t& value){
		if(this->n_chars + sizeof(int32_t) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		int32_t* p = reinterpret_cast<int32_t*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(int32_t);
		return *this;
	}

	inline self_type& operator+=(const double& value){
		if(this->n_chars + sizeof(double) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		double* p = reinterpret_cast<double*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(double);
		return *this;
	}

	inline self_type& operator+=(const uint64_t& value){
		if(this->n_chars + sizeof(uint64_t) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		uint64_t* p = reinterpret_cast<uint64_t*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(uint64_t);
		return *this;
	}

	inline self_type& operator+=(const int64_t& value){
		if(this->n_chars + sizeof(int64_t) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		int64_t* p = reinterpret_cast<int64_t*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(int64_t);
		return *this;
	}

	/*
	inline self_type& operator+=(const size_t& value){
		if(this->n_chars + sizeof(size_t) >= this->width)
			this->resize(std::max(this->width + 1000, this->width*2));

		size_t* p = reinterpret_cast<size_t*>(&this->buffer[this->n_chars]);
		*p = value;
		this->n_chars += sizeof(size_t);
		return *this;
	}
	*/

	inline self_type& operator+=(const std::string& value){
		if(this->n_chars + value.size() + sizeof(uint8_t) >= this->width){
			uint64_t resize_to = std::max(this->n_chars + value.size() + sizeof(uint8_t) + 1000, this->width * 2);
			this->resize(resize_to);
		}

		for(uint32_t i = 0; i < value.size(); ++i){
			this->buffer[this->n_chars] = value[i];
			++this->n_chars;
		}

		return *this;
	}

	inline reference operator[](const uint64_t position){ return this->buffer[position]; }
	inline const_reference operator[](const uint64_t position) const{ return this->buffer[position]; }
	inline reference at(const uint64_t position){ return this->buffer[position]; }
	inline const_reference at(const uint64_t position) const{ return this->buffer[position]; }
	inline pointer data(void){ return(this->buffer); }
	inline const_pointer data(void) const{ return(this->buffer); }

	void read(char* target, const uint32_t n_length){
		memcpy(target, &this->buffer[this->iterator_position_], n_length);
		this->iterator_position_ += n_length;
	}

private:
	friend self_type& operator>>(self_type& data, uint8_t& target){
		target = *reinterpret_cast<uint8_t*>(&data.buffer[data.iterator_position_++]);
		return(data);
	}

	friend self_type& operator>>(self_type& data, uint16_t& target){
		target = *reinterpret_cast<uint16_t*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(uint16_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, uint32_t& target){
		target = *reinterpret_cast<uint32_t*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(uint32_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, uint64_t& target){
		target = *reinterpret_cast<uint64_t*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(uint64_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, int8_t& target){
		target = *reinterpret_cast<int8_t*>(&data.buffer[data.iterator_position_++]);
		return(data);
	}

	friend self_type& operator>>(self_type& data, int16_t& target){
		target = *reinterpret_cast<int16_t*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(int16_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, int32_t& target){
		target = *reinterpret_cast<int32_t*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(int32_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, int64_t& target){
		target = *reinterpret_cast<int64_t*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(int64_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, float& target){
		target = *reinterpret_cast<float*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(float);
		return(data);
	}

	friend self_type& operator>>(self_type& data, double& target){
		target = *reinterpret_cast<double*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(double);
		return(data);
	}

	/*
	friend self_type& operator>>(self_type& data, size_type& target){
		target = *reinterpret_cast<size_type*>(&data.buffer[data.iterator_position_]);
		data.iterator_position_ += sizeof(size_type);
		return(data);
	}
	*/

	friend std::ostream& operator<<(std::ostream& out, const self_type& data){
		out.write(data.data(), data.size());
		return(out);
	}

public:
	bool    owns_data_;
	uint64_t     n_chars;
	uint64_t     width;
	uint64_t     iterator_position_;
	pointer buffer;
};

void SerializeString(const std::string& string, io::BasicBuffer& buffer);
void DeserializeString(std::string& string, io::BasicBuffer& buffer);

template <class T>
static void SerializePrimitive(const T& value, io::BasicBuffer& buffer){
	buffer += value;
}

template <class T>
static void DeserializePrimitive(T& value, io::BasicBuffer& buffer){
	buffer >> value;
}

} /* namespace IO */
} /* namespace Tomahawk */

#endif /* BASICBUFFER_H_ */
