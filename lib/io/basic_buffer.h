#ifndef BASICBUFFER__H_
#define BASICBUFFER__H_

#include <cstddef>
#include <iostream>
#include <inttypes.h>
#include <stdint.h>
#include <cassert>

#include "support/helpers.h"
#include "containers/components/generic_iterator.h"

namespace tachyon {
namespace io{

struct BasicBuffer{
public:
    typedef BasicBuffer       self_type;
    typedef char              value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef std::ptrdiff_t    difference_type;
    typedef std::size_t       size_type;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
	BasicBuffer();
	BasicBuffer(const uint64_t size);
	BasicBuffer(char* target, const size_t length);
	BasicBuffer(const uint64_t size, char* target);
	BasicBuffer(const self_type& other);
	virtual ~BasicBuffer();
	self_type& operator=(const self_type& other);
	self_type& operator=(self_type&& other) noexcept;

	inline reference back(void){ return(this->buffer_[this->n_chars_-1]); }
	inline reference front(void){ return(this->buffer_[0]); }
	inline const_reference back(void) const{ return(this->buffer_[this->n_chars_-1]); }
	inline const_reference front(void) const{ return(this->buffer_[0]); }

	// Iterator
	inline iterator begin(){ return iterator(&this->buffer_[0]); }
	inline iterator end()  { return iterator(&this->buffer_[this->n_chars_]); }
	inline const_iterator begin()  const{ return const_iterator(&this->buffer_[0]); }
	inline const_iterator end()    const{ return const_iterator(&this->buffer_[this->n_chars_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->buffer_[0]); }
	inline const_iterator cend()   const{ return const_iterator(&this->buffer_[this->n_chars_]); }

	inline void set(const size_t size){
		this->n_chars_ = 0;
		this->width_ = size;
		if(this->buffer_ != nullptr)
			delete [] this->buffer_;

		this->buffer_ = new char[size];
	}

	inline void set(const size_t size, char* target){
		this->n_chars_ = 0;
		this->width_ = size;
		this->buffer_ = target;
	}

	inline void set(char* target){
		this->n_chars_ = 0;
		this->width_ = 0;
		this->buffer_ = target;
	}

	inline void reset(){ this->n_chars_ = 0; this->iterator_position_ = 0; }
	inline void resetIterator(){ this->iterator_position_ = 0; }
	inline void move(const uint64_t to){ this->n_chars_ = to; }
	inline const uint64_t& size(void) const{ return this->n_chars_; }
	inline const uint64_t& capacity(void) const{ return this->width_; }

	void resize(const uint64_t new_size);
	void resize(const self_type& other);

	void Add(const char* data, const uint32_t length);
	void AddReadble(const int8_t& value);
	void AddReadble(const int16_t& value);
	void AddReadble(const int32_t& value);
	void AddReadble(const uint8_t& value);
	void AddReadble(const uint16_t& value);
	void AddReadble(const uint32_t& value);
	void AddReadble(const uint64_t& value);
	void AddReadble(const float& value);
	void AddReadble(const double& value);
	void AddReadble(const std::string& value);
	self_type& operator+=(const self_type& other);
	self_type& operator+=(const char& value);
	self_type& operator+=(const int8_t& value);
	self_type& operator+=(const uint8_t& value);
	self_type& operator+=(const float& value);
	self_type& operator+=(const uint16_t& value);
	self_type& operator+=(const int16_t& value);
	self_type& operator+=(const uint32_t& value);
	self_type& operator+=(const int32_t& value);
	self_type& operator+=(const double& value);
	self_type& operator+=(const uint64_t& value);
	self_type& operator+=(const int64_t& value);
	self_type& operator+=(const std::string& value);

	inline reference operator[](const uint64_t position){ return this->buffer_[position]; }
	inline const_reference operator[](const uint64_t position) const{ return this->buffer_[position]; }
	inline reference at(const uint64_t position){ return this->buffer_[position]; }
	inline const_reference at(const uint64_t position) const{ return this->buffer_[position]; }
	inline pointer data(void){ return(this->buffer_); }
	inline const_pointer data(void) const{ return(this->buffer_); }

	void read(char* target, const uint32_t n_length){
		memcpy(target, &this->buffer_[this->iterator_position_], n_length);
		this->iterator_position_ += n_length;
	}

private:
	friend self_type& operator>>(self_type& data, uint8_t& target){
		target = *reinterpret_cast<uint8_t*>(&data.buffer_[data.iterator_position_++]);
		return(data);
	}

	friend self_type& operator>>(self_type& data, uint16_t& target){
		target = *reinterpret_cast<uint16_t*>(&data.buffer_[data.iterator_position_]);
		data.iterator_position_ += sizeof(uint16_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, uint32_t& target){
		target = *reinterpret_cast<uint32_t*>(&data.buffer_[data.iterator_position_]);
		data.iterator_position_ += sizeof(uint32_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, uint64_t& target){
		target = *reinterpret_cast<uint64_t*>(&data.buffer_[data.iterator_position_]);
		data.iterator_position_ += sizeof(uint64_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, int8_t& target){
		target = *reinterpret_cast<int8_t*>(&data.buffer_[data.iterator_position_++]);
		return(data);
	}

	friend self_type& operator>>(self_type& data, int16_t& target){
		target = *reinterpret_cast<int16_t*>(&data.buffer_[data.iterator_position_]);
		data.iterator_position_ += sizeof(int16_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, int32_t& target){
		target = *reinterpret_cast<int32_t*>(&data.buffer_[data.iterator_position_]);
		data.iterator_position_ += sizeof(int32_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, int64_t& target){
		target = *reinterpret_cast<int64_t*>(&data.buffer_[data.iterator_position_]);
		data.iterator_position_ += sizeof(int64_t);
		return(data);
	}

	friend self_type& operator>>(self_type& data, float& target){
		target = *reinterpret_cast<float*>(&data.buffer_[data.iterator_position_]);
		data.iterator_position_ += sizeof(float);
		return(data);
	}

	friend self_type& operator>>(self_type& data, double& target){
		target = *reinterpret_cast<double*>(&data.buffer_[data.iterator_position_]);
		data.iterator_position_ += sizeof(double);
		return(data);
	}

	friend std::ostream& operator<<(std::ostream& out, const self_type& data){
		out.write(data.data(), data.size());
		return(out);
	}

public:
	bool         owns_data_;
	uint64_t     n_chars_;
	uint64_t     width_;
	uint64_t     iterator_position_;
	pointer      buffer_;
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

#endif /* BASICBUFFER__H_ */
