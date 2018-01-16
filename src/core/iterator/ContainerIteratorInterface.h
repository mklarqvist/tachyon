#ifndef CONTAINERITERATORINTERFACE_H_
#define CONTAINERITERATORINTERFACE_H_

#include <cassert>

namespace Tachyon{
namespace Iterator{

#define BCF_BYTE_MISSING  0x80
#define BCF_BYTE_EOV      0x81
#define BCF_SHORT_MISSING 0x8000
#define BCF_SHORT_EOV     0x8001
#define BCF_INT_MISSING   0x80000000
#define BCF_INT_EOV       0x80000001

/**
 * Interface for iterators over arbitrary containers
 *
 * The interface has to provide the callable functions
 * used in the derived classes. Since we have no way of
 * knowing what the return type will be a priori we
 * provide templated set/get functions that are called
 * from the derived class.
 *
 * In order to provide this abstract functionality we
 * need to know the byte-width of the return type and
 * what (if any) missing and end-of-vector (EOV) values
 * we should interpret. This functionality is provided by
 * calling the mandatory setup() function. This is
 * functionality is reminiscent of a factory pattern
 * but implicitly overloads itself rather than returning
 * a new class instance.
 */
class ContainerIteratorDataInterface{
protected:
	typedef IO::BasicBuffer buffer_type;

private:
	typedef ContainerIteratorDataInterface self_type;

	// Function pointers
	typedef void (self_type::*toStringFunctionDefinition)(std::ostream& stream, const U32& stride) const;
	typedef void (self_type::*toStringBufferFunctionDefinition)(buffer_type& buffer, const U32& stride) const;

public:
	ContainerIteratorDataInterface(const buffer_type& buffer) :
		position(0),
		n_entries(0),
		type_size(1),
		__missing_value(BCF_BYTE_MISSING),
		__end_of_vector_value(BCF_BYTE_EOV),
		__source_sign(0),
		__source_type(Core::TACHYON_CORE_TYPE(0)),
		toStringFunction(nullptr),
		toStringBufferFunction(nullptr),
		buffer(buffer)
	{

	}

	virtual ~ContainerIteratorDataInterface(){}

	inline const S32& size(void) const{ return(this->n_entries); }
	inline const U64& size_data(void) const{ return(this->buffer.size()); }

	/**<
	 * Returns a value of type T at the current position
	 * @return
	 */
	template <class T>
	inline const T& current(void) const{ return(*reinterpret_cast<const T* const>(&this->buffer.data[this->position*this->type_size])); }

	template <class T>
	inline const T* const currentAt(void) const{ return(reinterpret_cast<const T* const>(&this->buffer.data[this->position*this->type_size])); }

	template <class T>
	inline const T& first(void) const{ return(*reinterpret_cast<const T* const>(&this->buffer.data[0])); }

	template <class T>
	inline const T& last(void) const{
		if(this->n_entries - 1 < 0) return(this->first<T>());
		return(*reinterpret_cast<const T* const>(&this->buffer.data[(this->n_entries - 1)*this->type_size]));
	}

	template <class T>
	inline const T& operator[](const U32& p) const{ return(*reinterpret_cast<const T* const>(&this->buffer.data[p*this->type_size])); }

	template <class T>
	inline const T* at(const U32& p) const{ return( reinterpret_cast<const T*>(&this->buffer.data[p*this->type_size])); }

	void currentPointer(const void*& p) const{ p = reinterpret_cast<const void*>(&this->buffer.data[this->position*this->type_size]); }

	// Dangerous functions
	virtual const U32 getCurrentStride(void) const =0;

	template <class T>
	std::vector<T> getVector(void) const;


	/**<
	 *
	 * @param type
	 * @param signedness
	 */
	void setType(const U16& type, const BYTE signedness){
		if(type == Core::YON_TYPE_BOOLEAN){
			this->type_size = 0;
		} else if(type == Core::YON_TYPE_CHAR){
			this->type_size              = sizeof(char);
			this->__missing_value        = BCF_BYTE_MISSING;
			this->__end_of_vector_value  = BCF_BYTE_EOV;
			this->toStringFunction       = &self_type::__toStringNoSeparator<char>;
			this->toStringBufferFunction = &self_type::__toStringNoSeparator<char>;
		} else if(type == Core::YON_TYPE_8B){
			this->type_size = sizeof(BYTE);
			this->__missing_value = BCF_BYTE_MISSING;
			this->__end_of_vector_value = BCF_BYTE_EOV;
			if(signedness){
				this->toStringFunction = &self_type::__toStringSignedSmall<SBYTE, BYTE>;
				this->toStringBufferFunction = &self_type::__toStringSignedSmall<SBYTE, BYTE>;
			}
			else {
				this->toStringFunction = &self_type::__toStringUnsignedSmall<BYTE>;
				this->toStringBufferFunction = &self_type::__toStringUnsignedSmall<BYTE>;
			}
		} else if(type == Core::YON_TYPE_16B){
			this->type_size = sizeof(U16);
			this->__missing_value = BCF_SHORT_MISSING;
			this->__end_of_vector_value = BCF_SHORT_EOV;
			if(signedness){
				this->toStringFunction = &self_type::__toStringSigned<S16, U16>;
				this->toStringBufferFunction = &self_type::__toStringSigned<S16, U16>;
			}
			else {
				this->toStringFunction = &self_type::__toStringUnsigned<U16>;
				this->toStringBufferFunction = &self_type::__toStringUnsigned<U16>;
			}
		} else if(type == Core::YON_TYPE_32B){
			this->type_size = sizeof(U32);
			this->__missing_value = BCF_INT_MISSING;
			this->__end_of_vector_value = BCF_INT_EOV;
			if(signedness){
				this->toStringFunction = &self_type::__toStringSigned<S32, U32>;
				this->toStringBufferFunction = &self_type::__toStringSigned<S32, U32>;
			}
			else {
				this->toStringFunction = &self_type::__toStringUnsigned<U32>;
				this->toStringBufferFunction = &self_type::__toStringUnsigned<U32>;
			}
		} else if(type == Core::YON_TYPE_FLOAT){
			this->type_size = sizeof(float);
			this->__missing_value = 0;
			this->__end_of_vector_value = 0;
			this->toStringFunction = &self_type::__toStringFloat<float>;
			this->toStringBufferFunction = &self_type::__toStringFloat<float>;
		} else if(type == Core::YON_TYPE_DOUBLE){
			this->type_size = sizeof(double);
			this->__missing_value = 0;
			this->__end_of_vector_value = 0;
			this->toStringFunction = &self_type::__toStringFloat<double>;
			this->toStringBufferFunction = &self_type::__toStringFloat<double>;
		} else if(type == Core::YON_TYPE_64B){
			this->type_size = sizeof(U64);
			this->__missing_value = 0;
			this->__end_of_vector_value = 0;
			this->toStringFunction = &self_type::__toStringUnsigned<U64>;
			this->toStringBufferFunction = &self_type::__toStringUnsigned<U64>;
		} else {
				std::cerr << Helpers::timestamp("ERROR") << std::endl;
				exit(1);
		}

		if(type == Core::YON_TYPE_BOOLEAN){
			this->n_entries = 0;
		} else {
			this->n_entries = buffer.pointer / this->type_size;
			assert(buffer.pointer % this->type_size == 0);
		}

		this->__source_type = Core::TACHYON_CORE_TYPE(type);
		this->__source_sign = signedness;
	}

	// Todo:
	// Comparison operators
	// Must be overloaded
	virtual const bool operator>(const U32& cmp)  const{ return false; }
	virtual const bool operator<(const U32& cmp)  const{ return false; }
	virtual const bool operator==(const U32& cmp) const{ return false; }
	virtual const bool operator>=(const U32& cmp) const{ return false; }
	virtual const bool operator<=(const U32& cmp) const{ return false; }

	// Increment/decrement operators
	inline void operator++(){
		if(this->position + 1 == this->n_entries) return;
		++this->position;
	}

	inline void operator--(){
		if(this->position == 0) return;
		--this->position;
	}

	inline void operator+=(const U32& p){
		if(this->position + p >= this->n_entries){
			//this->position = this->n_entries - 1;
			return;
		}

		this->position += p;
	}

	inline void operator-=(const U32& p){
		if(this->position - p < 0){
			this->position = 0;
			return;
		}
		this->position -= p;
	}

	/**<
	 *
	 * @param stream
	 * @param stride
	 */
	inline void toString(std::ostream& stream, const U32& stride) const{ (this->*toStringFunction)(stream, stride); }
	inline void toString(buffer_type& buffer, const U32& stride) const{ (this->*toStringBufferFunction)(buffer, stride); }

private:
	/* ****************************************
	*  Internal convertion functions to VCF strings
	******************************************/
	/**<
	 *
	 * @param stream
	 * @param stride
	 */
	template <class T>
	void __toStringNoSeparator(std::ostream& stream, const U32& stride) const{
		if(stride == 1){
			stream.put(this->current<T>());
		} else {
			stream.write(this->currentAt<T>(), stride);
		}
	}

	template <class T>
	void __toStringNoSeparator(buffer_type& buffer, const U32& stride) const{
		if(stride == 1){
			buffer += this->current<T>();
		} else {
			buffer.Add(this->currentAt<T>(), stride);
		}
	}

	/**<
	 *
	 * @param stream
	 * @param stride
	 */
	template <class T>
	void __toStringUnsigned(std::ostream& stream, const U32& stride) const{
		if(stride == 1){
			stream << this->current<T>();
		} else {
			const T* const r = this->currentAt<T>();

			stream << r[0];
			for(U32 i = 1; i < stride; ++i){
				stream.put(',');
				stream << r[i];
			}
		}
	}

	template <class T>
	void __toStringUnsigned(buffer_type& buffer, const U32& stride) const{
		if(stride == 1){
			buffer.AddReadble(this->current<T>());
		} else {
			const T* const r = this->currentAt<T>();

			buffer.AddReadble(r[0]);
			for(U32 i = 1; i < stride; ++i){
				buffer += ',';
				buffer.AddReadble(r[i]);
			}
		}
	}

	/**<
	 *
	 * @param stream
	 * @param stride
	 */
	template <class T, class Y>
	void __toStringSigned(std::ostream& stream, const U32& stride) const{
		if(*(const Y* const)this->currentAt<T>() == this->__end_of_vector_value){
			stream.put('.');
			return;
		}

		if(stride == 1){
			if(*(const Y* const)this->currentAt<T>() == this->__missing_value){
				stream.put('.');
				return;
			}

			stream << this->current<T>();
		} else {
			const T* const r = this->currentAt<T>();
			const Y* const u = reinterpret_cast<const Y* const>(this->currentAt<T>());

			if(u[0] == this->__missing_value) stream.put('.');
			else if(u[0] == this->__end_of_vector_value) return;
			else stream << r[0];

			for(U32 i = 1; i < stride; ++i){
				if(u[i] == this->__missing_value) stream << ",.";
				else if(u[i] == this->__end_of_vector_value) return;
				else stream << ',' << r[i];
			}
		}
	}

	template <class T, class Y>
	void __toStringSigned(buffer_type& buffer, const U32& stride) const{
		if(*(const Y* const)this->currentAt<T>() == this->__end_of_vector_value){
			buffer += '.';
			return;
		}

		if(stride == 1){
			if(*(const Y* const)this->currentAt<T>() == this->__missing_value){
				buffer += '.';
				return;
			}

			buffer.AddReadble(this->current<T>());
		} else {
			const T* const r = this->currentAt<T>();
			const Y* const u = reinterpret_cast<const Y* const>(this->currentAt<T>());

			if(u[0] == this->__missing_value) buffer += '.';
			else if(u[0] == this->__end_of_vector_value) return;
			else buffer.AddReadble(r[0]);

			for(U32 i = 1; i < stride; ++i){
				if(u[i] == this->__missing_value) buffer += ",.";
				else if(u[i] == this->__end_of_vector_value) return;
				else {
					buffer += ',';
					buffer.AddReadble(r[i]);
				}
			}
		}
	}

	/**<
	 *
	 * @param stream
	 * @param stride
	 */
	template <class T>
	void __toStringUnsignedSmall(std::ostream& stream, const U32& stride) const{
		if(stride == 1){
			stream << (U32)this->current<T>();
		} else {
			const T* const r = this->currentAt<T>();

			stream << (U32)r[0];
			for(U32 i = 1; i < stride; ++i){
				stream.put(',');
				stream << (U32)r[i];
			}
		}
	}

	template <class T>
	void __toStringUnsignedSmall(buffer_type& buffer, const U32& stride) const{
		if(stride == 1){
			buffer.AddReadble(this->current<T>());
		} else {
			const T* const r = this->currentAt<T>();

			buffer.AddReadble(r[0]);
			for(U32 i = 1; i < stride; ++i){
				buffer += ',';
				buffer.AddReadble(r[i]);
			}
		}
	}

	/**<
	 *
	 * @param stream
	 * @param stride
	 */
	template <class T, class Y>
	void __toStringSignedSmall(std::ostream& stream, const U32& stride) const{
		if(*(const Y* const)this->currentAt<T>() == this->__end_of_vector_value){
			stream.put('.');
			return;
		}

		if(stride == 1){
			if(*(const Y* const)this->currentAt<T>() == this->__missing_value){
				stream.put('.');
				return;
			}

			stream << (S32)this->current<T>();
		} else {
			const T* const r = this->currentAt<T>();
			const Y* const u = reinterpret_cast<const Y* const>(this->currentAt<T>());


			if(u[0] == this->__missing_value) stream.put('.');
			else if(u[0] == this->__end_of_vector_value) return;
			else stream << (S32)r[0];

			for(U32 i = 1; i < stride; ++i){
				if(u[i] == this->__missing_value) stream << ",.";
				else if(u[i] == this->__end_of_vector_value) return;
				else stream << ',' << (S32)r[i];
			}
		}
	}

	template <class T, class Y>
	void __toStringSignedSmall(buffer_type& buffer, const U32& stride) const{
		if(*(const Y* const)this->currentAt<T>() == this->__end_of_vector_value){
			buffer += '.';
			return;
		}

		if(stride == 1){
			if(*(const Y* const)this->currentAt<T>() == this->__missing_value){
				buffer += '.';
				return;
			}

			buffer.AddReadble(this->current<T>());
		} else {
			const T* const r = this->currentAt<T>();
			const Y* const u = reinterpret_cast<const Y* const>(this->currentAt<T>());


			if(u[0] == this->__missing_value) buffer += '.';
			else if(u[0] == this->__end_of_vector_value) return;
			else buffer.AddReadble(r[0]);

			for(U32 i = 1; i < stride; ++i){
				if(u[i] == this->__missing_value) buffer += ",.";
				else if(u[i] == this->__end_of_vector_value) return;
				else {
					buffer += ',';
					buffer.AddReadble(r[i]);
				}
			}
		}
	}

	/**<
	 *
	 * @param stream
	 * @param stride
	 */
	template <class T>
	inline void __toStringFloat(std::ostream& stream, const U32& stride) const{
		this->__toStringUnsigned<T>(stream, stride);
	}

	template <class T>
	inline void __toStringFloat(buffer_type& buffer, const U32& stride) const{
		this->__toStringUnsigned<T>(buffer, stride);
	}

protected:
	S32 position;    // iterator position
	S32 n_entries;   // size
	BYTE type_size;  // sizeof(TYPE)

private:
	U32 __missing_value;
	U32 __end_of_vector_value;
	BYTE __source_sign;
	Core::TACHYON_CORE_TYPE __source_type;
	toStringFunctionDefinition toStringFunction;
	toStringBufferFunctionDefinition toStringBufferFunction;

protected:
	// Buffer reference
	const buffer_type& buffer; // buffer reference
};

template <class T>
class ContainerIteratorType : public ContainerIteratorDataInterface{
private:
	typedef ContainerIteratorDataInterface parent_type;
	typedef ContainerIteratorType self_type;
	typedef const T* const_pointer;
	typedef const T* const const_pointer_final;
	typedef const T& const_reference;

public:
	ContainerIteratorType(const buffer_type& buffer) :
		ContainerIteratorDataInterface(buffer)
	{
		//this->n_entries = buffer.pointer / sizeof(T);
		//assert(buffer.pointer % sizeof(T) == 0);
	}
	~ContainerIteratorType(){}

	/* ****************************************
	*  Basic iterator functions
	******************************************/
	inline const_reference current(void) const           { return(parent_type::current<T>());    }
	inline const_pointer   currentAt(void) const         { return(parent_type::currentAt<T>());  }
	inline const_reference first(void) const             { return(parent_type::first<T>());      }
	inline const_reference last(void) const              { return(parent_type::last<T>());       }
	inline const_reference operator[](const U32& p) const{ return(parent_type::operator[]<T>(p));}
	inline const_pointer   at(const U32& p) const        { return(parent_type::at<T>(p));        }
	inline const U32       getCurrentStride(void) const  {
		return((U32)*reinterpret_cast<const T* const>(&this->buffer.data[this->position*sizeof(T)]));
	}
};

// Special case for booleans
template <>
class ContainerIteratorType<void> : public ContainerIteratorDataInterface{
private:
	typedef ContainerIteratorType<void> self_type;

public:
	ContainerIteratorType(const buffer_type& buffer) :
		ContainerIteratorDataInterface(buffer)
	{
		//this->n_entries = 0;
	}
	~ContainerIteratorType(){}

	inline void* current(void) const{ return(reinterpret_cast<void*>(&this->buffer.data[this->position*this->type_size])); }
	inline void* first(void) const{ return(reinterpret_cast<void*>(&this->buffer.data[0])); }
	inline void* last(void) const{
		if(this->n_entries - 1 < 0) return(this->first());
		return(reinterpret_cast<void*>(&this->buffer.data[(this->n_entries - 1)*this->type_size]));
	}
	void currentPointer(const void*& p) const{ p = reinterpret_cast<const void*>(&this->buffer.data[this->position]); }

	// Dangerous functions
	inline const U32 getCurrentStride(void) const{ return(0); }
};

}
}

#endif /* CONTAINERITERATORINTERFACE_H_ */
