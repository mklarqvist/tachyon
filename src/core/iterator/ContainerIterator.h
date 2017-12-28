#ifndef CORE_ITERATOR_CONTAINERITERATOR_H_
#define CORE_ITERATOR_CONTAINERITERATOR_H_

namespace Tachyon{
namespace Core{
namespace Iterator{

/**
 * Interface for iterators over arbitrary containers
 *
 */
class ContainerIteratorDataInterface{
private:
	typedef ContainerIteratorDataInterface self_type;

protected:
	typedef IO::BasicBuffer buffer_type;
	// Function pointers
	typedef void (self_type::*toStringFunctionDefinition)(std::ostream& stream, const U32& stride);


public:
	ContainerIteratorDataInterface(const buffer_type& buffer) :
		position(0),
		n_entries(0),
		type_size(1),
		__missing_value(0x80),
		__end_of_vector_value(0x81),
		__source_sign(0),
		__source_type(Core::CORE_TYPE(0)),
		toStringFunction(nullptr),
		buffer(buffer)
	{

	}

	virtual ~ContainerIteratorDataInterface(){}

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


	/**<
	 *
	 * @param type
	 * @param signedness
	 */
	void setType(const U16& type, const BYTE signedness){
		if(type == Core::TYPE_BOOLEAN){
			this->type_size = 0;
		} else if(type == Core::TYPE_CHAR){
			this->type_size = sizeof(char);
			this->__missing_value = 0x80;
			this->__end_of_vector_value = 0x81;
			this->toStringFunction = &self_type::__toStringNoSeparator<char>;
		} else if(type == Core::TYPE_8B){
			this->type_size = sizeof(BYTE);
			this->__missing_value = 0x80;
			this->__end_of_vector_value = 0x81;
			if(signedness) this->toStringFunction = &self_type::__toStringSignedSmall<SBYTE, BYTE>;
			else this->toStringFunction = &self_type::__toStringUnsignedSmall<BYTE>;
		} else if(type == Core::TYPE_16B){
			this->type_size = sizeof(U16);
			this->__missing_value = 0x8000;
			this->__end_of_vector_value = 0x8001;
			if(signedness) this->toStringFunction = &self_type::__toStringSigned<S16, U16>;
			else this->toStringFunction = &self_type::__toStringUnsigned<U16>;
		} else if(type == Core::TYPE_32B){
			this->type_size = sizeof(U32);
			this->__missing_value = 0x80000000;
			this->__end_of_vector_value = 0x80000001;
			if(signedness) this->toStringFunction = &self_type::__toStringSigned<S32, U32>;
			else this->toStringFunction = &self_type::__toStringUnsigned<U32>;
		} else if(type == Core::TYPE_FLOAT){
			this->type_size = sizeof(float);
			this->__missing_value = 0;
			this->__end_of_vector_value = 0;
			this->toStringFunction = &self_type::__toStringFloat<float>;
		} else if(type == Core::TYPE_DOUBLE){
			this->type_size = sizeof(double);
			this->__missing_value = 0;
			this->__end_of_vector_value = 0;
			this->toStringFunction = &self_type::__toStringFloat<double>;
		} else if(type == Core::TYPE_64B){
			this->type_size = sizeof(U64);
			this->__missing_value = 0;
			this->__end_of_vector_value = 0;
			this->toStringFunction = &self_type::__toStringUnsigned<U64>;
		} else {
				std::cerr << Helpers::timestamp("ERROR") << std::endl;
				exit(1);
		}

		if(type == Core::TYPE_BOOLEAN){
			this->n_entries = 0;
		} else {
			this->n_entries = buffer.pointer / this->type_size;
			assert(buffer.pointer % this->type_size == 0);
		}

		this->__source_type = Core::CORE_TYPE(type);
		this->__source_sign = signedness;
	}

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

	// Virtual functions
	// Has to be overloaded in base class
	inline void toString(std::ostream& stream, const U32& stride){ (this->*toStringFunction)(stream, stride); }

private:
	template <class T>
	void __toStringNoSeparator(std::ostream& stream, const U32& stride){
		if(stride == 1){
			stream.put(this->current<T>());
		} else {
			stream.write(this->currentAt<T>(), stride);
		}
	}

	template <class T>
	void __toStringUnsigned(std::ostream& stream, const U32& stride){
		if(stride == 1){
			stream << this->current<T>();
		} else {
			const T* const r = this->currentAt<T>();
			for(U32 i = 0; i < stride - 1; ++i){
				stream << r[i] << ',';
			}
			stream << r[stride - 1];
		}
	}

	template <class T, class Y>
	void __toStringSigned(std::ostream& stream, const U32& stride){
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
			for(U32 i = 0; i < stride - 1; ++i){
				if(u[i] == this->__missing_value) stream << ".,";
				else if(u[i] == this->__end_of_vector_value) return;
				else stream << r[i] << ',';
			}
			if(u[stride - 1] == this->__missing_value) stream << ".";
			else stream << r[stride - 1];
		}
	}

	template <class T>
	void __toStringUnsignedSmall(std::ostream& stream, const U32& stride){
		if(stride == 1){
			stream << (U32)this->current<T>();
		} else {
			const T* const r = this->currentAt<T>();
			for(U32 i = 0; i < stride - 1; ++i){
				stream << (U32)r[i] << ',';
			}
			stream << (U32)r[stride - 1];
		}
	}

	template <class T, class Y>
	void __toStringSignedSmall(std::ostream& stream, const U32& stride){
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
			for(U32 i = 0; i < stride - 1; ++i){
				if(u[i] == this->__missing_value) stream << ".,";
				else if(u[i] == this->__end_of_vector_value) return;
				else stream << (S32)r[i] << ',';
			}
			if(u[stride - 1] == this->__missing_value) stream << ".";
			else stream << (S32)r[stride - 1];
		}
	}

	template <class T>
	void __toStringFloat(std::ostream& stream, const U32& stride){
		this->__toStringUnsigned<T>(stream, stride);
	}

public:
	S32 position;        // iterator position
	S32 n_entries;       // size
	BYTE type_size;      // sizeof(TYPE)

protected:
	U32 __missing_value;
	U32 __end_of_vector_value;
	BYTE __source_sign;
	Core::CORE_TYPE __source_type;
	toStringFunctionDefinition toStringFunction;


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

	inline const_reference current(void) const      { return(parent_type::current<T>()); }
	inline const_pointer   currentAt(void) const    { return(parent_type::currentAt<T>()); }
	inline const_reference first(void) const        { return(parent_type::first<T>()); }
	inline const_reference last(void) const         { return(parent_type::last<T>()); }
	inline const_reference operator[](const U32& p) const{ return(parent_type::operator[]<T>(p)); }
	inline const_pointer   at(const U32& p) const   { return(parent_type::at<T>(p)); }
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

	template <class T>
	inline const T& operator[](const U32& p) const{ return(*reinterpret_cast<const T* const>(&this->buffer.data[p*sizeof(T)])); }

	template <class T>
	inline const T* at(const U32& p) const{ return( reinterpret_cast<const T*>(&this->buffer.data[p*sizeof(T)])); }

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

class ContainerIterator{
private:
	typedef ContainerIterator              self_type;
	typedef ContainerIteratorDataInterface data_iterator_type;
	typedef ContainerIteratorDataInterface stride_iterator_type;

protected:
	typedef StreamContainer container_type;

public:
	explicit ContainerIterator(void) :
		position(0),
		hasStrideIteratorSet(false),
		container(nullptr),
		data_iterator(nullptr),
		stride_iterator(nullptr)
	{

	}

	ContainerIterator(const container_type& container) :
		position(0),
		hasStrideIteratorSet(false),
		container(&container),
		data_iterator(nullptr),
		stride_iterator(nullptr)
	{
		this->setup(container);
	}

	~ContainerIterator(){
		// Do not delete container. It's allocated externally
		delete this->data_iterator;
		delete this->stride_iterator;
	}

	/**< @brief Overloaded operator for setup() synonym
	 *
	 * @param container
	 */
	inline void operator()(const container_type& container){ this->setup(container); }

	/**< @brief Setup type-specific iterators for both data and stride
	 *
	 * @param container
	 */
	void setup(const container_type& container){
		// Store container reference
		this->container = &container;

		// Recycling is possible
		delete this->data_iterator;   this->data_iterator   = nullptr;
		delete this->stride_iterator; this->stride_iterator = nullptr;
		this->hasStrideIteratorSet = false;
		this->position = 0;

		// Factory
		if(container.header.controller.signedness == false){
			switch(container.header.controller.type){
			case(Core::TYPE_8B):     this->data_iterator = new ContainerIteratorType<BYTE>(this->container->buffer_data_uncompressed);   break;
			case(Core::TYPE_16B):    this->data_iterator = new ContainerIteratorType<U16>(this->container->buffer_data_uncompressed);    break;
			case(Core::TYPE_32B):    this->data_iterator = new ContainerIteratorType<U32>(this->container->buffer_data_uncompressed);    break;
			case(Core::TYPE_64B):    this->data_iterator = new ContainerIteratorType<U64>(this->container->buffer_data_uncompressed);    break;
			case(Core::TYPE_FLOAT):  this->data_iterator = new ContainerIteratorType<float>(this->container->buffer_data_uncompressed);  break;
			//case(Core::TYPE_DOUBLE): this->data_iterator = new ContainerIteratorType<double>(this->container->buffer_data_uncompressed); break;
			case(Core::TYPE_BOOLEAN):this->data_iterator = new ContainerIteratorType<void>(this->container->buffer_data_uncompressed);   break;
			default: std::cerr << Helpers::timestamp("ERROR") << "Illegal type" << std::endl; exit(1); break;
			}
		} else {
			switch(container.header.controller.type){
			case(Core::TYPE_CHAR):this->data_iterator = new ContainerIteratorType<char>(this->container->buffer_data_uncompressed); break;
			case(Core::TYPE_8B):  this->data_iterator = new ContainerIteratorType<char>(this->container->buffer_data_uncompressed); break;
			case(Core::TYPE_16B): this->data_iterator = new ContainerIteratorType<S16>(this->container->buffer_data_uncompressed);  break;
			case(Core::TYPE_32B): this->data_iterator = new ContainerIteratorType<S32>(this->container->buffer_data_uncompressed);  break;
			default: std::cerr << Helpers::timestamp("ERROR") << "Illegal type" << std::endl; exit(1); break;
			}
		}
		this->data_iterator->setType(this->container->header.controller.type, this->container->header.controller.signedness);

		// Construct this iterator if there is a mixed stride
		// otherwise we assume that the stride size is one (1)
		// Factory
		if(this->container->header.controller.mixedStride){
			this->hasStrideIteratorSet = true;
			switch(this->container->header_stride.controller.type){
				case(Core::TYPE_8B):  this->stride_iterator = new ContainerIteratorType<BYTE>(this->container->buffer_strides_uncompressed); break;
				case(Core::TYPE_16B): this->stride_iterator = new ContainerIteratorType<U16>(this->container->buffer_strides_uncompressed);  break;
				case(Core::TYPE_32B): this->stride_iterator = new ContainerIteratorType<U32>(this->container->buffer_strides_uncompressed);  break;
				case(Core::TYPE_64B): this->stride_iterator = new ContainerIteratorType<U64>(this->container->buffer_strides_uncompressed);  break;
				default: std::cerr << Helpers::timestamp("ERROR") << "Illegal stride type" << std::endl; exit(1); break;
			}

			// Stride iterator
			this->stride_iterator->setType(this->container->header_stride.controller.type, this->container->header_stride.controller.signedness);
		}
	}

	inline const bool isUniform(void) const{ return(this->container->header.controller.uniform); }
	inline const bool isUniformStrides(void) const{ return(this->container->header_stride.controller.uniform); }
	inline const bool hasStrides(void) const{ return(this->container->header.controller.mixedStride); }

	/**<
	 *
	 * @param stream
	 * @param field_name
	 * @return
	 */
	inline bool toString(std::ostream& stream, const std::string& field_name){
		if(field_name.size() == 0){
			std::cerr << "impossible length" << std::endl;
			return false;
		}

		// Inject key string
		stream.write(&field_name[0], field_name.size());
		// Inject an equal sign if the encoded type is not
		// a BOOLEAN
		if(this->container->header.controller.type == Core::TYPE_BOOLEAN)
			return true;

		stream.put('=');

		// Call toString function on iterators
		return(this->toString(stream));
	}

	/**< @brief Returns records from the data stream as a parsed string
	 * Records in the iterator return
	 *
	 * @param stream
	 * @return
	 */
	inline bool toString(std::ostream& stream){
		if(this->stride_iterator != nullptr){
			//std::cerr << "\ncurrent stride: " << this->stride_iterator->getCurrentStride() << std::endl;
			//std::cerr << "params: " << (int)this->container->header.controller.type << '\t' << (int)this->container->header.controller.signedness << std::endl;
			//std::cerr << "pos: " << this->data_iterator->position << std::endl;

			this->data_iterator->toString(stream, this->stride_iterator->getCurrentStride());
		}
		else {
			//std::cerr << "\ncurrent uniform stride: " << this->container->header.stride << std::endl;
			//std::cerr << "params: " << (int)this->container->header.controller.type << '\t' << (int)this->container->header.controller.signedness << std::endl;
			//std::cerr << "pos: " << this->data_iterator->position << std::endl;

			this->data_iterator->toString(stream, this->container->header.stride);
		}

		return(true);
	}

	// Increment operator
	void operator++(void){
		++this->position;
		// If mixed stride: update data with that number, increment stride
		if(this->hasStrideIteratorSet){
			*this->data_iterator += this->stride_iterator->getCurrentStride();
			++(*this->stride_iterator);
		}
		// Uniform stride size: update data with the fixed stride size
		else {
			*this->data_iterator += this->container->header.stride;
		}
	}

	void increment(const bool updateStride = true){
		++this->position;
		// If mixed stride: update data with that number, increment stride
		if(this->hasStrideIteratorSet){
			*this->data_iterator += this->stride_iterator->getCurrentStride();
			if(updateStride)
				++(*this->stride_iterator);
		}
		// Uniform stride size: update data with the fixed stride size
		else {
			*this->data_iterator += this->container->header.stride;
		}
	}

	inline void incrementStride(void){
		if(this->hasStrideIteratorSet) ++(*this->stride_iterator);
	}

	inline data_iterator_type& getDataIterator(void){ return(*this->data_iterator); }

public:
	U32 position;                          // iterator position
	bool hasStrideIteratorSet;             // flag triggered if mixed stride
	const container_type* container;       // reference container
	data_iterator_type* data_iterator;     // data iterator type
	stride_iterator_type* stride_iterator; // stride iterator type
};

}
}
}

#endif /* CORE_ITERATOR_CONTAINERITERATOR_H_ */
