#ifndef CORE_ITERATOR_CONTAINERITERATOR_H_
#define CORE_ITERATOR_CONTAINERITERATOR_H_

namespace Tachyon{
namespace Core{
namespace Iterator{

class ContainerIteratorDataInterface{
private:
	typedef ContainerIteratorDataInterface self_type;

protected:
	typedef IO::BasicBuffer buffer_type;

public:
	ContainerIteratorDataInterface(const buffer_type& buffer) :
		position(0),
		n_entries(0),
		type_size(1),
		buffer(buffer)
	{

	}

	virtual ~ContainerIteratorDataInterface(){}

	void setType(const U16& type){
		switch(type){
			case(Core::TYPE_BOOLEAN): this->type_size = 0; break;
			case(Core::TYPE_CHAR):
			case(Core::TYPE_8B):      this->type_size = 1; break;
			case(Core::TYPE_16B):     this->type_size = 2; break;
			case(Core::TYPE_32B):
			case(Core::TYPE_FLOAT):   this->type_size = 4; break;
			case(Core::TYPE_DOUBLE):
			case(Core::TYPE_64B):     this->type_size = 8; break;
			default:
				std::cerr << Helpers::timestamp("ERROR") << std::endl;
				exit(1);
				break;
		}

		if(type == Core::TYPE_BOOLEAN){
			this->n_entries = 0;
		} else {
			this->n_entries = buffer.pointer / this->type_size;
			assert(buffer.pointer % this->type_size == 0);
		}
	}

	// Comparison operators
	// Must be overloaded
	virtual const bool operator>(const U32& cmp) const{ return false; }
	virtual const bool operator<(const U32& cmp) const{ return false; }
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
	virtual const U32 getCurrentStride(void) const =0;
	virtual void toString(std::ostream& stream, const U32& stride) =0;
	virtual void toStringNoSeparator(std::ostream& stream, const U32& stride) =0;

public:
	S32 position;        // iterator position
	S32 n_entries;       // size
	BYTE type_size;      // sizeof(TYPE)
	const buffer_type& buffer; // buffer reference
};

template <class T>
class ContainerIteratorType : public ContainerIteratorDataInterface{
private:
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

	inline const_reference current(void) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[this->position*sizeof(T)])); }
	inline const_pointer   currentAt(void) const{ return(reinterpret_cast<const_pointer>(&this->buffer.data[this->position*sizeof(T)])); }
	inline const_reference first(void) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[0])); }
	inline const_reference last(void) const{
		if(this->n_entries - 1 < 0) return(this->first());
		return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[(this->n_entries - 1)*sizeof(T)]));
	}
	inline const_reference operator[](const U32& p) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[p*sizeof(T)])); }
	inline const_pointer   at(const U32& p) const{ return( reinterpret_cast<const_pointer>(&this->buffer.data[p*sizeof(T)])); }

	// Dangerous functions
	inline const U32 getCurrentStride(void) const{ return((U32)*reinterpret_cast<const_pointer_final>(&this->buffer.data[this->position*sizeof(T)])); }

	// Output functions
	void toStringNoSeparator(std::ostream& stream, const U32& stride){}
	void toString(std::ostream& stream, const U32& stride){
		if(stride == 1){
			stream << this->current();
		} else {
			const_pointer r = this->currentAt();
			for(U32 i = 0; i < stride - 1; ++i){
				stream << *r << ',';
				++r;
			}
			stream << *r;
		}
	}
};

template <>
class ContainerIteratorType<char> : public ContainerIteratorDataInterface{
private:
	typedef ContainerIteratorType<char> self_type;
	typedef const char* const_pointer;
	typedef const char* const const_pointer_final;
	typedef const char& const_reference;

public:
	ContainerIteratorType(const buffer_type& buffer) :
		ContainerIteratorDataInterface(buffer)
	{
		//this->n_entries = buffer.pointer / sizeof(char);
		//assert(buffer.pointer % sizeof(char) == 0);
	}
	~ContainerIteratorType(){}

	inline const_reference current(void) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[this->position])); }
	inline const_pointer   currentAt(void) const{ return(reinterpret_cast<const_pointer>(&this->buffer.data[this->position])); }
	inline const_reference first(void) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[0])); }
	inline const_reference last(void) const{
		if(this->n_entries - 1 < 0) return(this->first());
		return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[(this->n_entries - 1)]));
	}
	inline const_reference operator[](const U32& p) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[p])); }
	inline const_pointer   at(const U32& p) const{ return( reinterpret_cast<const_pointer>(&this->buffer.data[p])); }

	// Dangerous functions
	inline const U32 getCurrentStride(void) const{ return((U32)*reinterpret_cast<const_pointer_final>(&this->buffer.data[this->position])); }

	// Output functions
	void toString(std::ostream& stream, const U32& stride){
		if(stride == 1){
			stream << (S32)this->current();
		} else {
			const_pointer_final c = this->currentAt();
			for(U32 i = 0; i < stride - 1; ++i){
				stream << (S32)c[i] << ',';
			}
			stream << (S32)c[stride-1];
		}
	}

	void toStringNoSeparator(std::ostream& stream, const U32& stride){
		if(stride == 1){
			stream.put(this->current());
		} else {
			const_pointer_final c = this->currentAt();
			stream.write(this->currentAt(), stride);
		}
	}
};

template <>
class ContainerIteratorType<BYTE> : public ContainerIteratorDataInterface{
private:
	typedef ContainerIteratorType<BYTE> self_type;
	typedef const BYTE* const_pointer;
	typedef const BYTE* const const_pointer_final;
	typedef const BYTE& const_reference;

public:
	ContainerIteratorType(const buffer_type& buffer) :
		ContainerIteratorDataInterface(buffer)
	{
		//this->n_entries = buffer.pointer / sizeof(BYTE);
		//assert(buffer.pointer % sizeof(BYTE) == 0);
	}
	~ContainerIteratorType(){}

	inline const_reference current(void) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[this->position])); }
	inline const_pointer   currentAt(void) const{ return(reinterpret_cast<const_pointer>(&this->buffer.data[this->position])); }
	inline const_reference first(void) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[0])); }
	inline const_reference last(void) const{
		if(this->n_entries - 1 < 0) return(this->first());
		return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[(this->n_entries - 1)]));
	}
	inline const_reference operator[](const U32& p) const{ return(*reinterpret_cast<const_pointer_final>(&this->buffer.data[p])); }
	inline const_pointer   at(const U32& p) const{ return( reinterpret_cast<const_pointer>(&this->buffer.data[p])); }

	// Dangerous functions
	inline const U32 getCurrentStride(void) const{ return((U32)*reinterpret_cast<const_pointer_final>(&this->buffer.data[this->position])); }

	// Output functions
	void toStringNoSeparator(std::ostream& stream, const U32& stride){}
	void toString(std::ostream& stream, const U32& stride){
		if(stride == 1){
			stream << (U32)this->current();
		} else {
			const_pointer r = this->currentAt();
			for(U32 i = 0; i < stride - 1; ++i){
				stream << (U32)*r << ',';
				++r;
			}
			stream << (U32)*r;
		}
	}
};

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

	// Dangerous functions
	inline const U32 getCurrentStride(void) const{ return(0); }
	// Output
	void toStringNoSeparator(std::ostream& stream, const U32& stride){}
	inline void toString(std::ostream& stream, const U32& stride){}
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

		if(container.header.controller.signedness == false){
			switch(container.header.controller.type){
			case(Core::TYPE_8B):     this->data_iterator = new ContainerIteratorType<BYTE>(this->container->buffer_data_uncompressed);   break;
			case(Core::TYPE_16B):    this->data_iterator = new ContainerIteratorType<U16>(this->container->buffer_data_uncompressed);    break;
			case(Core::TYPE_32B):    this->data_iterator = new ContainerIteratorType<U32>(this->container->buffer_data_uncompressed);    break;
			case(Core::TYPE_64B):    this->data_iterator = new ContainerIteratorType<U64>(this->container->buffer_data_uncompressed);    break;
			case(Core::TYPE_FLOAT):  this->data_iterator = new ContainerIteratorType<float>(this->container->buffer_data_uncompressed);  break;
			case(Core::TYPE_DOUBLE): this->data_iterator = new ContainerIteratorType<double>(this->container->buffer_data_uncompressed); break;
			case(Core::TYPE_BOOLEAN):this->data_iterator = new ContainerIteratorType<void>(this->container->buffer_data_uncompressed);   break;
			default: std::cerr << Helpers::timestamp("ERROR") << "Illegal type" << std::endl; exit(1); break;
			}
		} else {
			switch(container.header.controller.type){
			case(Core::TYPE_CHAR):     this->data_iterator = new ContainerIteratorType<char>(this->container->buffer_data_uncompressed);   break;
			case(Core::TYPE_8B):  this->data_iterator = new ContainerIteratorType<char>(this->container->buffer_data_uncompressed); break;
			case(Core::TYPE_16B): this->data_iterator = new ContainerIteratorType<S16>(this->container->buffer_data_uncompressed);  break;
			case(Core::TYPE_32B): this->data_iterator = new ContainerIteratorType<S32>(this->container->buffer_data_uncompressed);  break;
			default: std::cerr << Helpers::timestamp("ERROR") << "Illegal type" << std::endl; exit(1); break;
			}
		}
		this->data_iterator->setType(this->container->header.controller.type);


		// Construct this iterator if there is a mixed stride
		// otherwise we assume that the stride size is one (1)
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
			this->stride_iterator->setType(this->container->header_stride.controller.type);
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
		if(this->container->header.controller.type != Core::TYPE_CHAR)
			return(this->toString(stream));
		else return(this->toStringNoSeparator(stream));
	}

	/**< @brief Returns records from the data stream as a parsed string
	 * Records in the iterator return
	 *
	 * @param stream
	 * @return
	 */
	inline bool toString(std::ostream& stream){
		if(this->stride_iterator != nullptr)
			this->data_iterator->toString(stream, this->stride_iterator->getCurrentStride());
		else
			this->data_iterator->toString(stream, this->container->header.stride);

		return(true);
	}

	inline bool toStringNoSeparator(std::ostream& stream){
		if(this->stride_iterator != nullptr)
			this->data_iterator->toStringNoSeparator(stream, this->stride_iterator->getCurrentStride());
		else
			this->data_iterator->toStringNoSeparator(stream, this->container->header.stride);

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
