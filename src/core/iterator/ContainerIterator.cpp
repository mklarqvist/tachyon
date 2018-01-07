#include "ContainerIterator.h"

namespace Tachyon{
namespace Iterator{

ContainerIterator::ContainerIterator(void) :
	position(0),
	hasStrideIteratorSet(false),
	status(YON_IT_START),
	container(nullptr),
	data_iterator(nullptr),
	stride_iterator(nullptr)
{

}

ContainerIterator::ContainerIterator(const container_type& container) :
	position(0),
	hasStrideIteratorSet(false),
	status(YON_IT_START),
	container(&container),
	data_iterator(nullptr),
	stride_iterator(nullptr)
{
	this->setup(container);
}

ContainerIterator::~ContainerIterator(){
	// Do not delete container. It's allocated externally
	delete this->data_iterator;
	delete this->stride_iterator;
}

const int ContainerIterator::setup(const container_type& container){
	// Store container reference
	this->container = &container;
	this->status = YON_IT_START;

	// Recycling is possible
	delete this->data_iterator;   this->data_iterator   = nullptr;
	delete this->stride_iterator; this->stride_iterator = nullptr;
	this->hasStrideIteratorSet = false;
	this->position = 0;

	// Factory
	if(container.header.controller.signedness == false){
		switch(container.header.controller.type){
		case(Core::YON_TYPE_8B):     this->data_iterator = new ContainerIteratorType<BYTE>(this->container->buffer_data_uncompressed);   break;
		case(Core::YON_TYPE_16B):    this->data_iterator = new ContainerIteratorType<U16>(this->container->buffer_data_uncompressed);    break;
		case(Core::YON_TYPE_32B):    this->data_iterator = new ContainerIteratorType<U32>(this->container->buffer_data_uncompressed);    break;
		case(Core::YON_TYPE_64B):    this->data_iterator = new ContainerIteratorType<U64>(this->container->buffer_data_uncompressed);    break;
		case(Core::YON_TYPE_FLOAT):  this->data_iterator = new ContainerIteratorType<float>(this->container->buffer_data_uncompressed);  break;
		//case(Core::TYPE_DOUBLE): this->data_iterator = new ContainerIteratorType<double>(this->container->buffer_data_uncompressed); break;
		case(Core::YON_TYPE_BOOLEAN):this->data_iterator = new ContainerIteratorType<void>(this->container->buffer_data_uncompressed);   break;
		default:
			this->status = YON_IT_ERROR_UNKNOWN_TYPE;
			return this->status;
		}
	} else {
		switch(container.header.controller.type){
		case(Core::YON_TYPE_CHAR):this->data_iterator = new ContainerIteratorType<char>(this->container->buffer_data_uncompressed); break;
		case(Core::YON_TYPE_8B):  this->data_iterator = new ContainerIteratorType<char>(this->container->buffer_data_uncompressed); break;
		case(Core::YON_TYPE_16B): this->data_iterator = new ContainerIteratorType<S16>(this->container->buffer_data_uncompressed);  break;
		case(Core::YON_TYPE_32B): this->data_iterator = new ContainerIteratorType<S32>(this->container->buffer_data_uncompressed);  break;
		default:
			this->status = YON_IT_ERROR_UNKNOWN_TYPE;
			return this->status;
		}
	}
	this->data_iterator->setType(this->container->header.controller.type, this->container->header.controller.signedness);

	// Construct this iterator if there is a mixed stride
	// otherwise we assume that the stride size is one (1)
	// Factory
	if(this->container->header.controller.mixedStride){
		this->hasStrideIteratorSet = true;
		switch(this->container->header_stride.controller.type){
			case(Core::YON_TYPE_8B):  this->stride_iterator = new ContainerIteratorType<BYTE>(this->container->buffer_strides_uncompressed); break;
			case(Core::YON_TYPE_16B): this->stride_iterator = new ContainerIteratorType<U16>(this->container->buffer_strides_uncompressed);  break;
			case(Core::YON_TYPE_32B): this->stride_iterator = new ContainerIteratorType<U32>(this->container->buffer_strides_uncompressed);  break;
			case(Core::YON_TYPE_64B): this->stride_iterator = new ContainerIteratorType<U64>(this->container->buffer_strides_uncompressed);  break;
			default:
				this->status = YON_IT_ERROR_UNKNOWN_TYPE;
				return this->status;
		}

		// Stride iterator
		this->stride_iterator->setType(this->container->header_stride.controller.type, this->container->header_stride.controller.signedness);
	}

	return(this->status);
}

const int ContainerIterator::toString(std::ostream& stream, const std::string& field_name){
	if(this->container == nullptr){
		this->status = YON_IT_ERROR_UNINITIALIZED;
		return this->status;
	}

	if(field_name.size() == 0){
		std::cerr << "impossible length" << std::endl;
		this->status = YON_IT_ERROR_GENERAL;
		return this->status;
	}

	// Inject key string
	stream.write(&field_name[0], field_name.size());
	// Inject an equal sign if the encoded type is not
	// a BOOLEAN
	if(this->container->header.controller.type == Core::YON_TYPE_BOOLEAN){
		this->status = YON_IT_GOOD;
		return this->status;
	}

	stream.put('=');

	// Call toString function on iterators
	return(this->toString(stream));
}

const int ContainerIterator::toString(buffer_type& buffer, const std::string& field_name){
	if(this->container == nullptr){
		this->status = YON_IT_ERROR_UNINITIALIZED;
		return this->status;
	}

	if(field_name.size() == 0){
		std::cerr << "impossible length" << std::endl;
		this->status = YON_IT_ERROR_GENERAL;
		return this->status;
	}

	// Inject key string
	buffer.Add(&field_name[0], field_name.size());
	// Inject an equal sign if the encoded type is not
	// a BOOLEAN
	if(this->container->header.controller.type == Core::YON_TYPE_BOOLEAN){
		this->status = YON_IT_GOOD;
		return this->status;
	}

	buffer += '=';

	// Call toString function on iterators
	return(this->toString(buffer));
}

/**< @brief Returns records from the data stream as a parsed string
 * Records in the iterator return
 *
 * @param stream
 * @return
 */
const int ContainerIterator::toString(std::ostream& stream){
	if(this->data_iterator == nullptr){
		this->status = YON_IT_ERROR_UNINITIALIZED;
		return this->status;
	}

	if(this->stride_iterator != nullptr){
		if(this->stride_iterator == nullptr){
			this->status = YON_IT_ERROR_UNINITIALIZED;
			return this->status;
		}
		this->data_iterator->toString(stream, this->stride_iterator->getCurrentStride());
	}
	else {
		this->data_iterator->toString(stream, this->container->header.stride);
	}

	return(YON_IT_GOOD);
}

const int ContainerIterator::toString(buffer_type& buffer){
	if(this->data_iterator == nullptr){
		this->status = YON_IT_ERROR_UNINITIALIZED;
		return this->status;
	}

	if(this->stride_iterator != nullptr){
		if(this->stride_iterator == nullptr){
			this->status = YON_IT_ERROR_UNINITIALIZED;
			return this->status;
		}
		this->data_iterator->toString(buffer, this->stride_iterator->getCurrentStride());
	}
	else {
		this->data_iterator->toString(buffer, this->container->header.stride);
	}

	return(YON_IT_GOOD);
}

// Increment operator
void ContainerIterator::operator++(void){
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

void ContainerIterator::increment(const bool updateStride){
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

}
}


