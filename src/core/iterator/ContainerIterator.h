#ifndef CORE_ITERATOR_CONTAINERITERATOR_H_
#define CORE_ITERATOR_CONTAINERITERATOR_H_

#include "../../io/BasicBuffer.h"
#include "../StreamContainer.h"
#include "IteratorStatus.h"
#include "ContainerIteratorInterface.h"

namespace Tachyon{
namespace Iterator{

class ContainerIterator{
private:
	typedef ContainerIterator              self_type;
	typedef ContainerIteratorDataInterface data_iterator_type;
	typedef ContainerIteratorDataInterface stride_iterator_type;
	typedef IO::BasicBuffer                buffer_type;

protected:
	typedef Core::StreamContainer container_type;

public:
	explicit ContainerIterator(void);
	ContainerIterator(const container_type& container);
	~ContainerIterator();

	/**< @brief Overloaded operator for setup() synonym
	 *
	 * @param container
	 */
	inline void operator()(const container_type& container){ this->setup(container); }

	/**< @brief Setup type-specific iterators for both data and stride
	 *
	 * @param container
	 */
	const int setup(const container_type& container);

	inline const bool isUniform(void){
		if(this->container == nullptr){
			this->status = YON_IT_ERROR_UNINITIALIZED;
			return false;
		}
		return(this->container->header.controller.uniform);
	}

	inline const bool isUniformStrides(void){
		if(this->container == nullptr){
			this->status = YON_IT_ERROR_UNINITIALIZED;
			return false;
		}
		return(this->container->header_stride.controller.uniform);
	}

	inline const bool hasStrides(void){
		if(this->container == nullptr){
			this->status = YON_IT_ERROR_UNINITIALIZED;
			return false;
		}
		return(this->container->header.controller.mixedStride);
	}

	/**<
	 *
	 * @param stream
	 * @param field_name
	 * @return
	 */
	const int toString(std::ostream& stream, const std::string& field_name);

	const int toString(buffer_type& buffer, const std::string& field_name);

	/**< @brief Returns records from the data stream as a parsed string
	 * Records in the iterator return
	 *
	 * @param stream
	 * @return
	 */
	const int toString(std::ostream& stream);

	const int toString(buffer_type& buffer);

	// Increment operator
	void operator++(void);
	void reset(void){
		this->getDataIterator()->reset();
		if(this->hasStrideIteratorSet)
			this->getStrideIterator()->reset();
	}

	void increment(const bool updateStride = true);
	void incrementSecondaryUsage(void);

	inline void incrementStride(void){
		if(this->hasStrideIteratorSet) ++(*this->stride_iterator);
	}

	inline data_iterator_type* getDataIterator(void) const{ return(this->data_iterator); }
	inline stride_iterator_type* getStrideIterator(void) const{ return(this->stride_iterator); }

public:
	U32 position;                           // iterator position
	bool hasStrideIteratorSet;              // flag triggered if mixed stride
	TACHYON_ITERATOR_STATUS status;         //
	const container_type*   container;      // reference container
	data_iterator_type*     data_iterator;  // data iterator type
	stride_iterator_type*   stride_iterator;// stride iterator type
};

}
}

#endif /* CORE_ITERATOR_CONTAINERITERATOR_H_ */
