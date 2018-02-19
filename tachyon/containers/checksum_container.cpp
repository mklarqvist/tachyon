#include "checksum_container.h"

namespace tachyon{
namespace containers{

ChecksumContainer::ChecksumContainer(void) :
		n_entries(0),
		n_capacity(1000),
		__entries(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{

}

ChecksumContainer::ChecksumContainer(const size_type capacity) :
		n_entries(0),
		n_capacity(capacity),
		__entries(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{

}

ChecksumContainer::~ChecksumContainer(){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		((this->__entries + i)->~DigitalDigestPair)();

	::operator delete[](static_cast<void*>(this->__entries));
}

}

}
