#include "containers/genotype_container.h"
#include "containers/primitive_container.h"
#include "containers/stride_container.h"
#include "support/enums.h"

namespace tachyon{
namespace containers{

GenotypeContainer::GenotypeContainer(void) :
	n_entries(0),
	__iterators(nullptr)
{

}

GenotypeContainer::~GenotypeContainer(){
	for(std::size_t i = 0; i < this->n_entries; ++i)
		(this->__iterators + i)->~GenotypeContainerInterface();

	::operator delete[](static_cast<void*>(this->__iterators));
}

}
}
