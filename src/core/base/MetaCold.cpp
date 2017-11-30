#include "MetaCold.h"

namespace Tachyon{
namespace Core{

MetaCold::MetaCold(void) :
		QUAL(0),
		n_allele(0),
		n_ID(0),
		ID(nullptr),
		alleles(nullptr)
	{}

MetaCold::~MetaCold(void){ delete [] this->alleles; }

bool MetaCold::write(const bcf_type& entry, stream_container& buffer){
	// Write out data
	// offset is
	buffer.buffer_data += entry.body->QUAL;
	buffer.buffer_data += (U16)entry.body->n_allele;

	// Write out ID
	buffer.buffer_data += (U16)entry.l_ID;
	buffer.buffer_data.Add(entry.ID, entry.l_ID);

	// Write out alleles
	for(U32 i = 0; i < entry.body->n_allele; ++i){
		// Write out allele
		buffer.buffer_data += (U16)entry.alleles[i].length;
		buffer.buffer_data.Add(entry.alleles[i].data, entry.alleles[i].length);
	}

	return true;
}

}
}
