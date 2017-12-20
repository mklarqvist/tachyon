#include "MetaCold.h"

namespace Tachyon{
namespace Core{

MetaCold::MetaCold(void) :
	l_body(0),
	QUAL(0),
	n_allele(0),
	n_ID(0),
	ID(nullptr),
	alleles(nullptr)
{}

MetaCold::~MetaCold(void){ delete [] this->alleles; }

bool MetaCold::write(const bcf_type& entry, stream_container& buffer){
	this->l_body = sizeof(float) + sizeof(U16) + sizeof(U16) + entry.l_ID;
	for(U32 i = 0; i < entry.body->n_allele; ++i){
		this->l_body += sizeof(U16) + entry.alleles[i].length;
	}
	this->l_body += sizeof(U32);

	const U64 start = buffer.buffer_data_uncompressed.pointer;
	buffer.buffer_data_uncompressed += this->l_body;
	//std::cerr << "adding: " << this->l_body << std::endl;

	// Write out data
	// offset is
	buffer.buffer_data_uncompressed += entry.body->QUAL;
	buffer.buffer_data_uncompressed += (U16)entry.body->n_allele;

	// Write out ID
	buffer.buffer_data_uncompressed += (U16)entry.l_ID;
	if((U16)entry.l_ID) buffer.buffer_data_uncompressed.Add(entry.ID, entry.l_ID);

	// Write out alleles
	for(U32 i = 0; i < entry.body->n_allele; ++i){
		// Write out allele
		buffer.buffer_data_uncompressed += (U16)entry.alleles[i].length;
		buffer.buffer_data_uncompressed.Add(entry.alleles[i].data, entry.alleles[i].length);
	}

	// Assert length is correct
	assert(buffer.buffer_data_uncompressed.pointer - start == this->l_body);

	return true;
}

}
}
