#include "MetaCold.h"

namespace tachyon{
namespace core{

MetaCold::MetaCold(void) :
	l_body(0),
	QUAL(0),
	n_allele(0),
	n_ID(0),
	ID(nullptr),
	alleles(nullptr)
{}

MetaCold::MetaCold(char* in) :
	l_body(*reinterpret_cast<U32*>(in)),
	QUAL(*reinterpret_cast<float*>(&in[sizeof(U32)])),
	n_allele(*reinterpret_cast<U16*>(&in[sizeof(U32) + sizeof(float)])),
	n_ID(*reinterpret_cast<U16*>(&in[sizeof(U32) + sizeof(float) + sizeof(U16)])),
	ID(&in[sizeof(U32) + sizeof(float) + 2*sizeof(U16)]),
	alleles(new allele_type[this->n_allele])
{
	U32 cumpos = sizeof(U32) + sizeof(float) + 2*sizeof(U16) + this->n_ID;
	for(U32 i = 0; i < this->n_allele; ++i){
		this->alleles[i](&in[cumpos]);
		cumpos += this->alleles[i].objectSize();
	}
}

MetaCold::MetaCold(const self_type& other) :
	l_body(other.l_body),
	QUAL(other.QUAL),
	n_allele(other.n_allele),
	n_ID(other.n_ID),
	ID(new char[other.n_ID]),
	alleles(new allele_type[other.n_allele])
{
	memcpy(this->ID, other.ID, other.n_ID);
	for(U32 i = 0; i < this->n_allele; ++i){
		this->alleles[i] = other.alleles[i];
	}
}

MetaCold& MetaCold::operator=(const self_type& other){
	this->l_body = other.l_body;
	this->QUAL = other.QUAL;
	this->n_allele = other.n_allele;
	this->n_ID = other.n_ID;
	delete [] this->alleles;
	this->alleles = new allele_type[other.n_allele];

	for(U32 i = 0; i < other.n_allele; ++i){
		this->alleles[i] = other.alleles[i];
	}

	return *this;
}

MetaCold::~MetaCold(void){ delete [] this->alleles; }

void MetaCold::operator()(char* in){
	const U16 prev_n_allele = this->n_allele;

	// Interpret body
	this->l_body   = *reinterpret_cast<U32*>(in);
	this->QUAL     = *reinterpret_cast<float*>(&in[sizeof(U32)]);
	this->n_allele = *reinterpret_cast<U16*>(&in[sizeof(U32) + sizeof(float)]);
	this->n_ID     = *reinterpret_cast<U16*>(&in[sizeof(U32) + sizeof(float) + sizeof(U16)]);

	if(this->n_ID) this->ID = &in[sizeof(U32) + sizeof(float) + 2*sizeof(U16)];
	else this->ID  = nullptr;

	// Only update if we need to
	// Otherwise keep overwriting
	if(prev_n_allele < this->n_allele){
		delete [] this->alleles;
		this->alleles = new allele_type[this->n_allele];
	}

	// Overload allele data
	U32 cumpos = sizeof(U32) + sizeof(float) + 2*sizeof(U16) + this->n_ID;
	for(U32 i = 0; i < this->n_allele; ++i){
		this->alleles[i](&in[cumpos]);
		cumpos += this->alleles[i].objectSize();
	}
}

bool MetaCold::write(const bcf_type& entry, stream_container& buffer){
	this->l_body = sizeof(float) + sizeof(U16) + sizeof(U16) + entry.l_ID;
	for(U32 i = 0; i < entry.body->n_allele; ++i){
		this->l_body += sizeof(U16) + entry.alleles[i].length;
	}
	this->l_body += sizeof(U32);

	const U64 start = buffer.buffer_data_uncompressed.size();
	buffer.buffer_data_uncompressed += this->l_body;

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
	assert(buffer.buffer_data_uncompressed.size() - start == this->l_body);

	return true;
}

std::vector<std::string> MetaCold::getAlleleStrings(void) const{
	std::vector<std::string> ret;
	for(U32 i = 0; i < this->n_allele; ++i)
		ret.push_back(this->alleles[i].toString());

	return(ret);
}

}
}
