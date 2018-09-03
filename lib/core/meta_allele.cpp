#include "meta_allele.h"

namespace tachyon{
namespace core{

MetaAllele::MetaAllele() :
	l_allele(0),
	allele(nullptr)
{}

// Ctor from packed byte
MetaAllele::MetaAllele(const char reference) :
	l_allele(1),
	allele(new char[this->l_allele])
{
	this->allele[0] = reference;
}

MetaAllele::MetaAllele(const std::string& reference) :
	l_allele(reference.size()),
	allele(new char[this->l_allele])
{
	memcpy(this->allele, reference.data(), reference.size());
}

// Ctor from buffer
MetaAllele::MetaAllele(const char* const in) :
	l_allele(*reinterpret_cast<const uint16_t* const>(in)),
	allele(new char[this->l_allele])
{
	memcpy(this->allele, &in[sizeof(uint16_t)], this->l_allele);
}

// Ctor directly from buffer object
MetaAllele::MetaAllele(const buffer_type& buffer, const uint32_t position) :
	l_allele(*reinterpret_cast<const uint16_t* const>(&buffer[position])),
	allele(new char[this->l_allele])
{
	memcpy(this->allele, &buffer[position + sizeof(uint16_t)], this->l_allele);
}

MetaAllele::~MetaAllele(void){
	delete [] this->allele;
}

MetaAllele::MetaAllele(const self_type& other) :
	l_allele(other.l_allele),
	allele(new char[other.l_allele])
{
	memcpy(this->allele, other.allele, other.l_allele);
}

MetaAllele::MetaAllele(self_type&& other) :
	l_allele(other.l_allele),
	allele(other.allele)
{
	other.allele = nullptr;
}

MetaAllele& MetaAllele::operator=(const self_type& other){
	this->l_allele = other.l_allele;
	delete [] this->allele;
	this->allele = new char[other.l_allele];
	memcpy(this->allele, other.allele, other.l_allele);
	return(*this);
}

void MetaAllele::operator()(const char* const in){
	this->l_allele = *reinterpret_cast<const uint16_t* const>(in);
	delete [] this->allele;
	this->allele = new char[this->l_allele];
	memcpy(this->allele, &in[sizeof(uint16_t)], this->l_allele);
}

void MetaAllele::operator()(const char* const in, const uint32_t length){
	this->l_allele = length;
	delete [] this->allele;
	this->allele = new char[this->l_allele];
	memcpy(this->allele, &in[0], this->l_allele);
}

void MetaAllele::operator()(const std::string& in){
	this->l_allele = in.size();
	delete [] this->allele;
	this->allele = new char[this->l_allele];
	memcpy(this->allele, in.data(), this->l_allele);
}

}
}
