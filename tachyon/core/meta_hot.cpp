#include "meta_hot.h"

namespace tachyon{
namespace core{

// ctor
MetaHot::MetaHot() :
	position(0),
	contigID(0)
{}

MetaHot::MetaHot(const self_type& other) :
	controller(other.controller),
	ref_alt(other.ref_alt),
	position(other.position),
	contigID(other.contigID)
{}

MetaHot::MetaHot(self_type&& other) noexcept :
	controller(other.controller),
	ref_alt(other.ref_alt),
	position(other.position),
	contigID(other.contigID)
{

}

MetaHot& MetaHot::operator=(const self_type& other) noexcept{
	this->controller = other.controller;
	this->ref_alt    = other.ref_alt;
	this->position   = other.position;
	this->contigID   = other.contigID;
	return(*this);
}

MetaHot& MetaHot::operator=(self_type&& other) noexcept{
	if (this == &other)
		return *this;

	this->controller = other.controller;
	this->ref_alt    = other.ref_alt;
	this->position   = other.position;
	this->contigID   = other.contigID;
	return(*this);
}

// dtor
MetaHot::~MetaHot(){}

// MetaHotRefAlt
MetaHotRefAlt::MetaHotRefAlt() : ref(0), alt(0){}
MetaHotRefAlt::~MetaHotRefAlt(){}

bool MetaHotRefAlt::setRef(const char& c){
	switch(c){
	case 'A': this->ref = constants::REF_ALT_A; break;
	case 'T': this->ref = constants::REF_ALT_T; break;
	case 'G': this->ref = constants::REF_ALT_G; break;
	case 'C': this->ref = constants::REF_ALT_C; break;
	default:
		std::cerr << utility::timestamp("ERROR", "ENCODING") << "Illegal SNV reference..." << std::endl;
		return false;
	}
	return true;
}

bool MetaHotRefAlt::setAlt(const char& c){
	switch(c){
	case 'A': this->alt = constants::REF_ALT_A; break;
	case 'T': this->alt = constants::REF_ALT_T; break;
	case 'G': this->alt = constants::REF_ALT_G; break;
	case 'C': this->alt = constants::REF_ALT_C; break;
	case 'N': this->alt = constants::REF_ALT_N; break;
	default:
		std::cerr << utility::timestamp("ERROR", "ENCODING") << "Illegal SNV alternative..." << std::endl;
		return false;
	}
	return true;
}

}
}
