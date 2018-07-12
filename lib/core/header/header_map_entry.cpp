#include "header_map_entry.h"

namespace tachyon{
namespace core{

HeaderMapEntry::HeaderMapEntry() :
	IDX(0),
	primitive_type(0)
{}

HeaderMapEntry::HeaderMapEntry(const std::string& id, const S32& idx) :
	IDX(idx),
	primitive_type(0),
	ID(id)
{}

HeaderMapEntry::HeaderMapEntry(const std::string& id, const S32& idx, const S32& primitive_type) :
	IDX(idx),
	primitive_type(primitive_type),
	ID(id)
{}

HeaderMapEntry::HeaderMapEntry(const std::string& id) :
	IDX(0),
	primitive_type(0),
	ID(id)
{}

HeaderMapEntry::HeaderMapEntry(const self_type& other) :
	IDX(other.IDX),
	primitive_type(other.primitive_type),
	ID(other.ID)
{}

HeaderMapEntry::HeaderMapEntry(self_type&& other) :
	IDX(other.IDX),
	primitive_type(other.primitive_type),
	ID(other.ID)
{}

HeaderMapEntry& HeaderMapEntry::operator=(const self_type& other){
	this->IDX = other.IDX;
	this->primitive_type = other.primitive_type;
	this->ID  = other.ID;
	return(*this);
}

HeaderMapEntry& HeaderMapEntry::operator=(HeaderMapEntry&& other){
	this->IDX = other.IDX;
	this->primitive_type = other.primitive_type;
	this->ID  = other.ID;
	return(*this);
}

}
}
