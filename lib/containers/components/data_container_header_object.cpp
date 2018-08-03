#include "data_container_header_object.h"

namespace tachyon{
namespace containers{

DataContainerHeaderObject::DataContainerHeaderObject() :
	stride(1),
	offset(0),
	cLength(0),
	uLength(0),
	eLength(0),
	global_key(-1)
{
	memset(&this->crc[0], 0, MD5_DIGEST_LENGTH);
}

DataContainerHeaderObject::DataContainerHeaderObject(const DataContainerHeaderObject& other) :
	controller(other.controller),
	stride(other.stride),
	offset(other.offset),
	cLength(other.cLength),
	uLength(other.uLength),
	eLength(other.eLength),
	global_key(other.global_key)
{
	memcpy(&this->crc[0], &other.crc[0], MD5_DIGEST_LENGTH);
}

DataContainerHeaderObject::DataContainerHeaderObject(DataContainerHeaderObject&& other) noexcept :
	controller(other.controller),
	stride(other.stride),
	offset(other.offset),
	cLength(other.cLength),
	uLength(other.uLength),
	eLength(other.eLength),
	global_key(other.global_key)
{
	memcpy(&this->crc[0], &other.crc[0], MD5_DIGEST_LENGTH);
}

 // copy assignment
DataContainerHeaderObject& DataContainerHeaderObject::operator=(const DataContainerHeaderObject& other){
	this->controller = other.controller;
	this->stride     = other.stride;
	this->offset     = other.offset;
	this->cLength    = other.cLength;
	this->uLength    = other.uLength;
	this->eLength    = other.eLength;
	memcpy(&this->crc[0], &other.crc[0], MD5_DIGEST_LENGTH);
	this->global_key = other.global_key;
	return *this;
}


/** Move assignment operator */
DataContainerHeaderObject& DataContainerHeaderObject::operator=(DataContainerHeaderObject&& other) noexcept{
	this->controller = other.controller;
	this->stride     = other.stride;
	this->offset     = other.offset;
	this->cLength    = other.cLength;
	this->uLength    = other.uLength;
	this->eLength    = other.eLength;
	memcpy(&this->crc[0], &other.crc[0], MD5_DIGEST_LENGTH);
	this->global_key = other.global_key;
	return *this;
}

DataContainerHeaderObject::~DataContainerHeaderObject(){ }

void DataContainerHeaderObject::reset(void){
	this->controller.clear();
	this->stride     = 1;
	this->offset     = 0;
	this->cLength    = 0;
	this->uLength    = 0;
	memset(&this->crc[0], 0, MD5_DIGEST_LENGTH);
	this->global_key = -1;
}

bool DataContainerHeaderObject::operator==(const self_type& other) const{
	if(this->stride     != other.stride)     return false;
	if(this->offset     != other.offset)     return false;
	if(this->cLength    != other.cLength)    return false;
	if(this->uLength    != other.uLength)    return false;
	if(this->eLength    != other.eLength)    return false;
	if(this->global_key != other.global_key) return false;
	if(this->controller != other.controller) return false;
	for(U32 i = 0; i < MD5_DIGEST_LENGTH; ++i)
		if(this->crc[i] != other.crc[i]) return false;

	return true;
}

SBYTE DataContainerHeaderObject::GetPrimitiveWidth(void) const{
	// We do not care about signedness here
	switch(this->controller.type){
	case(YON_TYPE_UNKNOWN):
	case(YON_TYPE_STRUCT): return(-1);
	case(YON_TYPE_BOOLEAN):
	case(YON_TYPE_CHAR):   return(sizeof(char));
	case(YON_TYPE_8B):     return(sizeof(BYTE));
	case(YON_TYPE_16B):    return(sizeof(U16));
	case(YON_TYPE_32B):    return(sizeof(U32));
	case(YON_TYPE_64B):    return(sizeof(U64));
	case(YON_TYPE_FLOAT):  return(sizeof(float));
	case(YON_TYPE_DOUBLE): return(sizeof(double));
	}
	return 0;
}

}
}
