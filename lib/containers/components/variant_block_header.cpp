#include "variant_block_header.h"

#include <cmath>
#include <cassert>


namespace tachyon{
namespace containers{

VariantBlockHeader::VariantBlockHeader() :
	l_offset_footer(0),
	block_hash(0),
	contigID(-1),
	minPosition(0),
	maxPosition(0),
	n_variants(0)
{}

VariantBlockHeader::~VariantBlockHeader(){}

std::ostream& operator<<(std::ostream& stream, const VariantBlockHeader& entry){
	stream.write(reinterpret_cast<const char*>(&entry.l_offset_footer),   sizeof(U32));
	stream.write(reinterpret_cast<const char*>(&entry.block_hash),        sizeof(U64));
	stream << entry.controller;
	stream.write(reinterpret_cast<const char*>(&entry.contigID),          sizeof(U32));
	stream.write(reinterpret_cast<const char*>(&entry.minPosition),       sizeof(S64));
	stream.write(reinterpret_cast<const char*>(&entry.maxPosition),       sizeof(S64));
	stream.write(reinterpret_cast<const char*>(&entry.n_variants),        sizeof(U32));

	return(stream);
}

std::ifstream& operator>>(std::ifstream& stream, VariantBlockHeader& entry){
	stream.read(reinterpret_cast<char*>(&entry.l_offset_footer),   sizeof(U32));
	stream.read(reinterpret_cast<char*>(&entry.block_hash),        sizeof(U64));
	stream >> entry.controller;
	stream.read(reinterpret_cast<char*>(&entry.contigID),          sizeof(U32));
	stream.read(reinterpret_cast<char*>(&entry.minPosition),       sizeof(S64));
	stream.read(reinterpret_cast<char*>(&entry.maxPosition),       sizeof(S64));
	stream.read(reinterpret_cast<char*>(&entry.n_variants),        sizeof(U32));

	return(stream);
}

void VariantBlockHeader::reset(void){
	this->l_offset_footer    = 0;
	this->block_hash         = 0;
	this->controller.clear();
	this->contigID           = -1;
	this->minPosition        = 0;
	this->maxPosition        = 0;
	this->n_variants         = 0;
}

}
}
