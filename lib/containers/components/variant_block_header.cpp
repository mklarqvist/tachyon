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
	stream.write(reinterpret_cast<const char*>(&entry.l_offset_footer),   sizeof(uint32_t));
	stream.write(reinterpret_cast<const char*>(&entry.block_hash),        sizeof(uint64_t));
	stream << entry.controller;
	stream.write(reinterpret_cast<const char*>(&entry.contigID),          sizeof(uint32_t));
	stream.write(reinterpret_cast<const char*>(&entry.minPosition),       sizeof(int64_t));
	stream.write(reinterpret_cast<const char*>(&entry.maxPosition),       sizeof(int64_t));
	stream.write(reinterpret_cast<const char*>(&entry.n_variants),        sizeof(uint32_t));

	return(stream);
}

std::ifstream& operator>>(std::ifstream& stream, VariantBlockHeader& entry){
	stream.read(reinterpret_cast<char*>(&entry.l_offset_footer),   sizeof(uint32_t));
	stream.read(reinterpret_cast<char*>(&entry.block_hash),        sizeof(uint64_t));
	stream >> entry.controller;
	stream.read(reinterpret_cast<char*>(&entry.contigID),          sizeof(uint32_t));
	stream.read(reinterpret_cast<char*>(&entry.minPosition),       sizeof(int64_t));
	stream.read(reinterpret_cast<char*>(&entry.maxPosition),       sizeof(int64_t));
	stream.read(reinterpret_cast<char*>(&entry.n_variants),        sizeof(uint32_t));

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
