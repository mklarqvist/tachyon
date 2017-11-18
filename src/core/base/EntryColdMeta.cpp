#include "EntryColdMeta.h"

namespace Tachyon{
namespace Core{

EntryColdMeta::EntryColdMeta(void) :
		l_body(0),
		QUAL(0),
		n_allele(0),
		n_ID(0),
		ID(nullptr),
		alleles(nullptr)
	{}

EntryColdMeta::~EntryColdMeta(void){
	delete [] this->alleles;
}

bool EntryColdMeta::write(const bcf_type& entry, stream_container& buffer){
	// Determine offset
	// Base length
	// offset + QUAL + n_alleles
	U32 offset = sizeof(U32) + sizeof(float) + sizeof(U16);
	const U32 buffer_start = buffer.buffer_data.pointer;

	// ID length
	if(entry.l_ID < 63) offset += sizeof(BYTE);
	else if(entry.l_ID < 256) offset += 2*sizeof(BYTE); // BYTE + BYTE
	else offset += sizeof(BYTE) + sizeof(U16); // BYTE + U16
	offset += entry.l_ID;

	// Allele length
	for(U32 i = 0; i < entry.body->n_allele; ++i){
		if(entry.alleles[i].length < 63) offset += sizeof(BYTE);
		else if(entry.alleles[i].length < 256) offset += 2*sizeof(BYTE); // BYTE + BYTE
		else offset += sizeof(BYTE) + sizeof(U16); // BYTE + U16
		offset += entry.alleles[i].length;
	}

	// Write out data
	// offset is
	buffer.buffer_data += offset;
	buffer.buffer_data += entry.body->QUAL;
	buffer.buffer_data += (U16)entry.body->n_allele;

	// Write out ID
	typed_value n_ID;
	if(entry.l_ID < 63){
		n_ID.type = typed_value::BYTE_TYPE;
		n_ID.length = entry.l_ID;
		buffer.buffer_data += n_ID;
	} else if(entry.l_ID < 256){
		n_ID.type = typed_value::BYTE_TYPE;
		n_ID.length = 63;
		buffer.buffer_data += n_ID;
		buffer.buffer_data += (BYTE)entry.l_ID;
	} else{
		n_ID.type = typed_value::U16_TYPE;
		n_ID.length = 63;
		buffer.buffer_data += n_ID;
		buffer.buffer_data += (U16)entry.l_ID;
	}
	buffer.buffer_data.Add(entry.ID, entry.l_ID);

	// Write out alleles
	for(U32 i = 0; i < entry.body->n_allele; ++i){
		// Write out allele
		typed_value n_ID;
		if(entry.alleles[i].length < 63){
			n_ID.type = typed_value::BYTE_TYPE;
			n_ID.length = entry.alleles[i].length;
			buffer.buffer_data += n_ID;
		} else if(entry.alleles[i].length < 256){
			n_ID.type = typed_value::BYTE_TYPE;
			n_ID.length = 63;
			buffer.buffer_data += n_ID;
			buffer.buffer_data += (BYTE)entry.alleles[i].length;
		} else{
			n_ID.type = typed_value::U16_TYPE;
			n_ID.length = 63;
			buffer.buffer_data += n_ID;
			buffer.buffer_data += (U16)entry.alleles[i].length;
		}
		buffer.buffer_data.Add(entry.alleles[i].data, entry.alleles[i].length);
	}

	assert((S32)buffer.buffer_data.pointer - (buffer_start + offset) == 0);

	return true;
}

}
}
