#include "permutation_manager.h"

namespace tachyon{
namespace algorithm{

PermutationManager::PermutationManager() :
	n_samples(0),
	u_length(0),
	c_length(0),
	crc(0)
{}

PermutationManager::PermutationManager(const U32 n_samples) :
	n_samples(n_samples),
	u_length(0),
	c_length(0),
	crc(0),
	PPA(sizeof(U32)*n_samples)
{}
PermutationManager::~PermutationManager(){ this->PPA.deleteAll(); }

void PermutationManager::setSamples(const U32 n_samples){
	this->n_samples = n_samples;
	this->PPA.reset();
	this->PPA.resize(sizeof(S32)*n_samples);

	for(U32 i = 0; i < this->n_samples; ++i)
		this->PPA += (U32)i;
}

void PermutationManager::reset(void){
	for(U32 i = 0; i < this->n_samples; ++i)
		(*this)[i] = i;

	this->PPA.n_chars = this->n_samples*sizeof(U32);
	this->u_length = 0;
	this->c_length = 0;
	this->crc = 0;
}

bool PermutationManager::generateCRC(void){
	// Checksum for main buffer
	U32 crc = crc32(0, NULL, 0);
	crc = crc32(crc, (Bytef*)this->PPA.buffer, this->n_samples*sizeof(U32));
	this->crc = crc;
	return true;
}

}
}
