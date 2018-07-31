#include "permutation_manager.h"

namespace tachyon{
namespace algorithm{

PermutationManager::PermutationManager() :
	n_samples(0)
{}

PermutationManager::PermutationManager(const U32 n_samples) :
	n_samples(n_samples),
	PPA(sizeof(U32)*n_samples)
{}

PermutationManager::~PermutationManager(){ }

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

	this->header.reset();
	this->PPA.n_chars = this->n_samples*sizeof(U32);
}

bool PermutationManager::generateCRC(void){
	exit(1);
}

}
}
