#include <bitset>
#include "radix_sort_gt.h"

namespace tachyon {
namespace algorithm {

RadixSortGT::RadixSortGT() :
	n_samples(0),
	position(0),
	GT_array(nullptr),
	bins(new U32*[9]),
	manager(nullptr)
{
	memset(&p_i, 0, sizeof(U32)*9);
}

RadixSortGT::~RadixSortGT(){
	delete [] this->GT_array;
	for(U32 i = 0; i < 9; ++i)
		delete [] this->bins[i];
}

void RadixSortGT::setSamples(const U64 n_samples){
	this->n_samples = n_samples;

	// Delete previous
	delete [] this->GT_array;

	// Set new
	this->GT_array = new BYTE[this->n_samples];

	// Reset
	for(U32 i = 0; i < 9; ++i){
		this->bins[i] = new U32[n_samples];
		memset(this->bins[i], 0, sizeof(U32)*n_samples);
	}

	memset(this->GT_array, 0, sizeof(BYTE)*n_samples);

	this->manager->setSamples(n_samples);
}

void RadixSortGT::reset(void){
	this->position = 0;
	memset(this->GT_array, 0, sizeof(BYTE)*n_samples);
	memset(&p_i, 0, sizeof(U32)*9);

	this->manager->reset();
}

bool RadixSortGT::build(const bcf_reader_type& reader){
	if(reader.size() == 0)
		return false;

	// Cycle over BCF entries
	for(U32 i = 0; i < reader.size(); ++i){
		// Have to have genotypes available
		if(reader[i].hasGenotypes == false)
			continue;

		if(reader[i].gt_support.hasEOV || reader[i].gt_support.ploidy != 2)
			continue;

		// Has to be biallelic
		// otherwise skip
		if(!reader[i].isBiallelic())
			continue;

		//std::cerr << "update: " << i << std::endl;
		if(!this->update(reader[i]))
			continue;

	}

	// Return TRUE if the number of parsed
	// entries is > 0
	return(true);
}

bool RadixSortGT::update(const bcf_entry_type& entry){
	// Check again because we might use it
	// iteratively at some point in time
	// i.e. not operating through the
	// build() function
	if(!entry.isBiallelic())
		return false;

	// Cycle over genotypes at this position
	// Ignore phasing at this stage
	//
	// Genotype encodings are thus:
	// 0/0 -> 0000b = 0 -> 0
	// 0/1 -> 0001b = 1 -> 3
	// 0/. -> 0010b = 2	-> 4
	// 1/0 -> 0100b = 4 -> 2
	// 1/1 -> 0101b = 5 -> 1
	// 1/. -> 0110b = 6 -> 5
	// ./0 -> 1000b = 8 -> 6
	// ./1 -> 1001b = 9 -> 7
	// ./. -> 1010b = 10 -> 8
	//
	// Update GT_array
	if(entry.formatID[0].primitive_type != 1){
		std::cerr << "unexpected primitive: " << (int)entry.formatID[0].primitive_type << std::endl;
		return false;
	}

	U32 internal_pos = entry.formatID[0].l_offset;
	U32 k = 0;
	//BYTE* debug = new BYTE[2*this->n_samples];
	for(U32 i = 0; i < 2*this->n_samples; i += 2, ++k){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&entry.data[internal_pos++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&entry.data[internal_pos++]);
		//if(*reinterpret_cast<const BYTE* const>(&fmt_type_value1) == 0x80 || *reinterpret_cast<const BYTE* const>(&fmt_type_value2) == 0x80){
		//	std::cerr << "has missing|eov in radix" << std::endl;
		//	return false;
		//}
		const BYTE packed = (bcf::BCF_UNPACK_GENOTYPE(fmt_type_value2) << 2) | bcf::BCF_UNPACK_GENOTYPE(fmt_type_value1);
		this->GT_array[k] = packed;
		//debug[i] = fmt_type_value1;
		//debug[i+1] = fmt_type_value2;
	}

	// Build PPA
	// 3^2 = 9 state radix sort over
	// states: alleles \in {00, 01, 11}
	// b entries in a TWK block B
	// This is equivalent to a radix sort
	// on the alphabet {0,1,...,8}
	U32 target_ID = 0;
	for(U32 j = 0; j < this->n_samples; ++j){
		// Determine correct bin
		switch(this->GT_array[(*this->manager)[j]]){
		case 0:  target_ID = 0; break;
		case 1:  target_ID = 3; break;
		case 2:  target_ID = 4; break;
		case 4:  target_ID = 2; break;
		case 5:  target_ID = 1; break;
		case 6:  target_ID = 5; break;
		case 8:  target_ID = 6; break;
		case 9:  target_ID = 7; break;
		case 10: target_ID = 8; break;
		default: std::cerr << "illegal in radix" << std::endl; exit(1);
		}

		// Update bin i at position i with ppa[j]
		this->bins[target_ID][this->p_i[target_ID]] = (*this->manager)[j];
		++this->p_i[target_ID];
	} // end loop over individuals at position i
	//delete [] debug;

	// Update PPA data
	// Copy data in sorted order
	U32 cum_pos = 0;
	for(U32 i = 0; i < 9; ++i){
		// Copy data in bin i to current position
		//std::cerr << i << '\t' << cum_pos*sizeof(U32) << '\t' << this->p_i[i] << '\t';
		memcpy(&this->manager->PPA[cum_pos*sizeof(U32)], this->bins[i], this->p_i[i]*sizeof(U32));

		// Update cumulative position and reset
		cum_pos += this->p_i[i];
		this->p_i[i] = 0;
	}
	//std::cerr << std::endl;
	// Make sure the cumulative position
	// equals the number of samples in the
	// dataset
	assert(cum_pos == this->n_samples);

	// Keep track of how many entries we've iterated over
	++this->position;

	return true;
}


} /* namespace IO */
} /* namespace Tachyon */
