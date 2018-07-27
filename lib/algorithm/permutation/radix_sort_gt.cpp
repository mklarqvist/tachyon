#include "radix_sort_gt.h"

#include <bitset>

namespace tachyon {
namespace algorithm {

RadixSortGT::RadixSortGT() :
	n_samples(0),
	position(0),
	GT_array(nullptr),
	bins(new U32*[9]),
	manager(nullptr)
{
	for(U32 i = 0; i < 9; ++i) bins[i] = nullptr;
	memset(&p_i, 0, sizeof(U32)*9);
}

RadixSortGT::~RadixSortGT(){
	delete [] this->GT_array;
	for(U32 i = 0; i < 9; ++i) delete [] this->bins[i];
	delete [] this->bins;
}

void RadixSortGT::SetSamples(const U64 n_samples){
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

	//this->manager->setSamples(n_samples);
	this->manager->setSamples(n_samples*2);
}

void RadixSortGT::reset(void){
	this->position = 0;
	memset(this->GT_array, 0, sizeof(BYTE)*n_samples);
	memset(&p_i, 0, sizeof(U32)*9);
	this->manager->reset();
}

bool RadixSortGT::Build(const vcf_container_type& vcf_container){
	// Allocate tetraploid worth of memory in the first instance.
	yon_radix_gt* gt_pattern = new yon_radix_gt[this->GetNumberSamples()];

	// Assume largest ALT allele count is 10 and ploidy is 4.
	//uint64_t* temp_bins = new uint64_t[pow(10, 4)];

	// In order to keep track of bins that are non-empty we use a
	// vector of pointers to the bins and a map from bins to the vector
	// offsets. The map uses the hash of the alleles as the key.
	std::vector< std::vector<yon_radix_gt*> > bin_used;
	std::unordered_map<uint64_t, uint32_t> bin_used_map;
	std::vector<uint64_t> bin_used_packed_integer;

	// Alphabet: 0, 1, missing, sentinel symbol (EOV). We know each allele
	// must map to the range [0, n_alleles + 2):
	for(U32 i = 0; i < vcf_container.sizeWithoutCarryOver(); ++i){
		const io::VcfGenotypeSummary g = vcf_container.GetGenotypeSummary(i, this->GetNumberSamples());
		// If the Vcf line has invariant genotypes then simply
		// continue as this line will have no effect.
		if(g.invariant)        continue;
		if(g.base_ploidy == 0) continue;

		// This functionality is only valid if ploidy is under or equal to
		// eight as the packed integer requires 1 byte / ploidy.
		if(g.base_ploidy > 8){
			std::cerr << utility::timestamp("ERROR","PBWT") << "Cannot permute genotypes with base ploidy > 8. Run without permutation activiated!" << std::endl;
			return false;
		}

		const bcf1_t* bcf = vcf_container[i];
		const uint8_t* gt = bcf->d.fmt[0].p;
		//std::cerr << bcf->pos << ": " << bcf->d.fmt[0].p_len << " base ploidy: " << (int)g.base_ploidy << std::endl;
		assert(bcf->d.fmt[0].p_len == sizeof(uint8_t) * g.base_ploidy * this->GetNumberSamples());

		// Base ploidy is equal to fmt.n in htslib.
		U32 gt_offset = 0;
		for(U32 s = 0; s < this->GetNumberSamples(); ++s){
			gt_pattern[s].n_ploidy = g.base_ploidy;
			gt_pattern[s].id = s;
			assert(g.base_ploidy < gt_pattern[s].n_allocated);
			for(U32 a = 0; a < g.base_ploidy; ++a, gt_offset++){
				gt_pattern[s].alleles[a] = (gt[gt_offset] >> 1);
			}
			const U64 hash_pattern = RadixSortGT::HashGenotypes(gt_pattern[s].alleles, gt_pattern[s].n_ploidy);
			//std::cerr << "Hash: " << hash_pattern << std::endl;
			std::unordered_map<uint64_t, uint32_t>::const_iterator it = bin_used_map.find(hash_pattern);
			if(it == bin_used_map.end()){
				bin_used_map[hash_pattern] = bin_used.size();
				bin_used.push_back(std::vector<yon_radix_gt*>());
				bin_used.back().push_back(&gt_pattern[s]);
				bin_used_packed_integer.push_back(gt_pattern[s].GetPackedInteger());
			} else {
				bin_used[it->second].push_back(&gt_pattern[s]);
			}
		}
		assert(gt_offset == bcf->d.fmt[0].p_len);

		std::cerr << "Patterns: " << bin_used_map.size() << " unique: " << bin_used.size() << std::endl;
		for(U32 i = 0; i < bin_used.size(); ++i){
			std::cerr << "Bin " << i << ": n_entries = " << bin_used[i].size() << " packed " << bin_used_packed_integer[i] << " -> " << *bin_used[i].front() << std::endl;
		}

		bin_used.clear();
		bin_used_map.clear();

		// Sort genotypes according to their temporary encoded values.
		//std::sort(&gt_pattern[0], &gt_pattern[this->GetNumberSamples()]);

		/*
		yon_radix_gt* ref = &gt_pattern[0];
		for(U32 k = 1; k < this->GetNumberSamples(); ++k){

		}
		*/

		// Debug
		/*
		yon_radix_gt* ref = &gt_pattern[0];
		U32 l_run = 1;
		for(U32 k = 1; k < this->GetNumberSamples(); ++k){
			if(gt_pattern[k] != *ref){
				std::cerr << i << "/" << vcf_container.sizeWithoutCarryOver() << ": " << *ref << "(" << l_run << ")" << std::endl;
				l_run = 1;
				ref = &gt_pattern[k];
				continue;
			}
			++l_run;
		}
		std::cerr << i << "/" << vcf_container.sizeWithoutCarryOver() << ": " << *ref << "(" << l_run << ")" << std::endl << std::endl;
		*/
	}

	delete [] gt_pattern;

	return true;
}

bool RadixSortGT::Build(const bcf_reader_type& reader){
	if(reader.size() == 0)
		return false;

	// Cycle over BCF entries
	for(U32 i = 0; i < reader.size(); ++i){
		if(!this->Update(reader[i]))
			continue;
	}
	return(true);
}

bool RadixSortGT::Update(const bcf_entry_type& entry){
	// Check again because we might use it
	// iteratively at some point in time
	// i.e. not operating through the
	// build() function
	// Have to have genotypes available
	if(entry.hasGenotypes == false)
		return false;

	if(entry.gt_support.hasEOV || entry.gt_support.ploidy != 2)
		return false;

	// Has to be biallelic
	// otherwise skip
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
		std::cerr << utility::timestamp("ERROR","PERMUTE") << "Illegal primitive: " << (int)entry.formatID[0].primitive_type << std::endl;
		exit(1);
	}

	U32 internal_pos = entry.formatID[0].l_offset;
	U32 k = 0;
	for(U32 i = 0; i < 2*this->n_samples; i += 2, ++k){
		const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&entry.data[internal_pos++]);
		const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&entry.data[internal_pos++]);
		const BYTE packed = (bcf::BCF_UNPACK_GENOTYPE(fmt_type_value2) << 2) | bcf::BCF_UNPACK_GENOTYPE(fmt_type_value1);
		this->GT_array[k] = packed;
	}

	// Build PPA
	// 3^2 = 9 state radix sort over
	// states: alleles \in {00, 01, 10}
	// b entries in a YON block B
	// This is equivalent to a radix sort
	// on the alphabet {0,1,...,8}
	U32 target_ID = 0;
	for(U32 j = 0; j < this->n_samples; ++j){
		// Determine correct bin
		switch(this->GT_array[(*this->manager)[j]]){
		case 0:  target_ID = 0; break; // 0000: Ref, ref
		case 1:  target_ID = 3; break; // 0001: Ref, alt
		case 2:  target_ID = 4; break; // 0010: Ref, Missing
		case 4:  target_ID = 2; break; // 0100: Alt, ref
		case 5:  target_ID = 1; break; // 0101: Alt, alt
		case 6:  target_ID = 5; break; // 0110: Alt, missing
		case 8:  target_ID = 6; break; // 1000: Missing, ref
		case 9:  target_ID = 7; break; // 1001: Missing, alt
		case 10: target_ID = 8; break; // 1010: Missing, missing
		default:
			std::cerr << utility::timestamp("ERROR","PERMUTE") << "Illegal state in radix sort..." << std::endl;
			exit(1);
		}

		// Update bin i at position i with ppa[j]
		this->bins[target_ID][this->p_i[target_ID]] = (*this->manager)[j];
		++this->p_i[target_ID];
	} // end loop over individuals at position i

	// Update PPA data
	// Copy data in sorted order
	U32 cum_pos = 0;
	for(U32 i = 0; i < 9; ++i){
		// Copy data in bin i to current position
		memcpy(&this->manager->PPA[cum_pos*sizeof(U32)], this->bins[i], this->p_i[i]*sizeof(U32));

		// Update cumulative position and reset
		cum_pos += this->p_i[i];
		this->p_i[i] = 0;
	}
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
