#include "radix_sort_gt.h"

#include <bitset>

namespace tachyon {
namespace algorithm {

RadixSortGT::RadixSortGT() :
	n_samples(0),
	position(0),
	GT_array(nullptr),
	bins(new U32*[9]),
	manager(nullptr),
	gt_pattern(nullptr)
{
	for(U32 i = 0; i < 9; ++i) bins[i] = nullptr;
	memset(&p_i, 0, sizeof(U32)*9);
	memset(this->gt_remap, 0, 256*sizeof(uint8_t));
}

RadixSortGT::~RadixSortGT(){
	delete [] this->GT_array;
	for(U32 i = 0; i < 9; ++i) delete [] this->bins[i];
	delete [] this->bins;

	delete [] this->gt_pattern;
}

void RadixSortGT::SetSamples(const U64 n_samples){
	this->n_samples = n_samples;

	/*
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
	*/

	//
	this->permutation_array.Allocate(n_samples);
	this->gt_pattern = new yon_radix_gt[n_samples];
}

void RadixSortGT::reset(void){
	this->position = 0;
	memset(this->GT_array, 0, sizeof(BYTE)*n_samples);
	memset(&p_i, 0, sizeof(U32)*9);
	this->manager->reset();
}

bool RadixSortGT::Build(const vcf_container_type& vcf_container){
	if(this->GetNumberSamples() == 0)
		return true;

	// Reset the permutation array to [0, n_samples).
	this->permutation_array.reset();

	// Allocate tetraploid worth of memory in the first instance.
	int32_t  largest_ploidy    = 0;
	uint32_t largest_n_alleles = 0;
	for(U32 i = 0; i < vcf_container.sizeWithoutCarryOver(); ++i){
		largest_ploidy    = std::max(vcf_container[i]->d.fmt[0].n, largest_ploidy);
		largest_n_alleles = std::max((uint32_t)vcf_container[i]->n_allele + 2, largest_n_alleles);
	}
	// Shift largest_n_alleles one bit to the left as this is how
	// alleles are stored in Vcf as the first bit encodes for
	// the phasing. The phasing may be unphased (0) or phased (1).
	// Add a count of one to the largest alleles as the count should
	// be fully inclusive [0, n_alleles].
	// Add a count of one to the shifted value to represent the phased
	// case as described above.
	const uint32_t largest_n_alleles_binary = (((largest_n_alleles + 1) << 1) + 1);
	const uint8_t  shift_size = ceil(log2(largest_n_alleles_binary));
	//std::cerr << "Largest ploidy: " << largest_ploidy << " largest alleles " << largest_n_alleles << "(" << largest_n_alleles_binary << ")" << std::endl;
	//std::cerr << "Max: " << (U32)shift_size << " bits -> " << floor(64 / shift_size) << " ploidy limit" << std::endl;
	assert(shift_size * largest_ploidy <= 64);

	// Ascertain that enough memory has been allocated.
	if(largest_ploidy >= this->gt_pattern[0].n_allocated){
		for(U32 i = 0; i < this->GetNumberSamples(); ++i)
			this->gt_pattern[i].resize(largest_ploidy + 3);
	}

	// Map a genotype such that missing and sentinel node symbol (EOV)
	// is stored in the back of the order.
	for(U32 i = 1; i <= largest_n_alleles; ++i) this->gt_remap[i] = i;
	this->gt_remap[0]  = largest_n_alleles - 1; // Missing value.
	this->gt_remap[64] = largest_n_alleles;     // Sentinel node symbol in unsigned space (129 >> 1 = 64).

	// In order to keep track of bins that are non-empty we use a
	// vector of pointers to the bins and a map from bins to the vector
	// offsets. The map uses the hash of the alleles as the key.
	std::vector< std::vector<yon_radix_gt*> > bin_used;
	std::unordered_map<uint64_t, uint32_t> bin_used_map;
	// The byte-packed integer is used to determine the relative sort
	// order of the bins.
	std::vector<uint64_t> bin_used_packed_integer;
	// At the end of constructing the vectors of potential genotypic
	// bins there is no guarantee at all that they are in order. Thus
	// we sort a tuple (bit-packed integer, and incremental order) and
	// sort on the bit-packed integer and merge bins with help of the
	// original array ordering.
	std::vector< std::pair<uint64_t, uint32_t> > sort_helper;

	// Recycle iterator objects because the constructors appears to be
	// relatively expensive.
	std::unordered_map<uint64_t, uint32_t>::const_iterator it;
	std::unordered_map<uint64_t, uint32_t>::const_iterator end;

	// Iterate over all available bcf1_t records in the container.
	for(U32 i = 0; i < vcf_container.sizeWithoutCarryOver(); ++i){
		const bcf1_t* bcf = vcf_container[i];
		const uint8_t* gt = bcf->d.fmt[0].p;
		const uint32_t base_ploidy = bcf->d.fmt[0].n;
		assert(bcf->d.fmt[0].p_len == sizeof(int8_t) * base_ploidy * this->GetNumberSamples());

		// Keep track of buffer position.
		U32 gt_offset = 0;

		// Iterate over all available samples.
		for(U32 s = 0; s < this->GetNumberSamples(); ++s){
			this->gt_pattern[s].n_ploidy = base_ploidy;
			this->gt_pattern[s].id = s;
			assert(base_ploidy < gt_pattern[s].n_allocated);
			// Iterate over the ploidy for this sample and update
			// the allele for that chromosome in the pattern helper
			// structure.
			for(U32 a = 0; a < base_ploidy; ++a, ++gt_offset){
				const uint8_t repacked = (this->gt_remap[gt[gt_offset] >> 1] << 1) | (gt[gt_offset] & 1);
				assert((repacked >> 1) <= largest_n_alleles);
				assert(repacked < largest_n_alleles_binary);
				this->gt_pattern[s].alleles[a] = repacked;
			}
		}
		assert(gt_offset == bcf->d.fmt[0].p_len);

		// Iterate over all encoded genotypes and assign them
		// to different bins according to their bitpacked values.
		for(U32 s = 0; s < this->GetNumberSamples(); ++s){
			// Hash the pattern of alleles
			const U64 hash_pattern = XXH64(this->gt_pattern[this->permutation_array[s]].alleles,
			                               sizeof(uint16_t) * this->gt_pattern[this->permutation_array[s]].n_ploidy,
			                               651232);
			// Update const_iterators for the hash mapper.
			it  = bin_used_map.find(hash_pattern);
			end = bin_used_map.cend();
			if(it == end){ // Case: does not exist in map.
				bin_used_map[hash_pattern] = bin_used.size();
				bin_used.push_back(std::vector<yon_radix_gt*>());
				bin_used.back().push_back(&this->gt_pattern[this->permutation_array[s]]);
				bin_used_packed_integer.push_back(this->gt_pattern[this->permutation_array[s]].GetPackedInteger(shift_size));
			} else { // Case: exist in map.
				bin_used[it->second].push_back(&this->gt_pattern[this->permutation_array[s]]);
			}
		}

		// Sort by bin value.
		for(U32 k = 0; k < bin_used.size(); ++k)
			sort_helper.push_back(std::pair<uint64_t,uint32_t>(bin_used_packed_integer[k], k));

		std::sort(sort_helper.begin(), sort_helper.end());
		//for(U32 k = 0; k < bin_used.size(); ++k)
		//	std::cerr << sort_helper[k].first << ": " << sort_helper[k].second << std::endl;

		U32 n_sample_c = 0;
		for(U32 s = 0; s < sort_helper.size(); ++s){
			for(U32 k = 0; k < bin_used[sort_helper[s].second].size(); ++k, ++n_sample_c){
				permutation_array[n_sample_c] = bin_used[sort_helper[s].second][k]->id;
			}
		}
		assert(n_sample_c == this->GetNumberSamples());

		//std::cerr << "Patterns: " << bin_used_map.size() << " unique: " << bin_used.size() << std::endl;
		//for(U32 i = 0; i < bin_used.size(); ++i){
		//	std::cerr << "Bin " << i << ": n_entries = " << bin_used[i].size() << " packed " << bin_used_packed_integer[i] << " -> " << *bin_used[i].front() << std::endl;
		//}

		bin_used.clear();
		bin_used_map.clear();
		bin_used_packed_integer.clear();
		sort_helper.clear();
	}

	//for(U32 i = 0; i < this->GetNumberSamples(); ++i)
	//	std::cerr << "," << permutation_array[i];
	//std::cerr << std::endl;

	//this->Debug(std::cout, vcf_container, permutation_array);

	return true;
}

void RadixSortGT::Debug(std::ostream& stream, const vcf_container_type& vcf_container, const yon_gt_ppa& ppa){
	for(U32 i = 0; i < vcf_container.sizeWithoutCarryOver(); ++i){
		const bcf1_t* bcf = vcf_container[i];
		const uint8_t* gt = bcf->d.fmt[0].p;
		const uint32_t base_ploidy = bcf->d.fmt[0].n;
		assert(bcf->d.fmt[0].p_len == sizeof(int8_t) * base_ploidy * this->GetNumberSamples());

		// Keep track of buffer position.
		U32 gt_offset = 0;

		stream << bcf->pos + 1 << "\t";
		// Iterate over all available samples.
		for(U32 s = 0; s < this->GetNumberSamples(); ++s){
			const uint8_t* gt_target = &gt[ppa[s] * sizeof(int8_t) * base_ploidy];

			stream << (U32)(gt_target[0] >> 1);
			// Iterate over the ploidy for this sample and update
			// the allele for that chromosome in the pattern helper
			// structure.
			for(U32 a = 1; a < base_ploidy; ++a){
				stream << "," << (U32)(gt_target[a] >> 1);
			}
			stream << " ";
		}
		stream << std::endl;
	}
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
