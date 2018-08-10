#include "genotype_sorter.h"
#include "third_party/xxhash/xxhash.h"

namespace tachyon {
namespace algorithm {

GenotypeSorter::GenotypeSorter() :
	n_samples(0),
	gt_pattern(nullptr)
{
	memset(this->gt_remap, 0, 256*sizeof(uint8_t));
}

GenotypeSorter::~GenotypeSorter(){
	delete [] this->gt_pattern;
}

void GenotypeSorter::SetSamples(const U64 n_samples){
	this->n_samples = n_samples;
	this->permutation_array.Allocate(n_samples);
	this->gt_pattern = new yon_radix_gt[n_samples];
}

void GenotypeSorter::reset(void){
	this->permutation_array.reset();
}

bool GenotypeSorter::Build(const vcf_container_type& vcf_container, io::VcfHeader& vcf_header){
	if(this->GetNumberSamples() == 0)
		return true;

	// Reset the permutation array to [0, n_samples).
	this->permutation_array.reset();

	// Allocate tetraploid worth of memory in the first instance.
	int32_t  largest_ploidy    = 0;
	uint32_t largest_n_alleles = 0;
	uint32_t n_valid_records   = 0;
	for(U32 i = 0; i < vcf_container.sizeWithoutCarryOver(); ++i){
		if(vcf_container[i]->n_fmt == 0) continue;

		// Perform these actions if FORMAT:GT data is available.
		const int& hts_format_key = vcf_container.at(i)->d.fmt[0].id; // htslib IDX value
		if(vcf_header.GetFormat(hts_format_key)->id != "GT"){
			continue;
		}

		++n_valid_records;
		largest_ploidy    = std::max(vcf_container[i]->d.fmt[0].n, largest_ploidy);
		largest_n_alleles = std::max((uint32_t)vcf_container[i]->n_allele + 2, largest_n_alleles);
	}

	// If there are no valid FORMAT:GT entries in the list then
	// return false.
	if(n_valid_records == 0)
		return false;

	// Shift largest_n_alleles one bit to the left as this is how
	// alleles are stored in Vcf as the first bit encodes for
	// the phasing. The phasing may be unphased (0) or phased (1).
	// Add a count of one to the largest alleles as the count should
	// be fully inclusive [0, n_alleles].
	// Add a count of one to the shifted value to represent the phased
	// case as described above.
	const uint32_t largest_n_alleles_binary = (((largest_n_alleles + 1) << 1) + 1);
	const uint8_t  shift_size = ceil(log2(largest_n_alleles_binary));
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
		if(vcf_container[i]->n_fmt == 0) continue;

		// Perform these actions if FORMAT:GT data is available.
		const int& hts_format_key = vcf_container.at(i)->d.fmt[0].id; // htslib IDX value
		if(vcf_header.GetFormat(hts_format_key)->id != "GT"){
			continue;
		}

		// Setup.
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

		U32 n_sample_c = 0;
		for(U32 s = 0; s < sort_helper.size(); ++s){
			for(U32 k = 0; k < bin_used[sort_helper[s].second].size(); ++k, ++n_sample_c){
				permutation_array[n_sample_c] = bin_used[sort_helper[s].second][k]->id;
			}
		}
		assert(n_sample_c == this->GetNumberSamples());

		bin_used.clear();
		bin_used_map.clear();
		bin_used_packed_integer.clear();
		sort_helper.clear();
	}

	//this->Debug(std::cout, vcf_container, permutation_array);

	return true;
}

void GenotypeSorter::Debug(std::ostream& stream, const vcf_container_type& vcf_container, const yon_gt_ppa& ppa){
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
				stream << "|" << (U32)(gt_target[a] >> 1);
			}
			stream << "\t";
		}
		stream << std::endl;
	}
}

} /* namespace IO */
} /* namespace Tachyon */
