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

void GenotypeSorter::SetSamples(const uint64_t n_samples){
	delete [] this->gt_pattern;
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
	for(uint32_t i = 0; i < vcf_container.sizeWithoutCarryOver(); ++i){
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
		for(uint32_t i = 0; i < this->GetNumberSamples(); ++i)
			this->gt_pattern[i].resize(largest_ploidy + 3);
	}

	// Map a genotype such that missing and sentinel node symbol (EOV)
	// is stored in the back of the order.
	for(uint32_t i = 1; i <= largest_n_alleles; ++i) this->gt_remap[i] = i;
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
	for(uint32_t i = 0; i < vcf_container.sizeWithoutCarryOver(); ++i){
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
		uint32_t gt_offset = 0;

		// Iterate over all available samples.
		for(uint32_t s = 0; s < this->GetNumberSamples(); ++s){
			yon_radix_gt& target_pattern = this->gt_pattern[s];
			target_pattern.n_ploidy = base_ploidy;
			target_pattern.id = s;

			// Iterate over the ploidy for this sample and update
			// the allele for that chromosome in the pattern helper
			// structure.
			for(uint32_t a = 0; a < base_ploidy; ++a, ++gt_offset){
				const uint8_t repacked = (this->gt_remap[gt[gt_offset] >> 1] << 1) | (gt[gt_offset] & 1);
				assert((repacked >> 1) <= largest_n_alleles);
				assert(repacked < largest_n_alleles_binary);
				target_pattern.alleles[a] = repacked;
			}
		}
		assert(gt_offset == bcf->d.fmt[0].p_len);

		// Iterate over all encoded genotypes and assign them
		// to different bins according to their bitpacked values.
		for(uint32_t s = 0; s < this->GetNumberSamples(); ++s){
			// Hash the pattern of alleles
			yon_radix_gt& target_pattern = this->gt_pattern[this->permutation_array[s]];
			const uint64_t hash_pattern = XXH64(target_pattern.alleles,
			                               sizeof(uint16_t) * target_pattern.n_ploidy,
			                               651232);

			// Update const_iterators for the hash mapper.
			it  = bin_used_map.find(hash_pattern);
			end = bin_used_map.cend();
			if(it == end){ // Case: does not exist in map.
				bin_used_map[hash_pattern] = bin_used.size();
				bin_used.push_back(std::vector<yon_radix_gt*>());
				bin_used.back().push_back(&target_pattern);
				bin_used_packed_integer.push_back(target_pattern.GetPackedInteger(shift_size));
			} else { // Case: exist in map.
				bin_used[it->second].push_back(&target_pattern);
			}
		}

		// Sort by bin value.
		for(uint32_t k = 0; k < bin_used.size(); ++k)
			sort_helper.push_back(std::pair<uint64_t,uint32_t>(bin_used_packed_integer[k], k));

		std::sort(sort_helper.begin(), sort_helper.end());

		uint32_t n_sample_c = 0;
		for(uint32_t s = 0; s < sort_helper.size(); ++s){
			for(uint32_t k = 0; k < bin_used[sort_helper[s].second].size(); ++k, ++n_sample_c){
				this->permutation_array[n_sample_c] = bin_used[sort_helper[s].second][k]->id;
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

bool GenotypeSorter::Build(const yon1_vnt_t* rec, const uint32_t n_rec){
	if(this->GetNumberSamples() == 0)
		return true;

	// Reset the permutation array to [0, n_samples).
	this->permutation_array.reset();

	// Allocate tetraploid worth of memory in the first instance.
	int32_t  largest_ploidy    = 0;
	uint32_t largest_n_alleles = 0;
	uint32_t n_valid_records   = 0;
	for(uint32_t i = 0; i < n_rec; ++i){
		if(rec[i].controller.gt_available == false) continue;
		assert(rec[i].gt != nullptr);
		if(rec[i].gt->Expand() == false){
			std::cerr << utility::timestamp("ERROR") << "Failed to unpack Format:GT records into individual records..." << std::endl;
			return false;
		}

		++n_valid_records;
		largest_ploidy    = std::max((int32_t)rec[i].gt->m, largest_ploidy);
		largest_n_alleles = std::max((uint32_t)rec[i].n_alleles + 2, largest_n_alleles);
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
		for(uint32_t i = 0; i < this->GetNumberSamples(); ++i)
			this->gt_pattern[i].resize(largest_ploidy + 3);
	}

	// Map a genotype such that missing and sentinel node symbol (EOV)
	// is stored in the back of the order.
	for(uint32_t i = 1; i <= largest_n_alleles; ++i) this->gt_remap[i] = i;
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
	for(uint32_t i = 0; i < n_rec ; ++i){
		if(rec[i].controller.gt_available == false) continue;
		assert(rec[i].gt->d_exp != nullptr);

		// Setup.
		//const bcf1_t* bcf = vcf_container[i];
		//const uint8_t* gt = bcf->d.fmt[0].p;
		const uint32_t base_ploidy = rec[i].gt->m;

		// Iterate over all available samples.
		for(uint32_t s = 0; s < this->GetNumberSamples(); ++s){
			yon_radix_gt& target_pattern = this->gt_pattern[s];
			target_pattern.n_ploidy = base_ploidy;
			target_pattern.id = s;

			// Iterate over the ploidy for this sample and update
			// the allele for that chromosome in the pattern helper
			// structure.
			for(uint32_t a = 0; a < base_ploidy; ++a){
				const uint8_t repacked = (this->gt_remap[YON_GT_RCD_ALLELE_UNPACK(rec[i].gt->d_exp[s].allele[a])] << 1) | YON_GT_RCD_PHASE(rec[i].gt->d_exp[s].allele[a]);
				assert((repacked >> 1) <= largest_n_alleles);
				assert(repacked < largest_n_alleles_binary);
				target_pattern.alleles[a] = repacked;
			}
		}

		// Iterate over all encoded genotypes and assign them
		// to different bins according to their bitpacked values.
		for(uint32_t s = 0; s < this->GetNumberSamples(); ++s){
			// Hash the pattern of alleles
			yon_radix_gt& target_pattern = this->gt_pattern[this->permutation_array[s]];
			const uint64_t hash_pattern = XXH64(target_pattern.alleles,
										        sizeof(uint16_t) * target_pattern.n_ploidy,
										        651232);

			// Update const_iterators for the hash mapper.
			it  = bin_used_map.find(hash_pattern);
			end = bin_used_map.cend();
			if(it == end){ // Case: does not exist in map.
				bin_used_map[hash_pattern] = bin_used.size();
				bin_used.push_back(std::vector<yon_radix_gt*>());
				bin_used.back().push_back(&target_pattern);
				bin_used_packed_integer.push_back(target_pattern.GetPackedInteger(shift_size));
			} else { // Case: exist in map.
				bin_used[it->second].push_back(&target_pattern);
			}
		}

		// Sort by bin value.
		for(uint32_t k = 0; k < bin_used.size(); ++k)
			sort_helper.push_back(std::pair<uint64_t,uint32_t>(bin_used_packed_integer[k], k));

		std::sort(sort_helper.begin(), sort_helper.end());

		uint32_t n_sample_c = 0;
		for(uint32_t s = 0; s < sort_helper.size(); ++s){
			for(uint32_t k = 0; k < bin_used[sort_helper[s].second].size(); ++k, ++n_sample_c){
				this->permutation_array[n_sample_c] = bin_used[sort_helper[s].second][k]->id;
			}
		}
		assert(n_sample_c == this->GetNumberSamples());

		bin_used.clear();
		bin_used_map.clear();
		bin_used_packed_integer.clear();
		sort_helper.clear();
	}

	return true;
}

void GenotypeSorter::Debug(std::ostream& stream, const vcf_container_type& vcf_container, const yon_gt_ppa& ppa){
	for(uint32_t i = 0; i < vcf_container.sizeWithoutCarryOver(); ++i){
		const bcf1_t* bcf = vcf_container[i];
		const uint8_t* gt = bcf->d.fmt[0].p;
		const uint32_t base_ploidy = bcf->d.fmt[0].n;
		assert(bcf->d.fmt[0].p_len == sizeof(int8_t) * base_ploidy * this->GetNumberSamples());

		// Keep track of buffer position.
		uint32_t gt_offset = 0;

		stream << bcf->pos + 1 << "\t";
		// Iterate over all available samples.
		for(uint32_t s = 0; s < this->GetNumberSamples(); ++s){
			const uint8_t* gt_target = &gt[ppa[s] * sizeof(int8_t) * base_ploidy];

			stream << (uint32_t)(gt_target[0] >> 1);
			// Iterate over the ploidy for this sample and update
			// the allele for that chromosome in the pattern helper
			// structure.
			for(uint32_t a = 1; a < base_ploidy; ++a){
				stream << "|" << (uint32_t)(gt_target[a] >> 1);
			}
			stream << "\t";
		}
		stream << std::endl;
	}
}

} /* namespace IO */
} /* namespace Tachyon */
