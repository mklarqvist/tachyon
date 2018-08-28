#include <iostream>

#include "index.h"

namespace tachyon{
namespace index{

std::vector<VariantIndexEntry> Index::FindOverlap(const uint32_t& contig_id) const{
	if(contig_id > this->GetMetaIndex().size())
		return(std::vector<entry_type>());

	std::vector<entry_type> yon_blocks;
	for(int i = this->index_meta_[contig_id].start_block; i <= this->index_meta_[contig_id].end_block; ++i){
		yon_blocks.push_back(this->linear_[i]);
	}

	return(yon_blocks);
}

std::vector<VariantIndexEntry> Index::FindOverlap(const uint32_t& contig_id,
                                                  const uint64_t& start_pos,
                                                  const uint64_t& end_pos) const
{
	if(contig_id > this->GetMetaIndex().size())
		return(std::vector<entry_type>());

	if(this->GetMetaIndex().at(contig_id).n_blocks == 0)
		return(std::vector<entry_type>());


	// We also need to know possible overlaps in the quad-tree:
	// Seek from root to origin in quad-tree for potential overlapping bins with counts > 0.
	// Retrieve vector of bins that might contain the data of interest.
	// The PossibleBins function does not check if these bins are populated.
	std::vector<bin_type> possible_bins = this->index_[contig_id].PossibleBins(start_pos, end_pos);
	std::vector<uint32_t> yon_blocks;

	// Check if possible bins exists in the linear index
	for(int i = 0; i < possible_bins.size(); ++i){
		// Cycle over the YON blocks this bin have data mapping to
		for(uint32_t j = 0; j < possible_bins[i].size(); ++j){
			if(this->linear_[possible_bins[i][j]].minPosition < end_pos &&
			   this->linear_[possible_bins[i][j]].maxPosition > start_pos)
			{
				yon_blocks.push_back(possible_bins[i][j]);
			}
		}
	}

	// Return nothing if all empty
	if(yon_blocks.size() == 0)
		return(std::vector<entry_type>());

	// Sort to dedupe
	std::sort(yon_blocks.begin(), yon_blocks.end());

	// Dedupe
	std::vector<entry_type> yon_blocks_deduped;
	yon_blocks_deduped.push_back(this->linear_[yon_blocks[0]]);

	for(uint32_t i = 1; i < yon_blocks.size(); ++i){
		if(yon_blocks[i] != yon_blocks_deduped.back().blockID){
			yon_blocks_deduped.push_back(this->linear_[yon_blocks[i]]);
		}
	}

	// Debug
	/*
	std::cerr << "Final\n" << std::endl;
	for(uint32_t i = 0; i < yon_blocks_deduped.size(); ++i){
		yon_blocks_deduped[i].print(std::cerr);
		std::cerr << std::endl;
	}
	*/

	return(yon_blocks_deduped);
}

}
}
