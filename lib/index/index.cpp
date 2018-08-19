#include <iostream>

#include "index.h"

namespace tachyon{
namespace index{

bool Index::buildMetaIndex(void){
	// This criterion should never be satisfied
	if(this->index_.size() == 0)
		return false;

	for(uint32_t c = 0; c < this->index_.size(); ++c){ // foreach contig
		if(this->index_.linear_at(c).size() == 0){
			this->index_meta_ += entry_meta_type();
			continue;
		}

		entry_meta_type indexindex;
		indexindex(this->index_.linear_at(c)[0]); // Start reference
		for(uint32_t i = 1; i < this->index_.linear_at(c).size(); ++i){
			if(indexindex == this->index_.linear_at(c)[i]) // If the blocks share the same contig identifier
				indexindex += this->index_.linear_at(c)[i];
			else { // Otherwise push one entry onto the chain and start a new reference
				this->index_meta_ += indexindex;
				indexindex(this->index_.linear_at(c)[i]);
			}
		}

		this->index_meta_ += indexindex;
	}

	return true;
}

std::vector<IndexEntry> Index::findOverlap(const uint32_t& contig_id) const{
	if(contig_id > this->getMetaIndex().size())
		return(std::vector<entry_type>());

	std::vector<entry_type> yon_blocks;
	for(uint32_t i = 0; i < this->getIndex().linear_at(contig_id).size(); ++i){
		yon_blocks.push_back(this->getIndex().linear_at(contig_id).at(i));
	}

	return(yon_blocks);
}

std::vector<IndexEntry> Index::findOverlap(const uint32_t& contig_id, const uint64_t& start_pos, const uint64_t& end_pos) const{
	if(contig_id > this->getMetaIndex().size())
		return(std::vector<entry_type>());


	if(this->getMetaIndex().at(contig_id).n_blocks == 0){
		return(std::vector<entry_type>());
	}

	const uint32_t block_offset_start = this->getIndex().linear_at(contig_id).at(0).blockID;

	//std::cerr << "Linear index " << this->getIndex().linear_at(contig_id).size() << std::endl;
	//for(uint32_t i = 0; i < this->getIndex().linear_at(contig_id).size(); ++i){
	//	this->getIndex().linear_at(contig_id).at(i).print(std::cerr) << std::endl;
	//}
	//std::cerr << "block offset: " << block_offset_start << std::endl;

	// We also need to know possible overlaps in the quad-tree:
	// Seek from root to origin in quad-tree for potential overlapping bins with counts > 0

	// Retrieve vector of bins that might contain the data
	// The possibleBins function does not check if they exist
	std::vector<bin_type> possible_chunks = this->index_[contig_id].possibleBins(start_pos, end_pos);
	std::vector<uint32_t> yon_blocks;
	//std::cerr << "Possible chunks: " << possible_chunks.size() << std::endl;

	// Check if possible bins exists in the linear index
	for(uint32_t i = 0; i < possible_chunks.size(); ++i){
		// Cycle over the YON blocks this bin have data mapping to
		for(uint32_t j = 0; j < possible_chunks[i].size(); ++j){
			const uint32_t used_bins = possible_chunks[i][j] - block_offset_start;
			//std::cerr << i << "/" << possible_chunks[i].size() << ",raw bin: " << possible_chunks[i][j] << ", used bin: " << used_bins << "\t";
			//possible_chunks[i].print(std::cerr) << std::endl;
			//std::cerr << "Comparing: " << this->getIndex().linear_at(contig_id)[used_bins].minPosition << "<" << end_pos
			//		  << " and " << this->getIndex().linear_at(contig_id)[used_bins].maxPosition << ">" << start_pos << std::endl;

			// Check [a, b] overlaps with [x, y] iff b > x and a < y.
			// a = this->getIndex().linear_at(contig_id)[possible_bins[i][j]].minPosition;
			// b = this->getIndex().linear_at(contig_id)[possible_bins[i][j]].maxPosition;
			// x = start_pos;
			// y = end_pos;
			if(this->getIndex().linear_at(contig_id)[used_bins].minPosition < end_pos &&
			   this->getIndex().linear_at(contig_id)[used_bins].maxPosition > start_pos)
			{
				yon_blocks.push_back(used_bins);
				//std::cerr << "overlap" << std::endl;
			}
			//else {
			//	std::cerr << "no overlap" << std::endl;
			//}
		}
	}

	// Return nothing if all empty
	if(yon_blocks.size() == 0)
		return(std::vector<entry_type>());

	// Sort to dedupe
	std::sort(yon_blocks.begin(), yon_blocks.end());

	// Dedupe
	std::vector<entry_type> yon_blocks_deduped;
	yon_blocks_deduped.push_back(this->getIndex().linear_at(contig_id)[yon_blocks[0]]);

	for(uint32_t i = 1; i < yon_blocks.size(); ++i){
		if(yon_blocks[i] != yon_blocks_deduped.back().blockID - block_offset_start){
			yon_blocks_deduped.push_back(this->getIndex().linear_at(contig_id)[yon_blocks[i]]);
		}
	}

	// Debug
	//std::cerr << "Final\n" << std::endl;
	//for(uint32_t i = 0; i < yon_blocks_deduped.size(); ++i){
	//	yon_blocks_deduped[i].print(std::cerr);
	//	std::cerr << std::endl;
	//}

	return(yon_blocks_deduped);
}

}
}
