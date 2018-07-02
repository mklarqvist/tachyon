#ifndef INDEX_INDEX_H_
#define INDEX_INDEX_H_

#include <algorithm>

#include "index_meta_container.h"
#include "variant_index.h"
#include "variant_index_linear.h"

namespace tachyon{
namespace index{

class Index{
private:
	typedef Index                 self_type;
    typedef std::size_t           size_type;
	typedef VariantIndex          container_type;
	typedef IndexMetaContainer    container_meta_type;
	typedef IndexEntry            entry_type;
	typedef IndexIndexEntry       entry_meta_type;
	typedef VariantIndexBin       bin_type;

public:
	Index() : number_blocks(0){}
	Index(const self_type& other) : number_blocks(other.number_blocks), index_meta_(other.index_meta_), index_(other.index_){}
	~Index(){}

	/**<
	 * Builds the meta-index of index entries when
	 * the input data is sorted. This index corresponds to
	 * index entries all belonging to the same contig.
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool buildMetaIndex(void){
		// This criterion should never be satisfied
		if(this->index_.size() == 0)
			return false;

		for(U32 c = 0; c < this->index_.size(); ++c){ // foreach contig
			if(this->index_.linear_at(c).size() == 0){
				this->index_meta_ += entry_meta_type();
				continue;
			}

			entry_meta_type indexindex;
			indexindex(this->index_.linear_at(c)[0]); // Start reference
			for(U32 i = 1; i < this->index_.linear_at(c).size(); ++i){
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

	// Capacity
	inline const bool empty(void) const{ return(this->index_.empty()); }
	const size_t size(void) const{ return(this->index_.size()); }
	const size_t sizeMeta(void) const{ return(this->index_meta_.size()); }

	//inline void operator+=(const entry_type& entry){ this->index_ += entry; }
	//inline void operator+=(const entry_meta_type& entry){ this->index_meta_ += entry; }

	// Accessors
	inline container_type& getIndex(void){ return(this->index_); }
	inline const container_type& getIndex(void) const{ return(this->index_); }
	inline container_meta_type& getMetaIndex(void){ return(this->index_meta_); }
	inline const container_meta_type& getMetaIndex(void) const{ return(this->index_meta_); }

	// Overlap
	// Answer to the questions:
	// 1) Overlapping bins given a contig
	// 2) Overlapping bins given a contig and a single position
	// 3) Overlapping bins given a contig and a start and end position
	/**<
	 * Return interval of YON blocks overlapping target contig ID
	 * @param contig_id
	 * @return
	 */
	inline std::vector<entry_type> findOverlap(const U32& contig_id) const{
		if(contig_id > this->getMetaIndex().size())
			return(std::vector<entry_type>());

		// Todo
		return(std::vector<entry_type>());
	}

	/**<
	 * Return interval of YON blocks overlapping target tuple (contigID, position, position)
	 * @param contig_id
	 * @param position
	 * @return
	 */
	inline std::vector<entry_type> findOverlap(const U32& contig_id, const U64& position) const{
		return(this->findOverlap(contig_id, position, position));
	}

	/**<
	 * Return interval of YON blocks overlapping target tuple (contigID, start_pos, end_pos)
	 * @param contig_id
	 * @param start_pos
	 * @param end_pos
	 * @return
	 */
	std::vector<entry_type> findOverlap(const U32& contig_id, const U64& start_pos, const U64& end_pos) const{
		if(contig_id > this->getMetaIndex().size())
			return(std::vector<entry_type>());

		if(this->getMetaIndex().at(contig_id).n_blocks == 0)
			return(std::vector<entry_type>());

		// We also need to know possible overlaps in the quad-tree:
		// Seek from root to origin in quad-tree for potential overlapping bins with counts > 0

		// Retrieve vector of bins that might contain the data
		// The possibleBins function does not check if they exist
		std::vector<bin_type> possible_bins = this->index_[contig_id].possibleBins(start_pos, end_pos);
		std::vector<U32> yon_blocks;

		// Check if possible bins exists in the linear index
		for(U32 i = 0; i < possible_bins.size(); ++i){
			// Cycle over the YON blocks this bin have data mapping to
			for(U32 j = 0; j < possible_bins[i].size(); ++j){
				// Check [a, b] overlaps with [x, y] iff b > x and a < y.
				// a = this->getIndex().linear_at(contig_id)[possible_bins[i][j]].minPosition;
				// b = this->getIndex().linear_at(contig_id)[possible_bins[i][j]].maxPosition;
				// x = start_pos;
				// y = end_pos;
				if(this->getIndex().linear_at(contig_id)[possible_bins[i][j]].minPosition < end_pos &&
				   this->getIndex().linear_at(contig_id)[possible_bins[i][j]].maxPosition > start_pos)
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
		yon_blocks_deduped.push_back(this->getIndex().linear_at(contig_id)[yon_blocks[0]]);

		for(U32 i = 1; i < yon_blocks.size(); ++i){
			if(yon_blocks[i] != yon_blocks_deduped.back().blockID){
				yon_blocks_deduped.push_back(this->getIndex().linear_at(contig_id)[yon_blocks[i]]);
			}
		}

		// Debug
		//for(U32 i = 0; i < yon_blocks_deduped.size(); ++i){
		//	yon_blocks_deduped[i].print(std::cerr);
		//	std::cerr << std::endl;
		//}

		return(yon_blocks_deduped);
	}

	inline const U64& current_block_number(void) const{ return(this->number_blocks); }
	inline void operator++(void){ ++this->number_blocks; }

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.index_;
		stream << entry.index_meta_;
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream >> entry.index_;
		stream >> entry.index_meta_;
		return(stream);
	}

public:
	U64 number_blocks;
	container_meta_type index_meta_;
	container_type      index_;
};

}
}

#endif /* INDEX_INDEX_H_ */
