#ifndef INDEX_INDEX_H_
#define INDEX_INDEX_H_
#include <bitset>

#include "index_entry_container.h"
#include "index_meta_container.h"
#include "variant_index.h"

namespace tachyon{
namespace index{

class Index{
private:
	typedef Index                 self_type;
    typedef std::size_t           size_type;
	typedef IndexEntryContainer   container_type;
	typedef IndexMetaContainer    container_meta_type;
	typedef IndexEntry            entry_type;
	typedef IndexIndexEntry       entry_meta_type;
	typedef VariantIndex          variant_index_type;

public:
	Index(){}
	Index(const self_type& other) : index_(other.index_), index_meta_(other.index_meta_), variant_index_(other.variant_index_){}
	~Index(){}

	/**<
	 * Builds the meta-index of index entries when
	 * the input data is sorted. This index corresponds to
	 * index entries all belonging to the same contig.
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool buildSuperIndex(void){
		if(this->index_.size() == 0)
			return false;

		entry_meta_type indexindex;
		indexindex(this->index_[0]); // Start reference
		for(U32 i = 1; i < this->index_.size(); ++i){
			if(indexindex == this->index_[i]) // If the blocks share the same contig identifier
				indexindex += this->index_[i];
			else { // Otherwise push one entry onto the chain and start a new reference
				this->index_meta_ += indexindex;
				indexindex(this->index_[i]);
			}
		}

		if(this->index_meta_.size() == 0){
			this->index_meta_ += indexindex;
		}
		else if(indexindex != this->index_meta_.back()){
			this->index_meta_.back() = indexindex;
		}

		return true;
	}

	// Capacity
	inline const bool empty(void) const{ return(this->index_.empty()); }
	const size_t size(void) const{ return(this->index_.size()); }
	const size_t sizeMeta(void) const{ return(this->index_meta_.size()); }

	inline void operator+=(const entry_type& entry){ this->index_ += entry; }
	inline void operator+=(const entry_meta_type& entry){ this->index_meta_ += entry; }

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
	inline const entry_meta_type& findOverlap(const U32& contig_id) const{
		//if(contig_id > this->getMetaIndex().size()) return;
		return(this->getMetaIndex()[contig_id]);
	}

	/**<
	 * Return interval of YON blocks overlapping target tuple (contigID, position, position)
	 * @param contig_id
	 * @param position
	 * @return
	 */
	std::vector<entry_type> findOverlap(const U32& contig_id, const U64& position) const{
		if(contig_id > this->getMetaIndex().size()) std::vector<entry_type>();

		std::vector<entry_type> overlapping_blocks;
		for(U32 i = 0; i < this->getIndex().size(); ++i){
			// Interval overlap?
			// [a, b] overlaps with [x, y] iff b > x and a < y.
			if(position >= this->getIndex()[i].minPosition && position <= this->getIndex()[i].maxPosition){
				overlapping_blocks.push_back(this->getIndex()[i]);
			}
		}

		// We also need to know possible overlaps in the quad-tree:
		// Seek from root to origin in quad-tree for potential overlapping bins with counts > 0
		// this->variant_index_[contig_id].possibleBins(position, position);

		return(overlapping_blocks);
	}

	/**<
	 * Return interval of YON blocks overlapping target tuple (contigID, start_pos, end_pos)
	 * @param contig_id
	 * @param start_pos
	 * @param end_pos
	 * @return
	 */
	std::vector<entry_type> findOverlap(const U32& contig_id, const U64& start_pos, const U64& end_pos) const{
		if(contig_id > this->getMetaIndex().size()) std::vector<entry_type>();

		std::vector<entry_type> overlapping_blocks;
		for(U32 i = 0; i < this->getIndex().size(); ++i){
			// Interval overlap?
			// [a, b] overlaps with [x, y] iff b > x and a < y.
			if(end_pos >= this->getIndex()[i].minPosition && start_pos <= this->getIndex()[i].maxPosition){
				overlapping_blocks.push_back(this->getIndex()[i]);
			}
		}
		return(overlapping_blocks);
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.index_;
		stream << entry.index_meta_;
		stream << entry.variant_index_;
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream >> entry.index_;
		stream >> entry.index_meta_;
		stream >> entry.variant_index_;
		return(stream);
	}

public:
	container_type      index_;
	container_meta_type index_meta_;
	variant_index_type  variant_index_;
};

}
}

#endif /* INDEX_INDEX_H_ */
