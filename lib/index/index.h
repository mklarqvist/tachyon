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
	Index(const self_type& other) :
		number_blocks(other.number_blocks),
		index_meta_(other.index_meta_),
		index_(other.index_)
	{}
	~Index(){}

	/**<
	 * Builds the meta-index of index entries when
	 * the input data is sorted. This index corresponds to
	 * index entries all belonging to the same contig.
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool buildMetaIndex(void);

	// Capacity
	inline bool empty(void) const{ return(this->index_.empty()); }
	inline size_t size(void) const{ return(this->index_.size()); }
	inline size_t sizeMeta(void) const{ return(this->index_meta_.size()); }

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
	std::vector<entry_type> findOverlap(const U32& contig_id) const;

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
	std::vector<entry_type> findOverlap(const U32& contig_id, const U64& start_pos, const U64& end_pos) const;

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
