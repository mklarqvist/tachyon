#ifndef INDEX_INDEX_H_
#define INDEX_INDEX_H_

#include <algorithm>
#include <memory>

#include "variant_index_meta.h"
#include "variant_index_quad_tree.h"
#include "variant_index_linear.h"

namespace tachyon{
namespace index{

class Index {
public:
	typedef Index                 self_type;
    typedef std::size_t           size_type;
	typedef VariantIndexQuadTree  variant_quad_tree_type;
	typedef VariantIndexMeta      variant_meta_type;
	typedef VariantIndexEntry     entry_type;
	typedef VariantIndexMetaEntry entry_meta_type;
	typedef VariantIndexBin       bin_type;
	typedef YonContig             contig_type;
	typedef VariantIndexLinear    variant_linear_type;
	typedef VariantIndexEntry     linear_entry_type;

public:
	Index() : is_sorted(true){}
	Index(const self_type& other) :
		is_sorted(other.is_sorted),
		index_meta_(other.index_meta_),
		index_(other.index_)
	{}
	~Index(){}

	/**<
	 * Reduce operator for the VariantIndex. Used during parallel
	 * import from htslib Vcf files.
	 * @param other Other Index object.
	 * @return      Returns the reduced Index ojbect.
	 */
	Index& operator+=(const Index& other){
		this->index_ += other.index_;
		return(*this);
	}

	// Capacity
	inline bool empty(void) const{ return(this->index_.empty()); }
	inline size_t size(void) const{ return(this->index_.size()); }
	inline size_t sizeMeta(void) const{ return(this->index_meta_.size()); }

	inline uint64_t GetLinearSize(void) const{ return(this->linear_.size()); }

	std::ostream& Print(std::ostream& stream) const{
		for(int i = 0; i < this->index_meta_.size(); ++i){
			stream << "contig " << i << ". blocks: ";
			uint32_t n_blocks = 0;
			for(int j = 0; j < this->index_[i].size(); ++j){
				n_blocks += this->index_[i][j].size();
			}
			stream << n_blocks;
			this->index_meta_[i].Print(stream);
			stream << std::endl;
		}
		return(stream);
	}

	// Accessors
	inline variant_quad_tree_type& GetQuadIndex(void){ return(this->index_); }
	inline const variant_quad_tree_type& GetQuadIndex(void) const{ return(this->index_); }
	inline variant_meta_type& GetMetaIndex(void){ return(this->index_meta_); }
	inline const variant_meta_type& GetMetaIndex(void) const{ return(this->index_meta_); }
	inline variant_linear_type& GetLinearIndex(void){ return(this->linear_); }
	inline const variant_linear_type& GetLinearIndex(void) const{ return(this->linear_); }

	// Overlap functionality used to answer the questions:
	// 1) Overlapping bins given a contig
	// 2) Overlapping bins given a contig and a single position
	// 3) Overlapping bins given a contig and a start and end position
	/**<
	 * Return interval of YON blocks overlapping target contig ID
	 * @param contig_id
	 * @return
	 */
	std::vector<entry_type> FindOverlap(const uint32_t& contig_id) const;

	/**<
	 * Return interval of YON blocks overlapping target tuple (contigID, position, position)
	 * @param contig_id
	 * @param position
	 * @return
	 */
	inline std::vector<entry_type> FindOverlap(const uint32_t& contig_id, const uint64_t& position) const{
		return(this->FindOverlap(contig_id, position, position));
	}

	/**<
	 * Return interval of YON blocks overlapping target tuple (contigID, start_pos, end_pos)
	 * @param contig_id
	 * @param start_pos
	 * @param end_pos
	 * @return
	 */
	std::vector<entry_type> FindOverlap(const uint32_t& contig_id, const uint64_t& start_pos, const uint64_t& end_pos) const;

	/**<
	 * Wrapper function for adding a list of contigs to the index
	 * @param contigs
	 */
	void Setup(const std::vector<io::VcfContig>& contigs){
		this->index_.Add(contigs);
		this->index_meta_.reserve(contigs.size());
		for(int i = 0; i < contigs.size(); ++i)
			this->index_meta_[i].contigID = contigs[i].idx;
	}

	void Setup(const std::vector<YonContig>& contigs){
		this->index_.Add(contigs);
		this->index_meta_.reserve(contigs.size());
		for(int i = 0; i < contigs.size(); ++i)
			this->index_meta_[i].contigID = contigs[i].idx;
	}

	inline int32_t AddSorted(const uint32_t contig_id,
	                         const uint64_t fromPosition,
	                         const uint64_t toPosition,
	                         const uint32_t yon_block_id)
	{
		return(this->index_[contig_id].Add(fromPosition, toPosition, yon_block_id));
	}

	self_type& operator+=(const linear_entry_type& entry){
		this->linear_ += entry;

		// Update the meta index if the input data is sorted.
		if(this->is_sorted){
			if(this->index_meta_[entry.contig_id].n_variants == 0){
				this->index_meta_[entry.contig_id].byte_offset_begin = entry.byte_offset;
				this->index_meta_[entry.contig_id].minPosition = entry.min_position;
				this->index_meta_[entry.contig_id].start_block = entry.block_id;
			}

			this->index_meta_[entry.contig_id] += entry;
		}
		return(*this);
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		utility::SerializePrimitive(entry.is_sorted, stream);
		stream << entry.index_;
		stream << entry.index_meta_;
		stream << entry.linear_;

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		utility::DeserializePrimitive(entry.is_sorted, stream);
		stream >> entry.index_;
		stream >> entry.index_meta_;
		stream >> entry.linear_;

		return(stream);
	}

private:
	// Pimpl idiom
	//class IndexImpl;
	//std::unique_ptr<IndexImpl> mImpl;

	bool is_sorted;
	variant_meta_type      index_meta_;
	variant_quad_tree_type index_;
	variant_linear_type    linear_;
};

}
}

#endif /* INDEX_INDEX_H_ */
