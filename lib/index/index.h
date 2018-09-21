#ifndef INDEX_INDEX_H_
#define INDEX_INDEX_H_

#include <algorithm>
#include <memory>
#include <vector>

#include "core/header/variant_header.h"
#include "variant_index_entry.h"
#include "core/variant_record.h"

namespace tachyon {

class Index {
public:
	typedef Index             self_type;
    typedef std::size_t       size_type;
    typedef index::VariantIndexEntry entry_type;
	typedef YonContig         contig_type;

public:
	Index();
	Index(const self_type& other);
	~Index();

	/**<
	 * Reduce operator for the VariantIndex. Used during parallel
	 * import from htslib Vcf files.
	 * @param other Other Index object.
	 * @return      Returns the reduced Index object.
	 */
	Index& operator+=(const Index& other);

	// Capacity
	bool empty(void) const;
	size_t size(void) const;
	size_t sizeMeta(void) const;
	uint64_t GetLinearSize(void) const;

	/**<
	 * Debug function for printing out index data.
	 * @param stream Dst stream reference.
	 * @return       Returns the dst stream reference.
	 */
	std::ostream& Print(std::ostream& stream) const;

	/**<
	 * Return interval of YON blocks overlapping target tuple (contigID, position, position).
	 * @param contig_id Reference id as described in the header.
	 * @return          Returns a vector of index entries intersecting the supplied interval.
	 */
	std::vector<entry_type> FindOverlap(const uint32_t& contig_id) const;

	/**<
	 * Return interval of YON blocks overlapping target tuple (contigID, position, position).
	 * @param contig_id Reference id as described in the header.
	 * @param position  Start and end position in 0-base.
	 * @return          Returns a vector of index entries intersecting the supplied interval.
	 */
	inline std::vector<entry_type> FindOverlap(const uint32_t& contig_id, const uint64_t& position) const{
		return(this->FindOverlap(contig_id, position, position));
	}

	/**<
	 * Return interval of YON blocks overlapping target tuple (contigID, start_pos, end_pos).
	 * @param contig_id Reference id as described in the header.
	 * @param start_pos Start position in 0-base.
	 * @param end_pos   End position in 0-base.
	 * @return          Returns a vector of index entries intersecting the supplied interval.
	 */
	std::vector<entry_type> FindOverlap(const uint32_t& contig_id,
	                                    const uint64_t& start_pos,
	                                    const uint64_t& end_pos) const;

	/**<
	 * Wrapper function for using a list of contigs to construct internal
	 * index objects.
	 * @param contigs Src vector of contigs.
	 */
	void Setup(const std::vector<io::VcfContig>& contigs);
	void Setup(const std::vector<YonContig>& contigs);

	int32_t AddSorted(const uint32_t contig_id,
	                  const uint64_t from_position,
	                  const uint64_t to_position,
	                  const uint32_t yon_block_id);

	self_type& operator+=(const entry_type& entry);
	entry_type& operator[](const uint32_t block_id);
	const entry_type& operator[](const uint32_t block_id) const;

	inline entry_type& GetCurrent(void){ return(this->current_entry_); }
	inline const entry_type& GetCurrent(void) const{ return(this->current_entry_); }

	bool IndexRecord(const yon1_vnt_t& rcd, const uint32_t block_id);

	/**<
	 * Index a record using a htslib bcf1_t record. This function is used
	 * only during importing from htslib streams. Requires knowledge of the
	 * Info:End field to properly support indexing of indels/CNVs.
	 * @param record        Src pointer to a htslib bcf1_t record.
	 * @param block_id      Current incremental yon block id.
	 * @param info_end_key  Array offset of Info:End field in the htslib header. Set to -1 if unavaiable.
	 * @param rcd           Src yon1_vnt_t record providing additional information.
	 * @return              Returns TRUE upon success or FALSE otherwise.
	 */
	bool IndexRecord(const bcf1_t*  record,
	                 const uint32_t block_id,
	                 const int32_t  info_end_key,
	                 const yon1_vnt_t& rcd)
	{
		assert(record != nullptr);
		int32_t index_bin = -1;

		// Ascertain that the meta entry has been evaluated
		// prior to executing this function.
		if(rcd.n_alleles == 0){
			std::cerr << utility::timestamp("ERROR","IMPORT") << "The target meta record must be parsed prior to executing indexing functions..." << std::endl;
			return false;
		}

		int64_t end_position_used = record->pos;

		// The Info field END is used as the end position of an internal if it is available. This field
		// is usually only set for non-standard variants such as SVs or other special meaning records.
		if(info_end_key != -1){
			// Linear search for the END key: this is not optimal but is probably faster
			// than first constructing a hash table for each record.
			const int n_info_fields = record->n_info;

			// Iterate over available Info fields.
			for(uint32_t i = 0; i < n_info_fields; ++i){
				if(record->d.info[i].key == info_end_key){
					uint32_t end = 0;

					switch(record->d.info[i].type){
					case(BCF_BT_INT8):  end = *reinterpret_cast<int8_t*> (record->d.info[i].vptr); break;
					case(BCF_BT_INT16): end = *reinterpret_cast<int16_t*>(record->d.info[i].vptr); break;
					case(BCF_BT_INT32): end = *reinterpret_cast<int32_t*>(record->d.info[i].vptr); break;
					default:
						std::cerr << utility::timestamp("ERROR","INDEX") << "Illegal END primitive type: " << io::BCF_TYPE_LOOKUP[record->d.info[i].type] << std::endl;
						return false;
					}

					index_bin = this->AddSorted(rcd.rid, record->pos, end, block_id);
					break;
				}
			}
		}

		// If the END field cannot be found then we check if the variant is a
		if(index_bin == -1){
			int32_t longest = -1;
			// Iterate over available allele information and find the longest
			// SNV/indel length. The regex pattern ^[ATGC]{1,}$ searches for
			// simple SNV/indels.
			for(uint32_t i = 0; i < rcd.n_alleles; ++i){
				if(std::regex_match(rcd.alleles[i].allele, YON_REGEX_CANONICAL_BASES)){
					if(rcd.alleles[i].l_allele > longest)
						longest = rcd.alleles[i].l_allele;
				}
			}

			// Update the variant index with the target bin(s) found.
			if(longest > 1){
				index_bin = this->AddSorted(rcd.rid,
				                            record->pos,
				                            record->pos + longest,
				                            block_id);
				//index_bin = 0;
				end_position_used = record->pos + longest;
			}
			// In the cases of special-meaning alleles such as copy-number (e.g. <CN>)
			// or SV (e.g. A[B)) they are index according to their left-most value only.
			// This has the implication that they cannot be found by means of interval
			// intersection searches. If special-meaning variants were to be supported
			// in the index then many more blocks would have to be searched for each
			// query as the few special cases will dominate the many general cases. For
			// this reason special-meaning alleles are not completely indexed.
			else {
				index_bin = this->AddSorted(rcd.rid,
				                            record->pos,
				                            record->pos,
				                            block_id);
			}
		}

		if(index_bin > this->current_entry_.max_bin) this->current_entry_.max_bin = index_bin;
		if(index_bin < this->current_entry_.min_bin) this->current_entry_.min_bin = index_bin;
		if(end_position_used > this->current_entry_.max_position)
			this->current_entry_.max_position = end_position_used;

		// Update number of entries in block
		++this->current_entry_.n_variants;

		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry);
	friend std::istream& operator>>(std::istream& stream, self_type& entry);

private:
	bool is_sorted_;
	entry_type current_entry_;

	// Pimpl idiom
	class IndexImpl;
	std::unique_ptr<IndexImpl> mImpl;
};

}

#endif /* INDEX_INDEX_H_ */
