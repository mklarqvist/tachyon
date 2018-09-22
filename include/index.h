/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TACHYON_INDEX_H_
#define TACHYON_INDEX_H_

#include <algorithm>
#include <memory>
#include <vector>

#include "core/header/variant_header.h"
#include "index_record.h"
#include "variant_container.h"

namespace tachyon {

class Index {
public:
	typedef Index        self_type;
    typedef std::size_t  size_type;
    typedef yon1_idx_rec entry_type;
	typedef YonContig    contig_type;

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
	void Setup(const std::vector<VcfContig>& contigs);
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

	/**<
	 * Index all yon1_vnt_t records in a yon1_vc_t container. Iteratively
	 * invokes IndexRecord() on entries.
	 * @param vc       Src yon1_vc_t container.
	 * @param block_id Current incremental yon block id.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool IndexContainer(const yon1_vc_t& vc, const uint32_t block_id);

	/**<
	 * Index a yon1_vnt_t record.
	 * @param rcd      Src yon1_vnt_t record.
	 * @param block_id Current incremental yon block id.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
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
	                 const yon1_vnt_t& rcd);

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
