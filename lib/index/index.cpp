#include <iostream>

#include "index.h"
#include "variant_index_meta.h"
#include "variant_index_quad_tree.h"
#include "variant_index_linear.h"

namespace tachyon {

class yon_index_t::IndexImpl {
public:
	typedef index::VariantIndexQuadTree  variant_quad_tree_type;
	typedef index::VariantIndexMeta      variant_meta_type;
	typedef index::VariantIndexMetaEntry entry_meta_type;
	typedef index::VariantIndexBin       bin_type;
	typedef index::VariantIndexLinear    variant_linear_type;

public:
	IndexImpl() {}
	IndexImpl(const IndexImpl& other) {
		//index_meta_ += other.index_meta_;
		index_  += other.index_;
		//linear_ += other.linear_;
	}
	~IndexImpl() {}

	variant_meta_type      index_meta_;
	variant_quad_tree_type index_;
	variant_linear_type    linear_;
};

yon_index_t::yon_index_t() : mImpl(new IndexImpl), is_sorted_(true) {}

yon_index_t::yon_index_t(const self_type& other) :
	is_sorted_(other.is_sorted_),
	mImpl(new IndexImpl(*other.mImpl))
{
}

yon_index_t& yon_index_t::operator+=(const yon_index_t& other) {
	this->mImpl->index_ += other.mImpl->index_;
	return(*this);
}

yon_index_t::~yon_index_t() {}

bool yon_index_t::empty(void) const { return(this->mImpl->index_.empty()); }
size_t yon_index_t::size(void) const { return(this->mImpl->index_.size()); }
uint64_t yon_index_t::GetLinearSize(void) const { return(this->mImpl->linear_.size()); }

std::vector<yon1_idx_rec> yon_index_t::FindOverlap(const uint32_t& contig_id) const {
	if (contig_id > this->mImpl->index_meta_.size())
		return(std::vector<entry_type>());

	std::vector<entry_type> yon_blocks;
	for (int i = this->mImpl->index_meta_[contig_id].start_block; i <= this->mImpl->index_meta_[contig_id].end_block; ++i) {
		yon_blocks.push_back(this->mImpl->linear_[i]);
	}

	return(yon_blocks);
}

std::vector<yon1_idx_rec> yon_index_t::FindOverlap(const uint32_t& contig_id,
                                                   const uint64_t& start_pos,
                                                   const uint64_t& end_pos) const
{
	if (contig_id > this->mImpl->index_meta_.size())
		return(std::vector<entry_type>());

	if (this->mImpl->index_meta_.at(contig_id).n_blocks == 0)
		return(std::vector<entry_type>());


	// We also need to know possible overlaps in the quad-tree:
	// Seek from root to origin in quad-tree for potential overlapping bins with counts > 0.
	// Retrieve vector of bins that might contain the data of interest.
	// The PossibleBins function does not check if these bins are populated.
	std::vector<IndexImpl::bin_type> possible_bins = this->mImpl->index_[contig_id].PossibleBins(start_pos, end_pos);
	std::vector<uint32_t> yon_blocks;

	// Check if possible bins exists in the linear index
	for (int i = 0; i < possible_bins.size(); ++i) {
		// Cycle over the YON blocks this bin have data mapping to
		for (uint32_t j = 0; j < possible_bins[i].size(); ++j) {
			if (this->mImpl->linear_[possible_bins[i][j]].min_position < end_pos &&
			   this->mImpl->linear_[possible_bins[i][j]].max_position > start_pos)
			{
				yon_blocks.push_back(possible_bins[i][j]);
			}
		}
	}

	// Return nothing if all empty
	if (yon_blocks.size() == 0)
		return(std::vector<entry_type>());

	// Sort to dedupe
	std::sort(yon_blocks.begin(), yon_blocks.end());

	// Dedupe
	std::vector<entry_type> yon_blocks_deduped;
	yon_blocks_deduped.push_back(this->mImpl->linear_[yon_blocks[0]]);

	for (uint32_t i = 1; i < yon_blocks.size(); ++i) {
		if (yon_blocks[i] != yon_blocks_deduped.back().block_id) {
			yon_blocks_deduped.push_back(this->mImpl->linear_[yon_blocks[i]]);
		}
	}

	// Debug
	/*
	std::cerr << "Final\n" << std::endl;
	for (uint32_t i = 0; i < yon_blocks_deduped.size(); ++i) {
		yon_blocks_deduped[i].print(std::cerr);
		std::cerr << std::endl;
	}
	*/

	return(yon_blocks_deduped);
}

void yon_index_t::Setup(const std::vector<VcfContig>& contigs) {
	this->mImpl->index_.Add(contigs);
	this->mImpl->index_meta_.reserve(contigs.size());

	for (int i = 0; i < contigs.size(); ++i)
		this->mImpl->index_meta_[i].contig_id = contigs[i].idx;
}

void yon_index_t::Setup(const std::vector<YonContig>& contigs) {
	this->mImpl->index_.Add(contigs);
	this->mImpl->index_meta_.reserve(contigs.size());

	for (int i = 0; i < contigs.size(); ++i)
		this->mImpl->index_meta_[i].contig_id = contigs[i].idx;
}

int32_t yon_index_t::AddSorted(const uint32_t contig_id,
						 const uint64_t from_position,
						 const uint64_t to_position,
						 const uint32_t yon_block_id)
{
	return(this->mImpl->index_[contig_id].Add(from_position, to_position, yon_block_id));
}

yon_index_t& yon_index_t::operator+=(const entry_type& entry) {
	this->mImpl->linear_ += entry;

	// Update the meta index if the input data is sorted.
	if (this->is_sorted_) {
		if (this->mImpl->index_meta_[entry.contig_id].n_variants == 0) {
			this->mImpl->index_meta_[entry.contig_id].byte_offset_begin = entry.byte_offset;
			this->mImpl->index_meta_[entry.contig_id].min_position = entry.min_position;
			this->mImpl->index_meta_[entry.contig_id].start_block  = entry.block_id;
		}

		this->mImpl->index_meta_[entry.contig_id] += entry;
	}
	return(*this);
}

yon1_idx_rec& yon_index_t::operator[](const uint32_t block_id) { return(this->mImpl->linear_.at(block_id)); }
const yon1_idx_rec& yon_index_t::operator[](const uint32_t block_id) const { return(this->mImpl->linear_.at(block_id)); }


std::ostream& yon_index_t::Print(std::ostream& stream) const {
	for (int i = 0; i < this->mImpl->index_meta_.size(); ++i) {
		stream << "contig " << i << ". blocks: ";
		uint32_t n_blocks = 0;
		for (int j = 0; j < this->mImpl->index_[i].size(); ++j) {
			n_blocks += this->mImpl->index_[i][j].size();
		}
		stream << n_blocks;
		this->mImpl->index_meta_[i].Print(stream);
		stream << std::endl;
	}
	return(stream);
}

bool yon_index_t::IndexContainer(const yon1_vc_t& vc, const uint32_t block_id) {
	current_entry_.reset();
	current_entry_.block_id     = block_id;
	current_entry_.contig_id    = vc.block_.header.contig_id;
	current_entry_.min_position = vc.block_.header.min_position;
	current_entry_.max_position = vc.block_.header.max_position;
	current_entry_.n_variants   = vc.block_.header.n_variants;

	for (int i = 0; i < vc.size(); ++i) {
		if (this->IndexRecord(vc[i], block_id) == false) {
			std::cerr << "failed to index a record" << std::endl;
		}
	}

	return true;
}

bool yon_index_t::IndexRecord(const yon1_vnt_t& rcd, const uint32_t block_id) {
	int32_t index_bin = -1;
	const int32_t info_end_key = rcd.GetInfoOffset("END");

	// Ascertain that the meta entry has been evaluated
	// prior to executing this function.
	if (rcd.n_alleles == 0) {
		std::cerr << utility::timestamp("ERROR","IMPORT") << "The target meta record must be parsed prior to executing indexing functions..." << std::endl;
		return false;
	}

	int64_t end_position_used = rcd.pos;

	// The Info field END is used as the end position of an internal if it is available. This field
	// is usually only set for non-standard variants such as SVs or other special meaning records.
	if (info_end_key != -1) {
		if (rcd.info[info_end_key]->size() == 1) {
			const uint8_t word_width = rcd.info[info_end_key]->GetWordWidth();
			int64_t end_pos = 0;

			switch(word_width) {
			case(1): end_pos = reinterpret_cast<PrimitiveContainer<uint8_t>*>(rcd.info[info_end_key])->at(0);  break;
			case(2): end_pos = reinterpret_cast<PrimitiveContainer<uint16_t>*>(rcd.info[info_end_key])->at(0); break;
			case(4): end_pos = reinterpret_cast<PrimitiveContainer<uint32_t>*>(rcd.info[info_end_key])->at(0); break;
			case(8): end_pos = reinterpret_cast<PrimitiveContainer<uint64_t>*>(rcd.info[info_end_key])->at(0); break;
			default:
				std::cerr << utility::timestamp("ERROR", "INDEX") << "Unknown Info:End type: " << (int)word_width << std::endl;
				exit(1);
				break;
			}

			//std::cerr << "have end. " << "adding: (" << rcd.rid << "," << rcd.pos << "," << end_pos << ")"  << std::endl;

			index_bin = this->AddSorted(rcd.rid, rcd.pos, end_pos, block_id);

		} else {
			std::cerr << utility::timestamp("ERROR", "INDEX") << "Size of Info:End field is not 1..." << std::endl;
			index_bin = -1;
		}
	}

	// If the END field cannot be found then we check if the variant is a
	if (index_bin == -1) {
		int32_t longest = -1;
		// Iterate over available allele information and find the longest
		// SNV/indel length. The regex pattern ^[ATGC]{1,}$ searches for
		// simple SNV/indels.
		for (uint32_t i = 0; i < rcd.n_alleles; ++i) {
			if (std::regex_match(rcd.alleles[i].allele, YON_REGEX_CANONICAL_BASES)) {
				if (rcd.alleles[i].l_allele > longest)
					longest = rcd.alleles[i].l_allele;
			}
		}

		// Update the variant index with the target bin(s) found.
		if (longest > 1) {
			index_bin = this->AddSorted(rcd.rid, rcd.pos, rcd.pos + longest, block_id);
			//index_bin = 0;
			end_position_used = rcd.pos + longest;
		}
		// In the cases of special-meaning alleles such as copy-number (e.g. <CN>)
		// or SV (e.g. A[B)) they are index according to their left-most value only.
		// This has the implication that they cannot be found by means of interval
		// intersection searches. If special-meaning variants were to be supported
		// in the index then many more blocks would have to be searched for each
		// query as the few special cases will dominate the many general cases. For
		// this reason special-meaning alleles are not completely indexed.
		else {
			index_bin = this->AddSorted(rcd.rid, rcd.pos, rcd.pos, block_id);
		}
	}

	if (index_bin > this->current_entry_.max_bin) this->current_entry_.max_bin = index_bin;
	if (index_bin < this->current_entry_.min_bin) this->current_entry_.min_bin = index_bin;
	if (end_position_used > this->current_entry_.max_position)
		this->current_entry_.max_position = end_position_used;

	return true;
}

bool yon_index_t::IndexRecord(const bcf1_t*  record,
				 const uint32_t block_id,
				 const int32_t  info_end_key,
				 const yon1_vnt_t& rcd)
{
	assert(record != nullptr);
	int32_t index_bin = -1;

	// Ascertain that the meta entry has been evaluated
	// prior to executing this function.
	if (rcd.n_alleles == 0) {
		std::cerr << utility::timestamp("ERROR","IMPORT") << "The target meta record must be parsed prior to executing indexing functions..." << std::endl;
		return false;
	}

	int64_t end_position_used = record->pos;

	// The Info field END is used as the end position of an internal if it is available. This field
	// is usually only set for non-standard variants such as SVs or other special meaning records.
	if (info_end_key != -1) {
		// Linear search for the END key: this is not optimal but is probably faster
		// than first constructing a hash table for each record.
		const int n_info_fields = record->n_info;

		// Iterate over available Info fields.
		for (uint32_t i = 0; i < n_info_fields; ++i) {
			if (record->d.info[i].key == info_end_key) {
				uint32_t end = 0;

				switch(record->d.info[i].type) {
				case(BCF_BT_INT8):  end = *reinterpret_cast<int8_t*> (record->d.info[i].vptr); break;
				case(BCF_BT_INT16): end = *reinterpret_cast<int16_t*>(record->d.info[i].vptr); break;
				case(BCF_BT_INT32): end = *reinterpret_cast<int32_t*>(record->d.info[i].vptr); break;
				default:
					std::cerr << utility::timestamp("ERROR","INDEX") << "Illegal END primitive type: " << record->d.info[i].type << std::endl;
					return false;
				}

				index_bin = this->AddSorted(rcd.rid, record->pos, end, block_id);
				break;
			}
		}
	}

	// If the END field cannot be found then we check if the variant is a
	if (index_bin == -1) {
		int32_t longest = -1;
		// Iterate over available allele information and find the longest
		// SNV/indel length. The regex pattern ^[ATGC]{1,}$ searches for
		// simple SNV/indels.
		for (uint32_t i = 0; i < rcd.n_alleles; ++i) {
			if (std::regex_match(rcd.alleles[i].allele, YON_REGEX_CANONICAL_BASES)) {
				if (rcd.alleles[i].l_allele > longest)
					longest = rcd.alleles[i].l_allele;
			}
		}

		// Update the variant index with the target bin(s) found.
		if (longest > 1) {
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
			index_bin = this->AddSorted(rcd.rid, record->pos, record->pos, block_id);
		}
	}

	if (index_bin > this->current_entry_.max_bin) this->current_entry_.max_bin = index_bin;
	if (index_bin < this->current_entry_.min_bin) this->current_entry_.min_bin = index_bin;
	if (end_position_used > this->current_entry_.max_position)
		this->current_entry_.max_position = end_position_used;

	// Update number of entries in block
	++this->current_entry_.n_variants;

	return true;
}

std::ostream& operator<<(std::ostream& stream, const yon_index_t& entry) {
	utility::SerializePrimitive(entry.is_sorted_, stream);
	stream << entry.mImpl->index_;
	stream << entry.mImpl->index_meta_;
	stream << entry.mImpl->linear_;

	return(stream);
}

std::istream& operator>>(std::istream& stream, yon_index_t& entry) {
	utility::DeserializePrimitive(entry.is_sorted_, stream);
	stream >> entry.mImpl->index_;
	stream >> entry.mImpl->index_meta_;
	stream >> entry.mImpl->linear_;

	return(stream);
}

}
