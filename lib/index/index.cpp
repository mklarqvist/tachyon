#include <iostream>

#include "index.h"
#include "variant_index_meta.h"
#include "variant_index_quad_tree.h"
#include "variant_index_linear.h"

namespace tachyon {

class Index::IndexImpl {
public:
	typedef index::VariantIndexQuadTree  variant_quad_tree_type;
	typedef index::VariantIndexMeta      variant_meta_type;
	typedef index::VariantIndexMetaEntry entry_meta_type;
	typedef index::VariantIndexBin       bin_type;
	typedef index::VariantIndexLinear    variant_linear_type;

public:
	IndexImpl(){}
	IndexImpl(const IndexImpl& other){
		//index_meta_ += other.index_meta_;
		index_  += other.index_;
		//linear_ += other.linear_;
	}
	~IndexImpl(){}

	variant_meta_type      index_meta_;
	variant_quad_tree_type index_;
	variant_linear_type    linear_;
};

Index::Index() : mImpl(new IndexImpl), is_sorted_(true){}

Index::Index(const self_type& other) :
	is_sorted_(other.is_sorted_),
	mImpl(new IndexImpl(*other.mImpl))
{
}

Index& Index::operator+=(const Index& other){
	this->mImpl->index_ += other.mImpl->index_;
	return(*this);
}

Index::~Index(){}

bool Index::empty(void) const{ return(this->mImpl->index_.empty()); }
size_t Index::size(void) const{ return(this->mImpl->index_.size()); }
size_t Index::sizeMeta(void) const{ return(this->mImpl->index_meta_.size()); }
uint64_t Index::GetLinearSize(void) const{ return(this->mImpl->linear_.size()); }


std::vector<index::VariantIndexEntry> Index::FindOverlap(const uint32_t& contig_id) const{
	if(contig_id > this->mImpl->index_meta_.size())
		return(std::vector<entry_type>());

	std::vector<entry_type> yon_blocks;
	for(int i = this->mImpl->index_meta_[contig_id].start_block; i <= this->mImpl->index_meta_[contig_id].end_block; ++i){
		yon_blocks.push_back(this->mImpl->linear_[i]);
	}

	return(yon_blocks);
}

std::vector<index::VariantIndexEntry> Index::FindOverlap(const uint32_t& contig_id,
                                                  const uint64_t& start_pos,
                                                  const uint64_t& end_pos) const
{
	if(contig_id > this->mImpl->index_meta_.size())
		return(std::vector<entry_type>());

	if(this->mImpl->index_meta_.at(contig_id).n_blocks == 0)
		return(std::vector<entry_type>());


	// We also need to know possible overlaps in the quad-tree:
	// Seek from root to origin in quad-tree for potential overlapping bins with counts > 0.
	// Retrieve vector of bins that might contain the data of interest.
	// The PossibleBins function does not check if these bins are populated.
	std::vector<IndexImpl::bin_type> possible_bins = this->mImpl->index_[contig_id].PossibleBins(start_pos, end_pos);
	std::vector<uint32_t> yon_blocks;

	// Check if possible bins exists in the linear index
	for(int i = 0; i < possible_bins.size(); ++i){
		// Cycle over the YON blocks this bin have data mapping to
		for(uint32_t j = 0; j < possible_bins[i].size(); ++j){
			if(this->mImpl->linear_[possible_bins[i][j]].min_position < end_pos &&
			   this->mImpl->linear_[possible_bins[i][j]].max_position > start_pos)
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
	yon_blocks_deduped.push_back(this->mImpl->linear_[yon_blocks[0]]);

	for(uint32_t i = 1; i < yon_blocks.size(); ++i){
		if(yon_blocks[i] != yon_blocks_deduped.back().block_id){
			yon_blocks_deduped.push_back(this->mImpl->linear_[yon_blocks[i]]);
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

void Index::Setup(const std::vector<io::VcfContig>& contigs){
	this->mImpl->index_.Add(contigs);
	this->mImpl->index_meta_.reserve(contigs.size());
	for(int i = 0; i < contigs.size(); ++i)
		this->mImpl->index_meta_[i].contig_id = contigs[i].idx;
}

void Index::Setup(const std::vector<YonContig>& contigs){
	this->mImpl->index_.Add(contigs);
	this->mImpl->index_meta_.reserve(contigs.size());
	for(int i = 0; i < contigs.size(); ++i)
		this->mImpl->index_meta_[i].contig_id = contigs[i].idx;
}

int32_t Index::AddSorted(const uint32_t contig_id,
						 const uint64_t from_position,
						 const uint64_t to_position,
						 const uint32_t yon_block_id)
{
	return(this->mImpl->index_[contig_id].Add(from_position, to_position, yon_block_id));
}

Index& Index::operator+=(const entry_type& entry){
	this->mImpl->linear_ += entry;

	// Update the meta index if the input data is sorted.
	if(this->is_sorted_){
		if(this->mImpl->index_meta_[entry.contig_id].n_variants == 0){
			this->mImpl->index_meta_[entry.contig_id].byte_offset_begin = entry.byte_offset;
			this->mImpl->index_meta_[entry.contig_id].min_position = entry.min_position;
			this->mImpl->index_meta_[entry.contig_id].start_block  = entry.block_id;
		}

		this->mImpl->index_meta_[entry.contig_id] += entry;
	}
	return(*this);
}

index::VariantIndexEntry& Index::operator[](const uint32_t block_id){ return(this->mImpl->linear_.at(block_id)); }
const index::VariantIndexEntry& Index::operator[](const uint32_t block_id) const{ return(this->mImpl->linear_.at(block_id)); }


std::ostream& Index::Print(std::ostream& stream) const{
	for(int i = 0; i < this->mImpl->index_meta_.size(); ++i){
		stream << "contig " << i << ". blocks: ";
		uint32_t n_blocks = 0;
		for(int j = 0; j < this->mImpl->index_[i].size(); ++j){
			n_blocks += this->mImpl->index_[i][j].size();
		}
		stream << n_blocks;
		this->mImpl->index_meta_[i].Print(stream);
		stream << std::endl;
	}
	return(stream);
}

std::ostream& operator<<(std::ostream& stream, const Index& entry){
	utility::SerializePrimitive(entry.is_sorted_, stream);
	stream << entry.mImpl->index_;
	stream << entry.mImpl->index_meta_;
	stream << entry.mImpl->linear_;

	return(stream);
}

std::istream& operator>>(std::istream& stream, Index& entry){
	utility::DeserializePrimitive(entry.is_sorted_, stream);
	stream >> entry.mImpl->index_;
	stream >> entry.mImpl->index_meta_;
	stream >> entry.mImpl->linear_;

	return(stream);
}

}
