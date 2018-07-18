#include "variant_index_contig.h"

namespace tachyon{
namespace index{

VariantIndexContig::VariantIndexContig() :
	contigID_(0),
	l_contig_(0),
	l_contig_rounded_(0),
	n_bins_(0),
	n_capacity_(0),
	n_levels_(0),
	n_sites_(0),
	bins_cumsum_(nullptr),
	bins_(nullptr)
{

}

VariantIndexContig::VariantIndexContig(const U32 contigID, const U64 l_contig, const BYTE n_levels) :
	contigID_(contigID),
	l_contig_(l_contig),
	l_contig_rounded_(0),
	n_bins_(0),
	n_capacity_(0),
	n_levels_(n_levels),
	n_sites_(0),
	bins_cumsum_(nullptr),
	bins_(nullptr)
{
	this->l_contig_rounded_ = this->roundLengthClosestBase4_(this->l_contig_);
	if(this->n_levels_ != 0){
		this->calculateCumulativeSums_();
		this->n_capacity_ = this->bins_cumsum_[this->n_levels_] + 64;
		this->n_bins_     = this->bins_cumsum_[this->n_levels_];
		this->bins_       = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
		for(U32 i = 0; i < this->size(); ++i){
			new( &this->bins_[i] ) value_type(  );
			this->bins_[i].binID_ = i;
		}
	}
}

VariantIndexContig::VariantIndexContig(const self_type& other) :
	contigID_(other.contigID_),
	l_contig_(other.l_contig_),
	l_contig_rounded_(other.l_contig_rounded_),
	n_bins_(other.n_bins_),
	n_capacity_(other.n_capacity_),
	n_levels_(other.n_levels_),
	n_sites_(other.n_sites_),
	bins_cumsum_(nullptr),
	bins_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{
	this->calculateCumulativeSums_();
	for(U32 i = 0; i < this->size(); ++i)
		new( &this->bins_[i] ) value_type( other.bins_[i] );
}

void VariantIndexContig::operator=(const self_type& other){
	// Clean previous
	for(std::size_t i = 0; i < this->size(); ++i)
		(this->bins_ + i)->~VariantIndexBin();

	::operator delete[](static_cast<void*>(this->bins_));
	delete [] this->bins_cumsum_;

	this->contigID_ = other.contigID_;
	this->l_contig_ = other.l_contig_;
	this->l_contig_rounded_ = other.l_contig_rounded_;
	this->n_bins_ = other.n_bins_;
	this->n_capacity_ = other.n_capacity_;
	this->n_levels_ = other.n_levels_;
	this->bins_cumsum_ = nullptr;
	this->calculateCumulativeSums_();

	this->bins_ = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
	for(U32 i = 0; i < this->size(); ++i)
		new( &this->bins_[i] ) value_type( other.bins_[i] );
}


VariantIndexContig::~VariantIndexContig(){
	for(std::size_t i = 0; i < this->size(); ++i)
		(this->bins_ + i)->~VariantIndexBin();

	::operator delete[](static_cast<void*>(this->bins_));
}

const S32 VariantIndexContig::Add(const U64& fromPosition, const U64& toPosition, const U32& yon_block_id){
	for(S32 i = this->n_levels_; i != 0; --i){
		U32 binFrom = S64(fromPosition/(this->l_contig_rounded_ / pow(4,i)));
		U32 binTo   = S64(toPosition/(this->l_contig_rounded_ / pow(4,i)));
		/**
		 * If both ends of the interval map into the same chunk we know the interval is
		 * completely contained: in this case we deposit the interval there
		 **/
		if(binFrom == binTo){
			//if(i != this->n_levels_) std::cerr << fromPosition << "->" << toPosition << ", adding to " << binFrom << " level " << i << " cum : " << this->bins_cumsum_[i-1]+binFrom << "/" << this->size() << std::endl;
			++this->n_sites_;
			this->bins_[this->bins_cumsum_[i - 1]+binFrom].Add(yon_block_id);
			return(this->bins_cumsum_[i - 1]+binFrom);
		}
	}
	this->bins_[0].Add(yon_block_id);
	++this->n_sites_;
	return(0);
}

std::vector<VariantIndexBin> VariantIndexContig::possibleBins(const U64& from_position, const U64& to_position, const bool filter) const{
	std::vector<value_type> overlapping_chunks;
	//overlapping_chunks.push_back(this->at(0)); // level 0

	// If end position are out-of-bounds then trucated it to maximum
	// allowed value
	U64 used_to_posititon = to_position;
	if(used_to_posititon > this->l_contig_rounded_){
		//std::cerr << "out of bounds" << std::endl;
		//std::cerr << to_position << "->" << this->l_contig_rounded_ << std::endl;
		used_to_posititon = this->l_contig_rounded_;
	}

	for(S32 i = this->n_levels_; i != 0; --i){
		S64 binFrom = S64(from_position/(this->l_contig_rounded_ / pow(4,i)));
		S64 binTo   = S64(used_to_posititon/(this->l_contig_rounded_ / pow(4,i)));

		//std::cerr << i << "/" << (int)this->n_levels_ << ": level offset " << this->bins_cumsum_[i-1] << "; (from, to) " << binFrom << " -> " << binTo << " out of " << this->size() << std::endl;
		//std::cerr << "limit: " << this->bins_cumsum_[i] << std::endl;

		// Overlap from cumpos + (binFrom, binTo)
		// All these chunks could potentially hold intervals overlapping
		// the desired coordinates
		for(U32 j = binFrom; j <= binTo; ++j){
			if(filter == false)
				overlapping_chunks.push_back(this->at(this->bins_cumsum_[i - 1] + j));
			else {
				if(this->at(this->bins_cumsum_[i - 1] + j).size())
					overlapping_chunks.push_back(this->at(this->bins_cumsum_[i - 1] + j));
			}
		}
	}
	overlapping_chunks.push_back(this->at(0));


	return(overlapping_chunks);
}

}
}
