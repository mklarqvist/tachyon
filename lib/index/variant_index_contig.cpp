#include <iostream>

#include "variant_index_contig.h"

namespace tachyon{
namespace index{

VariantIndexContig::VariantIndexContig() :
	contig_id(0),
	l_contig(0),
	l_contig_rounded(0),
	n_bins(0),
	n_capacity(0),
	n_levels(0),
	n_sites(0),
	bins_cumsum(nullptr),
	bins(nullptr)
{

}

VariantIndexContig::VariantIndexContig(const uint32_t contigID, const uint64_t l_contig, const uint8_t n_levels) :
	contig_id(contigID),
	l_contig(l_contig),
	l_contig_rounded(0),
	n_bins(0),
	n_capacity(0),
	n_levels(n_levels),
	n_sites(0),
	bins_cumsum(nullptr),
	bins(nullptr)
{
	this->l_contig_rounded = this->RoundLengthClosestBase4(this->l_contig);
	if (this->n_levels != 0) {
		this->CalculateCumulativeSums();
		this->n_capacity = this->bins_cumsum[this->n_levels] + 64;
		this->n_bins     = this->bins_cumsum[this->n_levels];
		this->bins       = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
		for (uint32_t i = 0; i < this->size(); ++i) {
			new( &this->bins[i] ) value_type(  );
			this->bins[i].binID_ = i;
		}
	}
}

VariantIndexContig::VariantIndexContig(const self_type& other) :
	contig_id(other.contig_id),
	l_contig(other.l_contig),
	l_contig_rounded(other.l_contig_rounded),
	n_bins(other.n_bins),
	n_capacity(other.n_capacity),
	n_levels(other.n_levels),
	n_sites(other.n_sites),
	bins_cumsum(nullptr),
	bins(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{
	this->CalculateCumulativeSums();
	for (uint32_t i = 0; i < this->size(); ++i)
		new( &this->bins[i] ) value_type( other.bins[i] );
}

void VariantIndexContig::operator=(const self_type& other) {
	// Clean previous
	for (std::size_t i = 0; i < this->size(); ++i)
		(this->bins + i)->~VariantIndexBin();

	::operator delete[](static_cast<void*>(this->bins));
	delete [] this->bins_cumsum;

	this->contig_id   = other.contig_id;
	this->l_contig    = other.l_contig;
	this->l_contig_rounded = other.l_contig_rounded;
	this->n_bins      = other.n_bins;
	this->n_capacity  = other.n_capacity;
	this->n_levels    = other.n_levels;
	this->bins_cumsum = nullptr;

	this->CalculateCumulativeSums();

	this->bins = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
	for (uint32_t i = 0; i < this->size(); ++i)
		new( &this->bins[i] ) value_type( other.bins[i] );
}


VariantIndexContig::~VariantIndexContig() {
	for (std::size_t i = 0; i < this->size(); ++i)
		(this->bins + i)->~VariantIndexBin();

	::operator delete[](static_cast<void*>(this->bins));
}

int32_t VariantIndexContig::Add(const uint64_t fromPosition, const uint64_t toPosition, const uint32_t yon_block_id) {
	for (int32_t i = this->n_levels; i != 0; --i) {
		uint32_t binFrom = int64_t(fromPosition/(this->l_contig_rounded / pow(4,i)));
		uint32_t binTo   = int64_t(toPosition/(this->l_contig_rounded / pow(4,i)));
		/**
		 * If both ends of the interval map into the same chunk we know the interval is
		 * completely contained: in this case we deposit the interval there
		 **/
		if (binFrom == binTo) {
			//if (i != this->n_levels) std::cerr << fromPosition << "->" << toPosition << ", adding to " << binFrom << " level " << i << " cum : " << this->binscumsum[i-1]+binFrom << "/" << this->size() << std::endl;
			++this->n_sites;
			this->bins[this->bins_cumsum[i - 1] + binFrom].Add(yon_block_id);
			return(this->bins_cumsum[i - 1] + binFrom);
		}
	}
	this->bins[0].Add(yon_block_id);
	++this->n_sites;
	return(0);
}

std::vector<VariantIndexBin> VariantIndexContig::PossibleBins(const uint64_t& from_position, const uint64_t& to_position, const bool filter) const {
	std::vector<value_type> overlapping_chunks;
	//overlapping_chunks.push_back(this->at(0)); // level 0

	// If end position are out-of-bounds then trucated it to maximum
	// allowed value
	uint64_t used_to_posititon = to_position;
	if (used_to_posititon > this->l_contig_rounded) {
		//std::cerr << "out of bounds" << std::endl;
		//std::cerr << to_position << "->" << this->l_contig_rounded << std::endl;
		used_to_posititon = this->l_contig_rounded;
	}

	for (int32_t i = this->n_levels; i != 0; --i) {
		int64_t binFrom = int64_t(from_position/(this->l_contig_rounded / pow(4,i)));
		int64_t binTo   = int64_t(used_to_posititon/(this->l_contig_rounded / pow(4,i)));

		//std::cerr << i << "/" << (int)this->n_levels << ": level offset " << this->binscumsum[i-1] << "; (from, to) " << binFrom << " -> " << binTo << " out of " << this->size() << std::endl;
		//std::cerr << "limit: " << this->binscumsum[i] << std::endl;

		// Overlap from cumpos + (binFrom, binTo)
		// All these chunks could potentially hold intervals overlapping
		// the desired coordinates
		for (uint32_t j = binFrom; j <= binTo; ++j) {
			if (filter == false)
				overlapping_chunks.push_back(this->at(this->bins_cumsum[i - 1] + j));
			else {
				if (this->at(this->bins_cumsum[i - 1] + j).size())
					overlapping_chunks.push_back(this->at(this->bins_cumsum[i - 1] + j));
			}
		}
	}
	overlapping_chunks.push_back(this->at(0));


	return(overlapping_chunks);
}

}
}
