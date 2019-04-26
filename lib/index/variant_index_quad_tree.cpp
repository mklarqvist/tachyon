#include "variant_index_quad_tree.h"

namespace tachyon{
namespace index{

VariantIndexQuadTree::VariantIndexQuadTree() :
	n_contigs_(0),
	n_capacity_(1000),
	contigs_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{

}

VariantIndexQuadTree::VariantIndexQuadTree(const self_type& other) :
	n_contigs_(other.n_contigs_),
	n_capacity_(other.n_capacity_),
	contigs_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{
	for (uint32_t i = 0; i < this->size(); ++i) {
		new( &this->contigs_[i] ) value_type( other.contigs_[i] );
	}
}

VariantIndexQuadTree::~VariantIndexQuadTree() {
	for (std::size_t i = 0; i < this->size(); ++i) {
		(this->contigs_ + i)->~VariantIndexContig();
	}

	::operator delete[](static_cast<void*>(this->contigs_));
}

void VariantIndexQuadTree::resize(void) {
	pointer temp = this->contigs_;

	this->n_capacity_ *= 2;
	this->contigs_ = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));


	// Lift over values from old addresses
	for (uint32_t i = 0; i < this->size(); ++i) {
		new( &this->contigs_[i] ) value_type( temp[i] );
	}

	// Clear temp
	for (std::size_t i = 0; i < this->size(); ++i) {
		(temp + i)->~VariantIndexContig();
	}

	::operator delete[](static_cast<void*>(temp));
}

VariantIndexQuadTree& VariantIndexQuadTree::Add(const std::vector<VcfContig>& contigs) {
	while (this->size() + contigs.size() + 1 >= this->n_capacity_)
		this->resize();

	for (uint32_t i = 0; i < contigs.size(); ++i) {
		const uint64_t contig_length = contigs[i].n_bases;
		uint8_t n_levels = 7;
		// Safe-guard against using up excess memory in cases where no maximum
		// base-range is known.
		if(contigs[i].n_bases == std::numeric_limits<int32_t>::max())
			n_levels = 4;

		uint64_t bins_lowest = pow(4, n_levels);
		double used = ( bins_lowest - (contig_length % bins_lowest) ) + contig_length;

		if(used / bins_lowest < 2500) {
			for (int32_t i = n_levels; i != 0; --i) {
				if(used/pow(4,i) > 2500) {
					n_levels = i;
					break;
				}
			}
		}

		this->Add(i, contig_length, n_levels);
		//std::cerr << "contig: " << contigs[i].name << "(" << i << ")" << " -> " << contig_length << " levels: " << (int)n_levels << "->" << bins_lowest << "->" << used << std::endl;
		//std::cerr << "idx size:" << idx.size() << " at " << this->writer->index.variant_index_[i].size() << std::endl;
		//std::cerr << i << "->" << this->header->contigs[i].name << ":" << contig_length << " up to " << (uint64_t)used << " width (bp) lowest level: " << used/pow(4,n_levels) << "@level: " << (int)n_levels << std::endl;
	}
	return(*this);
}

VariantIndexQuadTree& VariantIndexQuadTree::Add(const std::vector<YonContig>& contigs) {
	while (this->size() + contigs.size() + 1 >= this->n_capacity_)
		this->resize();

	for (uint32_t i = 0; i < contigs.size(); ++i) {
		const uint64_t contig_length = contigs[i].n_bases;
		uint8_t n_levels = 7;
		uint64_t bins_lowest = pow(4,n_levels);
		double used = ( bins_lowest - (contig_length % bins_lowest) ) + contig_length;

		if(used / bins_lowest < 2500) {
			for (int32_t i = n_levels; i != 0; --i) {
				if(used/pow(4,i) > 2500) {
					n_levels = i;
					break;
				}
			}
		}

		this->Add(i, contig_length, n_levels);
		//std::cerr << "contig: " << this->header->contigs[i].name << "(" << i << ")" << " -> " << contig_length << " levels: " << (int)n_levels << std::endl;
		//std::cerr << "idx size:" << idx.size() << " at " << this->writer->index.variant_index_[i].size() << std::endl;
		//std::cerr << i << "->" << this->header->contigs[i].name << ":" << contig_length << " up to " << (uint64_t)used << " width (bp) lowest level: " << used/pow(4,n_levels) << "@level: " << (int)n_levels << std::endl;
	}
	return(*this);
}

std::ostream& operator<<(std::ostream& stream, const VariantIndexQuadTree& index) {
	stream.write(reinterpret_cast<const char*>(&index.n_contigs_), sizeof(std::size_t));
	std::size_t n_items_written = 0;
	for (uint32_t i = 0; i < index.size(); ++i) {
		if(index.contigs_[i].size_sites()) ++n_items_written;
	}
	stream.write(reinterpret_cast<const char*>(&n_items_written), sizeof(std::size_t));

	for (uint32_t i = 0; i < index.size(); ++i) {
		// Write if contig[i] contains data
		if(index.contigs_[i].size_sites())
			stream << index.contigs_[i];
	}

	return(stream);
}

std::istream& operator>>(std::istream& stream, VariantIndexQuadTree& index) {
	// Clear old data
	if(index.size()) {
		for (std::size_t i = 0; i < index.size(); ++i) {
			(index.contigs_ + i)->~VariantIndexContig();
		}

		::operator delete[](static_cast<void*>(index.contigs_));
	}

	stream.read(reinterpret_cast<char*>(&index.n_contigs_), sizeof(std::size_t));
	index.n_capacity_ = index.size() + 64;
	std::size_t n_items_written = 0;
	stream.read(reinterpret_cast<char*>(&n_items_written), sizeof(std::size_t));

	// Allocate new data
	index.contigs_ = static_cast<VariantIndexContig*>(::operator new[](index.capacity()*sizeof(VariantIndexContig)));
	for (uint32_t i = 0; i < index.size(); ++i) {
		new( &index.contigs_[i] ) VariantIndexContig( );
	}

	// Load data and update accordingly
	for (uint32_t i = 0; i < n_items_written; ++i) {
		VariantIndexContig temp;
		stream >> temp; // Read
		index.at(temp.GetContigID()) = temp;
	}

	return(stream);
}

}
}
