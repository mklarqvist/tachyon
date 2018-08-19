#include "variant_index.h"

namespace tachyon{
namespace index{

VariantIndex::VariantIndex() :
	n_contigs_(0),
	n_capacity_(1000),
	contigs_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)))),
	linear_(static_cast<linear_type*>(::operator new[](this->capacity()*sizeof(linear_type))))
{

}

VariantIndex::VariantIndex(const self_type& other) :
	n_contigs_(other.n_contigs_),
	n_capacity_(other.n_capacity_),
	contigs_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)))),
	linear_(static_cast<linear_type*>(::operator new[](this->capacity()*sizeof(linear_type))))
{
	for(uint32_t i = 0; i < this->size(); ++i){
		new( &this->contigs_[i] ) value_type( other.contigs_[i] );
		new( &this->linear_[i] )  linear_type( other.linear_[i] );
	}
}

VariantIndex::~VariantIndex(){
	for(std::size_t i = 0; i < this->size(); ++i){
		(this->contigs_ + i)->~VariantIndexContig();
		(this->linear_ + i)->~VariantIndexLinear();
	}

	::operator delete[](static_cast<void*>(this->contigs_));
	::operator delete[](static_cast<void*>(this->linear_));
}

void VariantIndex::resize(void){
	pointer temp = this->contigs_;
	linear_type* temp_linear = this->linear_;

	this->n_capacity_ *= 2;
	this->contigs_ = static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type)));
	this->linear_  = static_cast<linear_type*>(::operator new[](this->capacity()*sizeof(linear_type)));


	// Lift over values from old addresses
	for(uint32_t i = 0; i < this->size(); ++i){
		new( &this->contigs_[i] ) value_type( temp[i] );
		new( &this->linear_[i] )  linear_type( temp_linear[i] );
	}

	// Clear temp
	for(std::size_t i = 0; i < this->size(); ++i){
		(temp + i)->~VariantIndexContig();
		(temp_linear + i)->~VariantIndexLinear();
	}

	::operator delete[](static_cast<void*>(temp));
	::operator delete[](static_cast<void*>(temp_linear));
}

std::ostream& operator<<(std::ostream& stream, const VariantIndex& index){
	stream.write(reinterpret_cast<const char*>(&index.n_contigs_), sizeof(std::size_t));
	std::size_t n_items_written = 0;
	for(uint32_t i = 0; i < index.size(); ++i){
		if(index.contigs_[i].size_sites()) ++n_items_written;
	}
	stream.write(reinterpret_cast<const char*>(&n_items_written), sizeof(std::size_t));

	for(uint32_t i = 0; i < index.size(); ++i){
		// Write if contig[i] contains data
		if(index.contigs_[i].size_sites())
			stream << index.contigs_[i];
	}

	// Write linear index
	for(uint32_t i = 0; i < index.size(); ++i) stream << index.linear_[i];

	return(stream);
}

std::istream& operator>>(std::istream& stream, VariantIndex& index){
	// Clear old data
	if(index.size()){
		for(std::size_t i = 0; i < index.size(); ++i){
			(index.contigs_ + i)->~VariantIndexContig();
			(index.linear_ + i)->~VariantIndexLinear();
		}

		::operator delete[](static_cast<void*>(index.contigs_));
		::operator delete[](static_cast<void*>(index.linear_));
	}

	stream.read(reinterpret_cast<char*>(&index.n_contigs_), sizeof(std::size_t));
	index.n_capacity_ = index.size() + 64;
	std::size_t n_items_written = 0;
	stream.read(reinterpret_cast<char*>(&n_items_written), sizeof(std::size_t));

	// Allocate new data
	index.contigs_ = static_cast<VariantIndexContig*>(::operator new[](index.capacity()*sizeof(VariantIndexContig)));
	index.linear_ = static_cast<VariantIndexLinear*>(::operator new[](index.capacity()*sizeof(VariantIndexLinear)));
	for(uint32_t i = 0; i < index.size(); ++i) {
		new( &index.contigs_[i] ) VariantIndexContig( );
		new( &index.linear_[i] ) VariantIndexLinear( );
	}

	// Load data and update accordingly
	for(uint32_t i = 0; i < n_items_written; ++i){
		VariantIndexContig temp;
		stream >> temp; // Read
		index.at(temp.getContigID()) = temp;
	}

	// Load linear index
	for(uint32_t i = 0; i < index.size(); ++i) stream >> index.linear_[i];

	return(stream);
}

}
}
