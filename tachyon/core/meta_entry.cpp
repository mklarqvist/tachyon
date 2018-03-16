#include "meta_entry.h"

namespace tachyon{
namespace core{

MetaEntry::MetaEntry() :
	n_alleles(0),
	info_pattern_id(-1),
	filter_pattern_id(-1),
	format_pattern_id(-1),
	quality(NAN),
	contigID(0),
	position(0),
	alleles(nullptr)
{}

MetaEntry::MetaEntry(const bcf_entry_type& bcf_entry) :
	n_alleles(bcf_entry.body->n_allele),
	info_pattern_id(-1),
	filter_pattern_id(-1),
	format_pattern_id(-1),
	quality(bcf_entry.body->QUAL),
	contigID(bcf_entry.body->CHROM),
	position(bcf_entry.body->POS),
	name(bcf_entry.ID, bcf_entry.l_ID),
	alleles(static_cast<allele_type*>(::operator new[](this->n_alleles*sizeof(allele_type))))
{
	if(this->n_alleles == 2) this->controller.biallelic = true;
	this->controller.simple_snv = bcf_entry.isSimple();
	if(this->isSimpleSNV() || this->isReferenceNONREF())
		this->controller.simple_snv = true;

	for(U32 i = 0; i < this->n_alleles; ++i){
		new( &this->alleles[i] ) allele_type( );
		this->alleles[i](bcf_entry.alleles[i].data, bcf_entry.alleles[i].length);
	}
}

MetaEntry::MetaEntry(const bcf_entry_type& bcf_entry, const U64 position_offset) :
	n_alleles(bcf_entry.body->n_allele),
	info_pattern_id(-1),
	filter_pattern_id(-1),
	format_pattern_id(-1),
	quality(bcf_entry.body->QUAL),
	contigID(bcf_entry.body->CHROM),
	position(bcf_entry.body->POS - position_offset),
	name(bcf_entry.ID, bcf_entry.l_ID),
	alleles(static_cast<allele_type*>(::operator new[](this->n_alleles*sizeof(allele_type))))
{
	if(this->n_alleles == 2) this->controller.biallelic = true;
	this->controller.simple_snv = bcf_entry.isSimple();
	if(this->isSimpleSNV() || this->isReferenceNONREF())
		this->controller.simple_snv = true;

	for(U32 i = 0; i < this->n_alleles; ++i){
		new( &this->alleles[i] ) allele_type( );
		this->alleles[i](bcf_entry.alleles[i].data, bcf_entry.alleles[i].length);
	}
}

MetaEntry::~MetaEntry(){
	for(std::size_t i = 0; i < this->n_alleles; ++i)
		(this->alleles + i)->~MetaAllele();

	::operator delete[](static_cast<void*>(this->alleles));
};

}
}
