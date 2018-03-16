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

void MetaEntry::toVCFString(std::ostream& dest, const header_type& header) const{
	dest.write(&header.getContig(this->contigID).name[0], header.getContig(this->contigID).name.size()) << '\t';
	dest << this->position + 1 << '\t';

	if(this->name.size() == 0) dest.put('.');
	else dest.write(&this->name[0], this->name.size());
	dest.put('\t');
	if(this->n_alleles){
		dest.write(this->alleles[0].allele, this->alleles[0].l_allele);
		//dest << this->alleles[0].l_allele;
		dest.put('\t');
		dest.write(this->alleles[1].allele, this->alleles[1].l_allele);
		for(U32 i = 2; i < this->n_alleles; ++i){
			dest.put(',');
			dest.write(this->alleles[i].allele, this->alleles[i].l_allele);
		}
	} else dest << ".\t.\t";

	if(std::isnan(this->quality)) dest << "\t.\t";
	else {
		dest << '\t' << this->quality << '\t';
	}
}

void MetaEntry::toVCFString(buffer_type& dest, const header_type& header) const{
	return;

	/*
	dest.Add(&header.getContig(this->hot.contigID).name[0], header.getContig(this->hot.contigID).name.size());
	dest += '\t';
	//dest += std::to_string(blockPos + this->hot.position + 1);

	if(dest.size() + 100 >= dest.capacity()){
		//std::cerr << "resizing: " << dest.capacity() << "->" << dest.capacity()*2 << std::endl;
		dest.resize(dest.capacity()*2);
	}
	assert(dest.size() + 100 < dest.capacity());
	int ret = sprintf(&dest.buffer[dest.size()], "%llu", this->hot.position + 1);
	dest.n_chars += ret;
	dest += '\t';

	// If we have cold meta
	if(this->loaded_cold){
		if(this->cold.n_ID == 0) dest += '.';
		else dest.Add(this->cold.ID, this->cold.n_ID);
		dest += '\t';
		if(this->hot.controller.biallelic && this->hot.controller.simple_snv){
			dest += this->hot.ref_alt.getRef();
			dest += '\t';
			dest += this->hot.ref_alt.getAlt();
		}
		else {
			dest.Add(this->cold.alleles[0].allele, this->cold.alleles[0].l_allele);
			dest += '\t';
			U16 j = 1;
			for(; j < this->cold.n_allele - 1; ++j){
				dest.Add(this->cold.alleles[j].allele, this->cold.alleles[j].l_allele);
				dest += ',';
			}
			dest.Add(this->cold.alleles[j].allele, this->cold.alleles[j].l_allele);
		}
		if(std::isnan(this->cold.QUAL)) dest += "\t.\t";
		else {
			dest += '\t';
			ret = sprintf(&dest.buffer[dest.size()], "%g", this->cold.QUAL);
			dest.n_chars += ret;
			dest += '\t';
		}
	}
	*/
}

}
}
