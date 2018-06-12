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
	if(this->isBiallelicSNV())
		this->controller.simple_snv = true;

	for(U32 i = 0; i < this->n_alleles; ++i){
		new( &this->alleles[i] ) allele_type( );
		this->alleles[i](bcf_entry.alleles[i].data, bcf_entry.alleles[i].length);
	}
}

MetaEntry::MetaEntry(const bcf_entry_type& bcf_entry, const U64 position_offset) :
	n_alleles(bcf_entry.body->n_allele),
	//n_alleles(0),
	info_pattern_id(-1),
	filter_pattern_id(-1),
	format_pattern_id(-1),
	quality(bcf_entry.body->QUAL),
	contigID(bcf_entry.body->CHROM),
	position(bcf_entry.body->POS - position_offset),
	name(bcf_entry.ID, bcf_entry.l_ID),
	alleles(static_cast<allele_type*>(::operator new[](this->n_alleles*sizeof(allele_type))))
{
	for(U32 i = 0; i < this->n_alleles; ++i){
		new( &this->alleles[i] ) allele_type( );
		this->alleles[i](bcf_entry.alleles[i].data, bcf_entry.alleles[i].length);
	}

	if(this->n_alleles == 2) this->controller.biallelic = true;
	this->controller.simple_snv = bcf_entry.isSimple();
	if(this->isBiallelicSNV())
		this->controller.simple_snv = true;
}

MetaEntry::MetaEntry(const self_type& other) :
	controller(other.controller),
	n_alleles(other.n_alleles),
	info_pattern_id(other.info_pattern_id),
	filter_pattern_id(other.filter_pattern_id),
	format_pattern_id(other.format_pattern_id),
	quality(other.quality),
	contigID(other.contigID),
	position(other.position),
	alleles(static_cast<allele_type*>(::operator new[](this->n_alleles*sizeof(allele_type))))
{
	for(U32 i = 0; i < this->n_alleles; ++i)
		new( &this->alleles[i] ) allele_type( other.alleles[i] );
}

MetaEntry::~MetaEntry(){
	for(std::size_t i = 0; i < this->n_alleles; ++i)
		(this->alleles + i)->~MetaAllele();

	::operator delete[](static_cast<void*>(this->alleles));
};

const bool MetaEntry::usePackedRefAlt(void) const{
	if(this->isBiallelic() == false || this->isDiploid() == false)
		return false;

	if(std::regex_match(std::string(this->alleles[0].allele, this->alleles[0].l_allele), constants::YON_REGEX_PACKED_ALLELES) &&
	   std::regex_match(std::string(this->alleles[1].allele, this->alleles[1].l_allele), constants::YON_REGEX_PACKED_ALLELES)){
		return true;
	}
	return false;
}

const BYTE MetaEntry::packRefAltByte(void) const{
	assert(this->usePackedRefAlt());
	BYTE ref_alt = 0; // start out with empty

	if(this->alleles[0].l_allele == 9 && strncmp(this->alleles[0].allele, "<NON_REF>", 9) == 0){
		ref_alt ^= constants::REF_ALT_NON_REF << 4;
	} else {
		switch(this->alleles[0].allele[0]){
		case 'A': ref_alt ^= constants::REF_ALT_A << 4; break;
		case 'T': ref_alt ^= constants::REF_ALT_T << 4; break;
		case 'G': ref_alt ^= constants::REF_ALT_G << 4; break;
		case 'C': ref_alt ^= constants::REF_ALT_C << 4; break;
		case 'N': ref_alt ^= constants::REF_ALT_N << 4; break;
		case '.': ref_alt ^= constants::REF_ALT_MISSING << 4; break;
		default:
			std::cerr << utility::timestamp("ERROR", "BCF") << "Illegal SNV reference..." << std::endl;
			std::cerr << std::string(this->alleles[0].allele , this->alleles[0].l_allele) << std::endl;
			std::cerr << std::string(this->alleles[1].allele , this->alleles[1].l_allele) << std::endl;
			exit(1);
		}
	}

	if(this->alleles[1].l_allele == 9 && strncmp(this->alleles[1].allele, "<NON_REF>", 9) == 0){
		ref_alt ^= constants::REF_ALT_NON_REF << 0;
	} else {
		switch(this->alleles[1].allele[0]){
		case 'A': ref_alt ^= constants::REF_ALT_A << 0; break;
		case 'T': ref_alt ^= constants::REF_ALT_T << 0; break;
		case 'G': ref_alt ^= constants::REF_ALT_G << 0; break;
		case 'C': ref_alt ^= constants::REF_ALT_C << 0; break;
		case 'N': ref_alt ^= constants::REF_ALT_N << 0; break;
		case '.': ref_alt ^= constants::REF_ALT_MISSING << 0; break;
		default:
			std::cerr << utility::timestamp("ERROR", "BCF") << "Illegal SNV alt..." << std::endl;
			exit(1);
		}
	}
	return(ref_alt);
}

}
}
