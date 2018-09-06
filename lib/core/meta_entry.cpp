#include "meta_entry.h"

namespace tachyon{
namespace core{

MetaEntry::MetaEntry() :
	n_base_ploidy(0),
	n_alleles(0),
	info_pattern_id(-1),
	filter_pattern_id(-1),
	format_pattern_id(-1),
	quality(NAN),
	contigID(0),
	position(0),
	alleles(nullptr)
{}

MetaEntry::MetaEntry(const bcf1_t* record) :
	n_base_ploidy(0),
	n_alleles(record->n_allele),
	//n_alleles(0),
	info_pattern_id(-1),
	filter_pattern_id(-1),
	format_pattern_id(-1),
	quality(record->qual),
	contigID(record->rid),
	position(record->pos),
	name(record->d.id),
	alleles(static_cast<allele_type*>(::operator new[](this->n_alleles*sizeof(allele_type))))
{
	// Fix for the special case when ALT is not encoded
	if(this->n_alleles == 1){
		::operator delete[](static_cast<void*>(this->alleles));
		this->n_alleles = 2;
		this->alleles = static_cast<allele_type*>(::operator new[](this->n_alleles*sizeof(allele_type)));

		new( &this->alleles[0] ) allele_type( );
		this->alleles[0](std::string(record->d.allele[0]));
		new( &this->alleles[1] ) allele_type( );
		this->alleles[1].allele = new char[1];
		this->alleles[1].allele[0] = '.';
		this->alleles[1].l_allele = 1;
	} else {
		for(uint32_t i = 0; i < this->n_alleles; ++i){
			new( &this->alleles[i] ) allele_type( );
			this->alleles[i](std::string(record->d.allele[i]));
		}
	}

	if(this->n_alleles == 2) this->controller.biallelic = true;
	this->controller.simple_snv = (this->alleles[0].length() == 1 && this->alleles[1].length() == 1);
	if(this->IsBiallelicSNV())
		this->controller.simple_snv = true;
}

MetaEntry::MetaEntry(const bcf1_t* record, const uint64_t position_offset) :
	n_base_ploidy(0),
	n_alleles(record->n_allele),
	//n_alleles(0),
	info_pattern_id(-1),
	filter_pattern_id(-1),
	format_pattern_id(-1),
	quality(record->qual),
	contigID(record->rid),
	position(record->pos - position_offset),
	name(record->d.id),
	alleles(static_cast<allele_type*>(::operator new[](this->n_alleles*sizeof(allele_type))))
{
	// Fix for the special case when ALT is not encoded
	if(this->n_alleles == 1){
		::operator delete[](static_cast<void*>(this->alleles));
		this->n_alleles = 2;
		this->alleles = static_cast<allele_type*>(::operator new[](this->n_alleles*sizeof(allele_type)));

		new( &this->alleles[0] ) allele_type( );
		this->alleles[0](std::string(record->d.allele[0]));
		new( &this->alleles[1] ) allele_type( );
		this->alleles[1].allele = new char[1];
		this->alleles[1].allele[0] = '.';
		this->alleles[1].l_allele = 1;
	} else {
		for(uint32_t i = 0; i < this->n_alleles; ++i){
			new( &this->alleles[i] ) allele_type( );
			this->alleles[i](std::string(record->d.allele[i]));
		}
	}

	if(this->n_alleles == 2) this->controller.biallelic = true;
	this->controller.simple_snv = (this->alleles[0].length() == 1 && this->alleles[1].length() == 1);
	if(this->IsBiallelicSNV())
		this->controller.simple_snv = true;
}

MetaEntry::MetaEntry(const self_type& other) :
	n_base_ploidy(other.n_base_ploidy),
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
	for(uint32_t i = 0; i < this->n_alleles; ++i)
		new( &this->alleles[i] ) allele_type( other.alleles[i] );
}

MetaEntry::~MetaEntry(){
	for(std::size_t i = 0; i < this->n_alleles; ++i)
		(this->alleles + i)->~MetaAllele();

	::operator delete[](static_cast<void*>(this->alleles));
};

bool MetaEntry::UsePackedRefAlt(void) const{
	if(this->IsBiallelic() == false || this->IsDiploid() == false)
		return false;

	if(std::regex_match(std::string(this->alleles[0].allele, this->alleles[0].l_allele), constants::YON_REGEX_PACKED_ALLELES) &&
	   std::regex_match(std::string(this->alleles[1].allele, this->alleles[1].l_allele), constants::YON_REGEX_PACKED_ALLELES)){
		return true;
	}
	return false;
}

uint8_t MetaEntry::PackRefAltByte(void) const{
	assert(this->UsePackedRefAlt());
	uint8_t ref_alt = 0; // start out with empty

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

TACHYON_VARIANT_CLASSIFICATION_TYPE MetaEntry::ClassifyVariant(const uint32_t& allele) const {
	const int32_t ref_size = this->alleles[0].size();
	const int32_t l_diff   = ref_size - this->alleles[allele].size();

	if(this->alleles[0].allele[0] == '<' || this->alleles[allele].allele[0] == '<')
		return(YON_VARIANT_CLASS_SV);
	else if(l_diff == 0){
		if(ref_size == 1 && this->alleles[0].allele[0] != this->alleles[allele].allele[0]){
			if(this->alleles[allele].allele[0] == 'A' ||
			   this->alleles[allele].allele[0] == 'T' ||
			   this->alleles[allele].allele[0] == 'G' ||
			   this->alleles[allele].allele[0] == 'C')
			{
				return(YON_VARIANT_CLASS_SNP);
			}
			else return(YON_VARIANT_CLASS_UNKNOWN);
		}
		else if(ref_size != 1){
			uint32_t n_characters_identical = 0;
			const uint32_t length_shortest = ref_size < this->alleles[allele].size()
					                       ? ref_size
					                       : this->alleles[allele].size();

			for(uint32_t c = 0; c < length_shortest; ++c)
				n_characters_identical += (this->alleles[0].allele[c] == this->alleles[allele].allele[c]);

			if(n_characters_identical == 0) return(YON_VARIANT_CLASS_MNP);
			else return(YON_VARIANT_CLASS_CLUMPED);
		}
	} else {
		const uint32_t length_shortest = ref_size < this->alleles[allele].size()
		                               ? ref_size
		                               : this->alleles[allele].size();

		// Keep track of non-standard characters.
		uint32_t n_characters_non_standard = 0;

		// Iterate over available characters and check for non-standard
		// genetic characters (ATGC).
		for(uint32_t c = 0; c < length_shortest; ++c){
			n_characters_non_standard += (this->alleles[allele].allele[c] != 'A' &&
			                              this->alleles[allele].allele[c] != 'T' &&
			                              this->alleles[allele].allele[c] != 'C' &&
			                              this->alleles[allele].allele[c] != 'G');
		}

		// If non-standard characters are found then return as unknown
		// type. Otherwise, return classification as an indel.
		if(n_characters_non_standard) return(YON_VARIANT_CLASS_UNKNOWN);
		else return(YON_VARIANT_CLASS_INDEL);
	}
	return(YON_VARIANT_CLASS_UNKNOWN);
}

}
}
