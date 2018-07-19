#include "variant_header.h"

namespace tachyon{
namespace core{

VariantHeader::VariantHeader(void) :
	contigs(nullptr),
	samples(nullptr),
	info_fields(nullptr),
	format_fields(nullptr),
	filter_fields(nullptr),
	htable_contigs(nullptr),
	htable_samples(nullptr),
	htable_info_fields(nullptr),
	htable_format_fields(nullptr),
	htable_filter_fields(nullptr)
{}

VariantHeader::VariantHeader(const vcf_header_type& vcf_header) :
	contigs(nullptr),
	samples(nullptr),
	info_fields(nullptr),
	format_fields(nullptr),
	filter_fields(nullptr),
	htable_contigs(nullptr),
	htable_samples(nullptr),
	htable_info_fields(nullptr),
	htable_format_fields(nullptr),
	htable_filter_fields(nullptr)
{
	// Invoke copy operator
	*this = vcf_header;
}

VariantHeader::VariantHeader(const self_type& other) :
	header_magic(other.header_magic),
	contigs(new contig_type[this->header_magic.getNumberContigs()]),
	samples(new sample_type[this->header_magic.getNumberSamples()]),
	info_fields(new map_entry_type[this->header_magic.n_info_values]),
	format_fields(new map_entry_type[this->header_magic.n_format_values]),
	filter_fields(new map_entry_type[this->header_magic.n_filter_values]),
	htable_contigs(nullptr),
	htable_samples(nullptr),
	htable_info_fields(nullptr),
	htable_format_fields(nullptr),
	htable_filter_fields(nullptr)
{
	for(U32 i = 0; i < this->header_magic.getNumberContigs(); ++i)
		this->contigs[i] = other.contigs[i];

	for(U32 i = 0; i < this->header_magic.getNumberSamples(); ++i)
		this->samples[i] = other.samples[i];

	for(U32 i = 0; i < this->header_magic.n_info_values; ++i){
		this->info_fields[i] = other.info_fields[i];
	}

	for(U32 i = 0; i < this->header_magic.n_format_values; ++i){
		this->format_fields[i] = other.format_fields[i];
	}

	for(U32 i = 0; i < this->header_magic.n_filter_values; ++i){
		this->filter_fields[i] = other.filter_fields[i];
	}

	this->buildHashTables();
}

VariantHeader::~VariantHeader(){
	delete [] this->contigs;
	delete [] this->samples;
	delete [] this->info_fields;
	delete [] this->format_fields;
	delete [] this->filter_fields;
	delete this->htable_contigs;
	delete this->htable_samples;
	delete this->htable_info_fields;
	delete this->htable_format_fields;
	delete this->htable_filter_fields;
}

bool VariantHeader::getContig(const std::string& p, contig_type*& target) const{
	if(this->htable_contigs == nullptr) return false;
	S32* ret = nullptr;
	if(this->htable_contigs->GetItem(&p[0], &p, ret, p.size())){
		target = &this->contigs[*ret];
		return true;
	}
	return false;
}

bool VariantHeader::getSample(const std::string& p, sample_type*& target) const{
	if(this->htable_samples == nullptr) return false;
	S32* ret = nullptr;
	if(this->htable_samples->GetItem(&p[0], &p, ret, p.size())){
		target = &this->samples[*ret];
		return true;
	}
	return false;
}

bool VariantHeader::getInfoField(const std::string& p, map_entry_type*& target) const{
	if(this->htable_info_fields== nullptr) return false;
	S32* ret = nullptr;
	if(this->htable_info_fields->GetItem(&p[0], &p, ret, p.size())){
		target = &this->info_fields[*ret];
		return true;
	}
	return false;
}

bool VariantHeader::getFormatField(const std::string& p, map_entry_type*& target) const{
	if(this->htable_format_fields == nullptr) return false;
	S32* ret = nullptr;
	if(this->htable_format_fields->GetItem(&p[0], &p, ret, p.size())){
		target = &this->format_fields[*ret];
		return true;
	}
	return false;
}

bool VariantHeader::getFilterField(const std::string& p, map_entry_type*& target) const{
	if(this->htable_filter_fields == nullptr) return false;
	S32* ret = nullptr;
	if(this->htable_filter_fields->GetItem(&p[0], &p, ret, p.size())){
		target = &this->filter_fields[*ret];
		return true;
	}
	return false;
}

const core::HeaderMapEntry* VariantHeader::getInfoField(const std::string& p) const{
	if(this->htable_info_fields == nullptr) return nullptr;
	S32* ret = nullptr;
	if(this->htable_info_fields->GetItem(&p[0], &p, ret, p.size()))
		return(&this->info_fields[*ret]);

	return nullptr;
}

const core::HeaderMapEntry* VariantHeader::getFormatField(const std::string& p) const{
	if(this->htable_format_fields == nullptr) return nullptr;
	S32* ret = nullptr;
	if(this->htable_format_fields->GetItem(&p[0], &p, ret, p.size()))
		return(&this->format_fields[*ret]);

	return nullptr;
}

const core::HeaderMapEntry* VariantHeader::getFilterField(const std::string& p) const{
	if(this->htable_filter_fields == nullptr) return nullptr;
	S32* ret = nullptr;
	if(this->htable_filter_fields->GetItem(&p[0], &p, ret, p.size()))
		return(&this->filter_fields[*ret]);

	return nullptr;
}

bool VariantHeader::buildHashTables(void){
	if(this->header_magic.n_contigs){
		if(this->header_magic.n_contigs*2 < 5012){
			this->htable_contigs = new hash_table_type(5012);
		} else
			this->htable_contigs = new hash_table_type(this->header_magic.n_contigs*2);

		for(S32 i = 0; i < this->header_magic.n_contigs; ++i){
			this->htable_contigs->SetItem(&this->contigs[i].name[0], &this->contigs[i].name, i, this->contigs[i].name.size());
		}
	}

	if(this->header_magic.n_samples){
		if(this->header_magic.n_samples*2 < 5012){
			this->htable_samples = new hash_table_type(5012);
		} else
			this->htable_samples = new hash_table_type(this->header_magic.n_samples*2);

		for(S32 i = 0; i < this->header_magic.n_samples; ++i){
			this->htable_samples->SetItem(&this->samples[i].name[0], &this->samples[i].name, i, this->samples[i].name.size());
		}
	}

	if(this->header_magic.n_info_values){
		if(this->header_magic.n_info_values*2 < 5012){
			this->htable_info_fields = new hash_table_type(5012);
		} else
			this->htable_info_fields = new hash_table_type(this->header_magic.n_info_values*2);

		for(S32 i = 0; i < this->header_magic.n_info_values; ++i){
			this->htable_info_fields->SetItem(&this->info_fields[i].ID[0], &this->info_fields[i].ID, i, this->info_fields[i].ID.size());
		}
	}

	if(this->header_magic.n_format_values){
		if(this->header_magic.n_format_values*2 < 5012){
			this->htable_format_fields = new hash_table_type(5012);
		} else
			this->htable_format_fields = new hash_table_type(this->header_magic.n_format_values*2);

		for(S32 i = 0; i < this->header_magic.n_format_values; ++i){
			this->htable_format_fields->SetItem(&this->format_fields[i].ID[0], &this->format_fields[i].ID, i, this->format_fields[i].ID.size());
		}
	}

	if(this->header_magic.n_filter_values){
		if(this->header_magic.n_filter_values*2 < 5012){
			this->htable_filter_fields = new hash_table_type(5012);
		} else
			this->htable_filter_fields = new hash_table_type(this->header_magic.n_filter_values*2);

		for(S32 i = 0; i < this->header_magic.n_filter_values; ++i){
			this->htable_filter_fields->SetItem(&this->filter_fields[i].ID[0], &this->filter_fields[i].ID, i, this->filter_fields[i].ID.size());
		}
	}

	return true;
}

void VariantHeader::operator=(const vcf_header_type& vcf_header){
	this->header_magic.n_contigs       = vcf_header.contigs.size();
	this->header_magic.n_samples       = vcf_header.sampleNames.size();
	this->header_magic.n_info_values   = vcf_header.info_map.size();
	this->header_magic.n_format_values = vcf_header.format_map.size();
	this->header_magic.n_filter_values = vcf_header.filter_map.size();

	if(vcf_header.literal_lines.size()){
		this->literals += vcf_header.literal_lines[0];
		for(U32 i = 1; i < vcf_header.literal_lines.size(); ++i)
			this->literals += "\n" + vcf_header.literal_lines[i];
	}

	this->header_magic.l_literals     = this->literals.size();

	// Cleanup previous
	delete [] this->contigs;
	delete [] this->samples;
	delete [] this->info_fields;
	delete [] this->filter_fields;
	delete [] this->format_fields;

	this->contigs = new contig_type[this->header_magic.getNumberContigs()];
	for(U32 i = 0; i < this->header_magic.getNumberContigs(); ++i){
		this->contigs[i] = vcf_header.contigs[i];
		this->contigs[i].contigID = i;
	}

	this->samples = new sample_type[this->header_magic.getNumberSamples()];
	for(U32 i = 0; i < this->header_magic.getNumberSamples(); ++i)
		this->samples[i] = vcf_header.sampleNames[i];

	this->info_fields = new map_entry_type[this->header_magic.n_info_values];
	for(U32 i = 0; i < this->header_magic.n_info_values; ++i){
		this->info_fields[i] = vcf_header.info_map[i];
		this->info_fields[i].IDX = i; // update idx to local
	}

	this->format_fields = new map_entry_type[this->header_magic.n_format_values];
	for(U32 i = 0; i < this->header_magic.n_format_values; ++i){
		this->format_fields[i] = vcf_header.format_map[i];
		this->format_fields[i].IDX = i; // update idx to local
	}

	this->filter_fields = new map_entry_type[this->header_magic.n_filter_values];
	for(U32 i = 0; i < this->header_magic.n_filter_values; ++i){
		this->filter_fields[i] = vcf_header.filter_map[i];
		this->filter_fields[i].IDX = i; // update idx to local
	}

	this->buildHashTables();
}

std::ostream& VariantHeader::writeHeaderVCF(std::ostream& stream, const bool showFormat) const{
	stream << this->literals;
	if(this->literals.size()) stream.put('\n');
	if(showFormat){
		stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
		if(this->header_magic.n_samples){
			stream.put('\t');
			stream << this->samples[0].name;
			for(U32 i = 1; i < this->header_magic.n_samples; ++i){
				stream << "\t" << this->samples[i].name;
			}
		}
	} else {
		stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	}
	stream.put('\n');
	return(stream);
}

std::ostream& VariantHeader::write(std::ostream& stream){
	stream << this->header_magic;
	for(U32 i = 0; i < this->header_magic.n_contigs; ++i) stream << this->contigs[i];
	for(U32 i = 0; i < this->header_magic.n_samples; ++i) stream << this->samples[i];
	for(U32 i = 0; i < this->header_magic.n_info_values; ++i)   stream << this->info_fields[i];
	for(U32 i = 0; i < this->header_magic.n_format_values; ++i) stream << this->format_fields[i];
	for(U32 i = 0; i < this->header_magic.n_filter_values; ++i) stream << this->filter_fields[i];

	stream.write(&this->literals[0], this->literals.size());
	return(stream);
}

}
}
