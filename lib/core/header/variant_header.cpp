#include "variant_header.h"

namespace tachyon{

bool VariantHeader::BuildReverseMaps(void){
	this->contigs_reverse_map_.clear();
	this->info_fields_reverse_map_.clear();
	this->format_fields_reverse_map_.clear();
	this->filter_fields_reverse_map_.clear();

	for(uint32_t i = 0; i < this->contigs_.size(); ++i)       this->contigs_reverse_map_[this->contigs_[i].idx] = i;
	for(uint32_t i = 0; i < this->info_fields_.size(); ++i)   this->info_fields_reverse_map_[this->info_fields_[i].idx] = i;
	for(uint32_t i = 0; i < this->format_fields_.size(); ++i) this->format_fields_reverse_map_[this->format_fields_[i].idx] = i;
	for(uint32_t i = 0; i < this->filter_fields_.size(); ++i) this->filter_fields_reverse_map_[this->filter_fields_[i].idx] = i;

	return true;
}

bool VariantHeader::BuildMaps(void){
	this->info_fields_map_.clear();
	this->format_fields_map_.clear();
	this->filter_fields_map_.clear();
	this->contigs_map_.clear();

	for(uint32_t i = 0; i < this->contigs_.size(); ++i)       this->contigs_map_[this->contigs_[i].name] = i;
	for(uint32_t i = 0; i < this->info_fields_.size(); ++i)   this->info_fields_map_[this->info_fields_[i].id] = i;
	for(uint32_t i = 0; i < this->format_fields_.size(); ++i) this->format_fields_map_[this->format_fields_[i].id] = i;
	for(uint32_t i = 0; i < this->filter_fields_.size(); ++i) this->filter_fields_map_[this->filter_fields_[i].id] = i;
	for(uint32_t i = 0; i < this->samples_.size(); ++i)       this->samples_map_[this->samples_[i]] = i;

	return true;
}

bool VariantHeader::RecodeIndices(void){
	for(uint32_t i = 0; i < this->contigs_.size(); ++i)       this->contigs_[i].idx = i;
	for(uint32_t i = 0; i < this->info_fields_.size(); ++i)   this->info_fields_[i].idx = i;
	for(uint32_t i = 0; i < this->format_fields_.size(); ++i) this->format_fields_[i].idx = i;
	for(uint32_t i = 0; i < this->filter_fields_.size(); ++i) this->filter_fields_[i].idx = i;

	if(this->BuildMaps() == false) return false;
	if(this->BuildReverseMaps() == false) return false;
	return true;
}

bcf_hdr_t* VariantHeader::ConvertVcfHeaderLiterals(const bool add_format){
	std::string internal = this->literals_;
	internal += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if(this->samples_.size() && add_format){
		internal += "\tFORMAT\t";
		internal += this->samples_[0];
		for(size_t i = 1; i < this->samples_.size(); ++i)
			internal += "\t" + this->samples_[i];
	}
	internal += "\n";

	hts_vcf_header* hdr = bcf_hdr_init("r");
	int ret = bcf_hdr_parse(hdr, (char*)internal.c_str());
	if(ret != 0){
		std::cerr << "failed to get bcf header from literals" << std::endl;
		bcf_hdr_destroy(hdr);
		return(nullptr);
	}

	return(hdr);
}

bcf_hdr_t* VariantHeader::ConvertVcfHeader(const bool add_format){
	std::string internal = this->ToString(true);
	internal += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if(this->samples_.size() && add_format){
		internal += "\tFORMAT\t";
		internal += this->samples_[0];
		for(size_t i = 1; i < this->samples_.size(); ++i)
			internal += "\t" + this->samples_[i];
	}
	internal += "\n";

	hts_vcf_header* hdr = bcf_hdr_init("r");
	int ret = bcf_hdr_parse(hdr, (char*)internal.c_str());
	if(ret != 0){
		std::cerr << "failed to get bcf header from literals" << std::endl;
		bcf_hdr_destroy(hdr);
		return(nullptr);
	}

	return(hdr);
}

void VariantHeader::AddGenotypeAnnotationFields(void){
	//"NM","NPM","AN","HWE_P","AC","AF","AC_P","FS_A","F_PIC","HET","MULTI_ALLELIC"

	const YonInfo* info = this->GetInfo("NM");
	if(info == nullptr){
		YonInfo nm;
		nm.id = "NM";
		nm.number = "1";
		nm.type = "Integer";
		nm.yon_type = YON_VCF_HEADER_INTEGER;
		nm.description = "\"NM\"";
		nm.idx = this->info_fields_.size();
		this->literals_ += nm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(nm);
	}

	info = this->GetInfo("NPM");
	if(info == nullptr){
		YonInfo npm;
		npm.id = "NPM";
		npm.number = "1";
		npm.type = "Integer";
		npm.yon_type = YON_VCF_HEADER_INTEGER;
		npm.description = "\"NPM\"";
		npm.idx = this->info_fields_.size();
		this->literals_ += npm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(npm);
	}

	info = this->GetInfo("AN");
	if(info == nullptr){
		YonInfo npm;
		npm.id = "AN";
		npm.number = "1";
		npm.type = "Integer";
		npm.yon_type = YON_VCF_HEADER_INTEGER;
		npm.description = "\"Total number of alleles in called genotypes\"";
		npm.idx = this->info_fields_.size();
		this->literals_ += npm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(npm);
	}

	info = this->GetInfo("HWE_P");
	if(info == nullptr){
		YonInfo npm;
		npm.id = "HWE_P";
		npm.number = "1";
		npm.type = "Float";
		npm.yon_type = YON_VCF_HEADER_FLOAT;
		npm.description = "\"Hardy-Weinberg equilibrium P-value\"";
		npm.idx = this->info_fields_.size();
		this->literals_ += npm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(npm);
	}

	info = this->GetInfo("AC");
	if(info == nullptr){
		YonInfo npm;
		npm.id = "AC";
		npm.number = "A";
		npm.type = "Integer";
		npm.yon_type = YON_VCF_HEADER_INTEGER;
		npm.description = "\"Allele count in genotypes, for each REF and ALT allele, in the same order as listed with REF first\"";
		npm.idx = this->info_fields_.size();
		this->literals_ += npm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(npm);
	}

	info = this->GetInfo("AF");
	if(info == nullptr){
		YonInfo npm;
		npm.id = "AF";
		npm.number = "A";
		npm.type = "Float";
		npm.yon_type = YON_VCF_HEADER_FLOAT;
		npm.description = "\"Allele Frequency, for each REF and ALT allele, in the same order as listed with REF first. In the range [0, 1]\"";
		npm.idx = this->info_fields_.size();
		this->literals_ += npm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(npm);
	}

	info = this->GetInfo("AC_P");
	if(info == nullptr){
		YonInfo npm;
		npm.id = "AC_P";
		npm.number = "A";
		npm.type = "Integer";
		npm.yon_type = YON_VCF_HEADER_INTEGER;
		npm.description = "\"AC_P\"";
		npm.idx = this->info_fields_.size();
		this->literals_ += npm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(npm);
	}

	info = this->GetInfo("FS_A");
	if(info == nullptr){
		YonInfo npm;
		npm.id = "FS_A";
		npm.number = "A";
		npm.type = "Float";
		npm.yon_type = YON_VCF_HEADER_FLOAT;
		npm.description = "\"FS_A\"";
		npm.idx = this->info_fields_.size();
		this->literals_ += npm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(npm);
	}

	info = this->GetInfo("F_PIC");
	if(info == nullptr){
		YonInfo nm;
		nm.id = "F_PIC";
		nm.number = "1";
		nm.type = "Float";
		nm.yon_type = YON_VCF_HEADER_FLOAT;
		nm.description = "\"F_PIC\"";
		nm.idx = this->info_fields_.size();
		this->literals_ += nm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(nm);
	}

	info = this->GetInfo("HET");
	if(info == nullptr){
		YonInfo nm;
		nm.id = "HET";
		nm.number = "1";
		nm.type = "Float";
		nm.yon_type = YON_VCF_HEADER_FLOAT;
		nm.description = "\"Heterozygosity at this locus calculated as number of 0/1 or 1/0 genotypes divided by all non-missing genotypes\"";
		nm.idx = this->info_fields_.size();
		this->literals_ += nm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(nm);
	}

	info = this->GetInfo("MULTI_ALLELIC");
	if(info == nullptr){
		YonInfo nm;
		nm.id = "MULTI_ALLELIC";
		nm.number = "1";
		nm.type = "Flag";
		nm.yon_type = YON_VCF_HEADER_FLAG;
		nm.description = "\"Flag indicating if a site is multi-allelic (>1 ALT alleles)\"";
		nm.idx = this->info_fields_.size();
		this->literals_ += nm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(nm);
	}

	// temp
	info = this->GetInfo("FIS");
	if(info == nullptr){
		YonInfo nm;
		nm.id = "FIS";
		nm.number = "1";
		nm.type = "Float";
		nm.yon_type = YON_VCF_HEADER_FLOAT;
		nm.description = "\"FIS\"";
		nm.idx = this->info_fields_.size();
		this->literals_ += nm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(nm);
	}

	info = this->GetInfo("FST");
	if(info == nullptr){
		YonInfo nm;
		nm.id = "FST";
		nm.number = "1";
		nm.type = "Float";
		nm.yon_type = YON_VCF_HEADER_FLOAT;
		nm.description = "\"FST\"";
		nm.idx = this->info_fields_.size();
		this->literals_ += nm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(nm);
	}

	info = this->GetInfo("FIT");
	if(info == nullptr){
		YonInfo nm;
		nm.id = "FIT";
		nm.number = "1";
		nm.type = "Float";
		nm.yon_type = YON_VCF_HEADER_FLOAT;
		nm.description = "\"FIT\"";
		nm.idx = this->info_fields_.size();
		this->literals_ += nm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(nm);
	}

	this->BuildMaps();
	this->BuildReverseMaps();
}


void VariantHeader::AddGenotypeAnnotationFields(const std::vector<std::string>& group_names){
	for(int i = 0; i < group_names.size(); ++i){
		const YonInfo* info = this->GetInfo(group_names[i] + "_NM");
		if(info == nullptr){
			YonInfo nm;
			nm.id = group_names[i] + "_NM";
			nm.number = "1";
			nm.type = "Integer";
			nm.yon_type = YON_VCF_HEADER_INTEGER;
			nm.description = "\"NM\"";
			nm.idx = this->info_fields_.size();
			this->literals_ += nm.ToVcfString(false) + "\n";
			this->info_fields_.push_back(nm);
		}

		info = this->GetInfo(group_names[i] + "_NPM");
		if(info == nullptr){
			YonInfo npm;
			npm.id = group_names[i] + "_NPM";
			npm.number = "1";
			npm.type = "Integer";
			npm.yon_type = YON_VCF_HEADER_INTEGER;
			npm.description = "\"NPM\"";
			npm.idx = this->info_fields_.size();
			this->literals_ += npm.ToVcfString(false) + "\n";
			this->info_fields_.push_back(npm);
		}

		info = this->GetInfo(group_names[i] + "_AN");
		if(info == nullptr){
			YonInfo npm;
			npm.id = group_names[i] + "_AN";
			npm.number = "1";
			npm.type = "Integer";
			npm.yon_type = YON_VCF_HEADER_INTEGER;
			npm.description = "\"Total number of alleles in called genotypes\"";
			npm.idx = this->info_fields_.size();
			this->literals_ += npm.ToVcfString(false) + "\n";
			this->info_fields_.push_back(npm);
		}

		info = this->GetInfo(group_names[i] + "_HWE_P");
		if(info == nullptr){
			YonInfo npm;
			npm.id = group_names[i] + "_HWE_P";
			npm.number = "1";
			npm.type = "Float";
			npm.yon_type = YON_VCF_HEADER_FLOAT;
			npm.description = "\"Hardy-Weinberg equilibrium P-value\"";
			npm.idx = this->info_fields_.size();
			this->literals_ += npm.ToVcfString(false) + "\n";
			this->info_fields_.push_back(npm);
		}

		info = this->GetInfo(group_names[i] + "_AC");
		if(info == nullptr){
			YonInfo npm;
			npm.id = group_names[i] + "_AC";
			npm.number = "A";
			npm.type = "Integer";
			npm.yon_type = YON_VCF_HEADER_INTEGER;
			npm.description = "\"Allele count in genotypes, for each REF and ALT allele, in the same order as listed with REF first\"";
			npm.idx = this->info_fields_.size();
			this->literals_ += npm.ToVcfString(false) + "\n";
			this->info_fields_.push_back(npm);
		}

		info = this->GetInfo(group_names[i] + "_AF");
		if(info == nullptr){
			YonInfo npm;
			npm.id = group_names[i] + "_AF";
			npm.number = "A";
			npm.type = "Float";
			npm.yon_type = YON_VCF_HEADER_FLOAT;
			npm.description = "\"Allele Frequency, for each REF and ALT allele, in the same order as listed with REF first. In the range [0, 1]\"";
			npm.idx = this->info_fields_.size();
			this->literals_ += npm.ToVcfString(false) + "\n";
			this->info_fields_.push_back(npm);
		}

		info = this->GetInfo(group_names[i] + "_AC_P");
		if(info == nullptr){
			YonInfo npm;
			npm.id = group_names[i] + "_AC_P";
			npm.number = "A";
			npm.type = "Integer";
			npm.yon_type = YON_VCF_HEADER_INTEGER;
			npm.description = "\"AC_P\"";
			npm.idx = this->info_fields_.size();
			this->literals_ += npm.ToVcfString(false) + "\n";
			this->info_fields_.push_back(npm);
		}

		info = this->GetInfo(group_names[i] + "_FS_A");
		if(info == nullptr){
			YonInfo npm;
			npm.id = group_names[i] + "_FS_A";
			npm.number = "A";
			npm.type = "Float";
			npm.yon_type = YON_VCF_HEADER_FLOAT;
			npm.description = "\"FS_A\"";
			npm.idx = this->info_fields_.size();
			this->literals_ += npm.ToVcfString(false) + "\n";
			this->info_fields_.push_back(npm);
		}

		info = this->GetInfo(group_names[i] + "_F_PIC");
		if(info == nullptr){
			YonInfo nm;
			nm.id = group_names[i] + "_F_PIC";
			nm.number = "1";
			nm.type = "Float";
			nm.yon_type = YON_VCF_HEADER_FLOAT;
			nm.description = "\"F_PIC\"";
			nm.idx = this->info_fields_.size();
			this->literals_ += nm.ToVcfString(false) + "\n";
			this->info_fields_.push_back(nm);
		}

		info = this->GetInfo(group_names[i] + "_HET");
		if(info == nullptr){
			YonInfo nm;
			nm.id = group_names[i] + "_HET";
			nm.number = "1";
			nm.type = "Float";
			nm.yon_type = YON_VCF_HEADER_FLOAT;
			nm.description = "\"Heterozygosity at this locus calculated as number of 0/1 or 1/0 genotypes divided by all non-missing genotypes\"";
			nm.idx = this->info_fields_.size();
			this->literals_ += nm.ToVcfString(false) + "\n";
			this->info_fields_.push_back(nm);
		}

		info = this->GetInfo(group_names[i] + "_MULTI_ALLELIC");
		if(info == nullptr){
			YonInfo nm;
			nm.id = group_names[i] + "_MULTI_ALLELIC";
			nm.number = "1";
			nm.type = "Flag";
			nm.yon_type = YON_VCF_HEADER_FLAG;
			nm.description = "\"Flag indicating if a site is multi-allelic (>1 ALT alleles)\"";
			nm.idx = this->info_fields_.size();
			this->literals_ += nm.ToVcfString(false) + "\n";
			this->info_fields_.push_back(nm);
		}

		this->BuildMaps();
		this->BuildReverseMaps();
	}
}

std::ostream& VariantHeader::PrintVcfHeader(std::ostream& stream) const{
	stream << this->literals_;
	stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if(this->samples_.size()){
		stream << "\tFORMAT\t";
		stream << this->samples_[0];
		for(size_t i = 1; i < this->samples_.size(); ++i)
			stream << "\t" + this->samples_[i];
	}
	stream << "\n";
	return(stream);
}


std::string VariantHeader::ToString(const bool is_bcf) const{
	std::string string = "##fileformat=VCFv4.1\n";
	uint32_t idx = 0;
	for(uint32_t i = 0; i < this->contigs_.size(); ++i)       string += this->contigs_[i].ToVcfString(is_bcf) + "\n";
	for(uint32_t i = 0; i < this->structured_extra_fields_.size(); ++i) string += this->structured_extra_fields_[i].ToVcfString() + "\n";
	for(uint32_t i = 0; i < this->filter_fields_.size(); ++i) string += this->filter_fields_[i].ToVcfString(idx++) + "\n";
	for(uint32_t i = 0; i < this->info_fields_.size(); ++i)   string += this->info_fields_[i].ToVcfString(idx++) + "\n";
	for(uint32_t i = 0; i < this->format_fields_.size(); ++i) string += this->format_fields_[i].ToVcfString(idx++) + "\n";
	for(uint32_t i = 0; i < this->extra_fields_.size(); ++i)  string += this->extra_fields_[i].ToVcfString() + "\n";
	return(string);
}

std::ostream& operator<<(std::ostream& stream, const VariantHeader& header){
	utility::SerializeString(header.fileformat_string_, stream);
	utility::SerializeString(header.literals_, stream);

	size_t l_helper = header.samples_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(uint32_t i = 0; i < header.samples_.size(); ++i) utility::SerializeString(header.samples_[i], stream);

	l_helper = header.contigs_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(uint32_t i = 0; i < header.contigs_.size(); ++i) stream << header.contigs_[i];

	l_helper = header.info_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(uint32_t i = 0; i < header.info_fields_.size(); ++i) stream << header.info_fields_[i];

	l_helper = header.format_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(uint32_t i = 0; i < header.format_fields_.size(); ++i) stream << header.format_fields_[i];

	l_helper = header.filter_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(uint32_t i = 0; i < header.filter_fields_.size(); ++i) stream << header.filter_fields_[i];

	l_helper = header.structured_extra_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(uint32_t i = 0; i < header.structured_extra_fields_.size(); ++i) stream << header.structured_extra_fields_[i];

	l_helper = header.extra_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(uint32_t i = 0; i < header.extra_fields_.size(); ++i) stream << header.extra_fields_[i];

	return(stream);
}

std::istream& operator>>(std::istream& stream, VariantHeader& header){
	utility::DeserializeString(header.fileformat_string_, stream);
	utility::DeserializeString(header.literals_, stream);

	size_t l_helper;
	utility::DeserializePrimitive(l_helper, stream);
	header.samples_.resize(l_helper);
	for(uint32_t i = 0; i < header.samples_.size(); ++i) utility::DeserializeString(header.samples_[i], stream);

	utility::DeserializePrimitive(l_helper, stream);
	header.contigs_.resize(l_helper);
	for(uint32_t i = 0; i < header.contigs_.size(); ++i) stream >> header.contigs_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.info_fields_.resize(l_helper);
	for(uint32_t i = 0; i < header.info_fields_.size(); ++i) stream >> header.info_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.format_fields_.resize(l_helper);
	for(uint32_t i = 0; i < header.format_fields_.size(); ++i) stream >> header.format_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.filter_fields_.resize(l_helper);
	for(uint32_t i = 0; i < header.filter_fields_.size(); ++i) stream >> header.filter_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.structured_extra_fields_.resize(l_helper);
	for(uint32_t i = 0; i < header.structured_extra_fields_.size(); ++i) stream >> header.structured_extra_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.extra_fields_.resize(l_helper);
	for(uint32_t i = 0; i < header.extra_fields_.size(); ++i) stream >> header.extra_fields_[i];

	header.BuildMaps();
	header.BuildReverseMaps();

	return stream;
}

io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VariantHeader& header){
	io::SerializeString(header.fileformat_string_, buffer);
	io::SerializeString(header.literals_, buffer);

	uint32_t l_helper = header.samples_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(uint32_t i = 0; i < header.samples_.size(); ++i) io::SerializeString(header.samples_[i], buffer);

	l_helper = header.contigs_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(uint32_t i = 0; i < header.contigs_.size(); ++i) buffer << header.contigs_[i];

	l_helper = header.info_fields_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(uint32_t i = 0; i < header.info_fields_.size(); ++i) buffer << header.info_fields_[i];

	l_helper = header.format_fields_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(uint32_t i = 0; i < header.format_fields_.size(); ++i) buffer << header.format_fields_[i];

	l_helper = header.filter_fields_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(uint32_t i = 0; i < header.filter_fields_.size(); ++i) buffer << header.filter_fields_[i];

	l_helper = header.structured_extra_fields_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(uint32_t i = 0; i < header.structured_extra_fields_.size(); ++i) buffer << header.structured_extra_fields_[i];

	l_helper = header.extra_fields_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(uint32_t i = 0; i < header.extra_fields_.size(); ++i) buffer << header.extra_fields_[i];

	return(buffer);
}

io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VariantHeader& header){
	io::DeserializeString(header.fileformat_string_, buffer);
	io::DeserializeString(header.literals_, buffer);

	uint32_t l_helper;
	io::DeserializePrimitive(l_helper, buffer);
	header.samples_.resize(l_helper);
	for(uint32_t i = 0; i < header.samples_.size(); ++i) io::DeserializeString(header.samples_[i], buffer);

	io::DeserializePrimitive(l_helper, buffer);
	header.contigs_.resize(l_helper);
	for(uint32_t i = 0; i < header.contigs_.size(); ++i)       buffer >> header.contigs_[i];

	io::DeserializePrimitive(l_helper, buffer);
	header.info_fields_.resize(l_helper);
	for(uint32_t i = 0; i < header.info_fields_.size(); ++i)   buffer >> header.info_fields_[i];

	io::DeserializePrimitive(l_helper, buffer);
	header.format_fields_.resize(l_helper);
	for(uint32_t i = 0; i < header.format_fields_.size(); ++i) buffer >> header.format_fields_[i];

	io::DeserializePrimitive(l_helper, buffer);
	header.filter_fields_.resize(l_helper);
	for(uint32_t i = 0; i < header.filter_fields_.size(); ++i) buffer >> header.filter_fields_[i];

	io::DeserializePrimitive(l_helper, buffer);
	header.structured_extra_fields_.resize(l_helper);
	for(uint32_t i = 0; i < header.structured_extra_fields_.size(); ++i) buffer >> header.structured_extra_fields_[i];

	io::DeserializePrimitive(l_helper, buffer);
	header.extra_fields_.resize(l_helper);
	for(uint32_t i = 0; i < header.extra_fields_.size(); ++i) buffer >> header.extra_fields_[i];

	header.BuildMaps();
	header.BuildReverseMaps();

	return buffer;
}

}
