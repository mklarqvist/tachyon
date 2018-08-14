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
		nm.description = "NM";
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
		npm.description = "NPM";
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
		npm.description = "AN";
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
		npm.description = "HWE_P";
		npm.idx = this->info_fields_.size();
		this->literals_ += npm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(npm);
	}

	info = this->GetInfo("AC");
	if(info == nullptr){
		YonInfo npm;
		npm.id = "AC";
		npm.number = ".";
		npm.type = "Integer";
		npm.yon_type = YON_VCF_HEADER_INTEGER;
		npm.description = "AC";
		npm.idx = this->info_fields_.size();
		this->literals_ += npm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(npm);
	}

	info = this->GetInfo("AF");
	if(info == nullptr){
		YonInfo npm;
		npm.id = "AF";
		npm.number = ".";
		npm.type = "Float";
		npm.yon_type = YON_VCF_HEADER_FLOAT;
		npm.description = "AF";
		npm.idx = this->info_fields_.size();
		this->literals_ += npm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(npm);
	}

	info = this->GetInfo("AC_P");
	if(info == nullptr){
		YonInfo npm;
		npm.id = "AC_P";
		npm.number = ".";
		npm.type = "Integer";
		npm.yon_type = YON_VCF_HEADER_INTEGER;
		npm.description = "AC_P";
		npm.idx = this->info_fields_.size();
		this->literals_ += npm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(npm);
	}

	info = this->GetInfo("FS_A");
	if(info == nullptr){
		YonInfo npm;
		npm.id = "FS_A";
		npm.number = ".";
		npm.type = "Float";
		npm.yon_type = YON_VCF_HEADER_FLOAT;
		npm.description = "FS_A";
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
		nm.description = "F_PIC";
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
		nm.description = "HET";
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
		nm.description = "MULTI_ALLELIC";
		nm.idx = this->info_fields_.size();
		this->literals_ += nm.ToVcfString(false) + "\n";
		this->info_fields_.push_back(nm);
	}

	this->BuildMaps();
	this->BuildReverseMaps();
}

std::ostream& operator<<(std::ostream& stream, const VariantHeader& header){
	utility::SerializeString(header.fileformat_string_, stream);
	utility::SerializeString(header.literals_, stream);

	size_t l_helper = header.samples_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(U32 i = 0; i < header.samples_.size(); ++i) utility::SerializeString(header.samples_[i], stream);

	l_helper = header.contigs_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(U32 i = 0; i < header.contigs_.size(); ++i) stream << header.contigs_[i];

	l_helper = header.info_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(U32 i = 0; i < header.info_fields_.size(); ++i) stream << header.info_fields_[i];

	l_helper = header.format_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(U32 i = 0; i < header.format_fields_.size(); ++i) stream << header.format_fields_[i];

	l_helper = header.filter_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(U32 i = 0; i < header.filter_fields_.size(); ++i) stream << header.filter_fields_[i];

	l_helper = header.structured_extra_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(U32 i = 0; i < header.structured_extra_fields_.size(); ++i) stream << header.structured_extra_fields_[i];

	l_helper = header.extra_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for(U32 i = 0; i < header.extra_fields_.size(); ++i) stream << header.extra_fields_[i];

	return(stream);
}

std::istream& operator>>(std::istream& stream, VariantHeader& header){
	utility::DeserializeString(header.fileformat_string_, stream);
	utility::DeserializeString(header.literals_, stream);

	size_t l_helper;
	utility::DeserializePrimitive(l_helper, stream);
	header.samples_.resize(l_helper);
	for(U32 i = 0; i < header.samples_.size(); ++i) utility::DeserializeString(header.samples_[i], stream);

	utility::DeserializePrimitive(l_helper, stream);
	header.contigs_.resize(l_helper);
	for(U32 i = 0; i < header.contigs_.size(); ++i) stream >> header.contigs_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.info_fields_.resize(l_helper);
	for(U32 i = 0; i < header.info_fields_.size(); ++i) stream >> header.info_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.format_fields_.resize(l_helper);
	for(U32 i = 0; i < header.format_fields_.size(); ++i) stream >> header.format_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.filter_fields_.resize(l_helper);
	for(U32 i = 0; i < header.filter_fields_.size(); ++i) stream >> header.filter_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.structured_extra_fields_.resize(l_helper);
	for(U32 i = 0; i < header.structured_extra_fields_.size(); ++i) stream >> header.structured_extra_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.extra_fields_.resize(l_helper);
	for(U32 i = 0; i < header.extra_fields_.size(); ++i) stream >> header.extra_fields_[i];

	header.BuildMaps();
	header.BuildReverseMaps();

	return stream;
}

io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VariantHeader& header){
	io::SerializeString(header.fileformat_string_, buffer);
	io::SerializeString(header.literals_, buffer);

	uint32_t l_helper = header.samples_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(U32 i = 0; i < header.samples_.size(); ++i) io::SerializeString(header.samples_[i], buffer);

	l_helper = header.contigs_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(U32 i = 0; i < header.contigs_.size(); ++i) buffer << header.contigs_[i];

	l_helper = header.info_fields_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(U32 i = 0; i < header.info_fields_.size(); ++i) buffer << header.info_fields_[i];

	l_helper = header.format_fields_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(U32 i = 0; i < header.format_fields_.size(); ++i) buffer << header.format_fields_[i];

	l_helper = header.filter_fields_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(U32 i = 0; i < header.filter_fields_.size(); ++i) buffer << header.filter_fields_[i];

	l_helper = header.structured_extra_fields_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(U32 i = 0; i < header.structured_extra_fields_.size(); ++i) buffer << header.structured_extra_fields_[i];

	l_helper = header.extra_fields_.size();
	io::SerializePrimitive(l_helper, buffer);
	for(U32 i = 0; i < header.extra_fields_.size(); ++i) buffer << header.extra_fields_[i];

	return(buffer);
}

io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VariantHeader& header){
	io::DeserializeString(header.fileformat_string_, buffer);
	io::DeserializeString(header.literals_, buffer);

	uint32_t l_helper;
	io::DeserializePrimitive(l_helper, buffer);
	header.samples_.resize(l_helper);
	for(U32 i = 0; i < header.samples_.size(); ++i) io::DeserializeString(header.samples_[i], buffer);

	io::DeserializePrimitive(l_helper, buffer);
	header.contigs_.resize(l_helper);
	for(U32 i = 0; i < header.contigs_.size(); ++i)       buffer >> header.contigs_[i];

	io::DeserializePrimitive(l_helper, buffer);
	header.info_fields_.resize(l_helper);
	for(U32 i = 0; i < header.info_fields_.size(); ++i)   buffer >> header.info_fields_[i];

	io::DeserializePrimitive(l_helper, buffer);
	header.format_fields_.resize(l_helper);
	for(U32 i = 0; i < header.format_fields_.size(); ++i) buffer >> header.format_fields_[i];

	io::DeserializePrimitive(l_helper, buffer);
	header.filter_fields_.resize(l_helper);
	for(U32 i = 0; i < header.filter_fields_.size(); ++i) buffer >> header.filter_fields_[i];

	io::DeserializePrimitive(l_helper, buffer);
	header.structured_extra_fields_.resize(l_helper);
	for(U32 i = 0; i < header.structured_extra_fields_.size(); ++i) buffer >> header.structured_extra_fields_[i];

	io::DeserializePrimitive(l_helper, buffer);
	header.extra_fields_.resize(l_helper);
	for(U32 i = 0; i < header.extra_fields_.size(); ++i) buffer >> header.extra_fields_[i];

	header.BuildMaps();
	header.BuildReverseMaps();

	return buffer;
}

}
