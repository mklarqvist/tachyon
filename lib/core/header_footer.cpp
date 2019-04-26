#include "header_footer.h"
#include "utility.h"

namespace tachyon {

yon_vnt_hdr_t::yon_vnt_hdr_t(const yon_vnt_hdr_t& other) :
	fileformat_string_(other.fileformat_string_),
	literals_(other.literals_),
	samples_(other.samples_),
	contigs_(other.contigs_),
	info_fields_(other.info_fields_),
	format_fields_(other.format_fields_),
	filter_fields_(other.filter_fields_),
	structured_extra_fields_(other.structured_extra_fields_),
	extra_fields_(other.extra_fields_)
{
	this->BuildMaps();
	this->BuildReverseMaps();
}

void yon_vnt_hdr_t::AddSample(const std::string& sample_name) {
	if (this->samples_.size() == 0) {
		this->samples_.push_back(sample_name);
		this->samples_map_[sample_name] = 0;
		return;
	}

	if (this->samples_map_.find(sample_name) == this->samples_map_.end()) {
		this->samples_map_[sample_name] = this->samples_.size();
		this->samples_.push_back(sample_name);
	} else {
		std::cerr << utility::timestamp("ERROR", "HEADER") << "Illegal: duplicated sample name: " << sample_name << std::endl;
	}
}

const YonContig* yon_vnt_hdr_t::GetContig(const std::string& name) const {
	map_type::const_iterator it = this->contigs_map_.find(name);
	if (it != this->contigs_map_.end()) return(&this->contigs_[it->second]);
	return(nullptr);
}

const YonContig* yon_vnt_hdr_t::GetContig(const int& idx) const {
	map_reverse_type::const_iterator it = this->contigs_reverse_map_.find(idx);
	if (it != this->contigs_reverse_map_.end()) return(&this->contigs_[it->second]);
	return(nullptr);
}

const YonInfo* yon_vnt_hdr_t::GetInfo(const std::string& name) const {
	map_type::const_iterator it = this->info_fields_map_.find(name);
	if (it != this->info_fields_map_.end()) return(&this->info_fields_[it->second]);
	return(nullptr);
}

const YonInfo* yon_vnt_hdr_t::GetInfo(const int& idx) const {
	map_reverse_type::const_iterator it = this->info_fields_reverse_map_.find(idx);
	if (it != this->info_fields_reverse_map_.end()) return(&this->info_fields_[it->second]);
	return(nullptr);
}

const YonFormat* yon_vnt_hdr_t::GetFormat(const std::string& name) const {
	map_type::const_iterator it = this->format_fields_map_.find(name);
	if (it != this->format_fields_map_.end()) return(&this->format_fields_[it->second]);
	return(nullptr);
}

const YonFormat* yon_vnt_hdr_t::GetFormat(const int& idx) const {
	map_reverse_type::const_iterator it = this->format_fields_reverse_map_.find(idx);
	if (it != this->format_fields_reverse_map_.end()) return(&this->format_fields_[it->second]);
	return(nullptr);
}

const VcfFilter* yon_vnt_hdr_t::GetFilter(const std::string& name) const {
	map_type::const_iterator it = this->filter_fields_map_.find(name);
	if (it != this->filter_fields_map_.end()) return(&this->filter_fields_[it->second]);
	return(nullptr);
}

const VcfFilter* yon_vnt_hdr_t::GetFilter(const int& idx) const {
	map_reverse_type::const_iterator it = this->filter_fields_reverse_map_.find(idx);
	if (it != this->filter_fields_reverse_map_.end()) return(&this->filter_fields_[it->second]);
	return(nullptr);
}

const std::string* yon_vnt_hdr_t::GetSample(const std::string& name) const {
	map_type::const_iterator it = this->samples_map_.find(name);
	if (it != this->samples_map_.end()) return(&this->samples_[it->second]);
	return(nullptr);
}

int32_t yon_vnt_hdr_t::GetSampleId(const std::string& name) const {
	map_type::const_iterator it = this->samples_map_.find(name);
	if (it != this->samples_map_.end()) return(it->second);
	return(-1);
}

bool yon_vnt_hdr_t::BuildReverseMaps(void) {
	this->contigs_reverse_map_.clear();
	this->info_fields_reverse_map_.clear();
	this->format_fields_reverse_map_.clear();
	this->filter_fields_reverse_map_.clear();

	for (uint32_t i = 0; i < this->contigs_.size(); ++i)       this->contigs_reverse_map_[this->contigs_[i].idx] = i;
	for (uint32_t i = 0; i < this->info_fields_.size(); ++i)   this->info_fields_reverse_map_[this->info_fields_[i].idx] = i;
	for (uint32_t i = 0; i < this->format_fields_.size(); ++i) this->format_fields_reverse_map_[this->format_fields_[i].idx] = i;
	for (uint32_t i = 0; i < this->filter_fields_.size(); ++i) this->filter_fields_reverse_map_[this->filter_fields_[i].idx] = i;

	return true;
}

bool yon_vnt_hdr_t::BuildMaps(void) {
	this->info_fields_map_.clear();
	this->format_fields_map_.clear();
	this->filter_fields_map_.clear();
	this->contigs_map_.clear();

	for (uint32_t i = 0; i < this->contigs_.size(); ++i)       this->contigs_map_[this->contigs_[i].name] = i;
	for (uint32_t i = 0; i < this->info_fields_.size(); ++i)   this->info_fields_map_[this->info_fields_[i].id] = i;
	for (uint32_t i = 0; i < this->format_fields_.size(); ++i) this->format_fields_map_[this->format_fields_[i].id] = i;
	for (uint32_t i = 0; i < this->filter_fields_.size(); ++i) this->filter_fields_map_[this->filter_fields_[i].id] = i;
	for (uint32_t i = 0; i < this->samples_.size(); ++i)       this->samples_map_[this->samples_[i]] = i;

	return true;
}

bool yon_vnt_hdr_t::RecodeIndices(void) {
	for (uint32_t i = 0; i < this->contigs_.size(); ++i)       this->contigs_[i].idx = i;
	for (uint32_t i = 0; i < this->info_fields_.size(); ++i)   this->info_fields_[i].idx = i;
	for (uint32_t i = 0; i < this->format_fields_.size(); ++i) this->format_fields_[i].idx = i;
	for (uint32_t i = 0; i < this->filter_fields_.size(); ++i) this->filter_fields_[i].idx = i;

	if (this->BuildMaps() == false) return false;
	if (this->BuildReverseMaps() == false) return false;
	return true;
}



void yon_vnt_hdr_t::AddGenotypeAnnotationFields(void) {
	//"NM","NPM","AN","HWE_P","AC","AF","AC_P","FS_A","F_PIC","HET","MULTI_ALLELIC"

	const YonInfo* info = this->GetInfo("NM");
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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
	if (info == nullptr) {
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


void yon_vnt_hdr_t::AddGenotypeAnnotationFields(const std::vector<std::string>& group_names) {
	for (int i = 0; i < group_names.size(); ++i) {
		const YonInfo* info = this->GetInfo(group_names[i] + "_NM");
		if (info == nullptr) {
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
		if (info == nullptr) {
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
		if (info == nullptr) {
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
		if (info == nullptr) {
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
		if (info == nullptr) {
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
		if (info == nullptr) {
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
		if (info == nullptr) {
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
		if (info == nullptr) {
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
		if (info == nullptr) {
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
		if (info == nullptr) {
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
		if (info == nullptr) {
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

std::ostream& yon_vnt_hdr_t::PrintVcfHeader(std::ostream& stream) const {
	stream << this->literals_;
	stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
	if (this->samples_.size()) {
		stream << "\tFORMAT\t";
		stream << this->samples_[0];
		for (size_t i = 1; i < this->samples_.size(); ++i)
			stream << "\t" + this->samples_[i];
	}
	stream << "\n";
	return(stream);
}


std::string yon_vnt_hdr_t::ToString(const bool is_bcf) const {
	std::string string = "##fileformat=VCFv4.1\n";
	uint32_t idx = 0;
	for (uint32_t i = 0; i < this->contigs_.size(); ++i)       string += this->contigs_[i].ToVcfString(is_bcf) + "\n";
	for (uint32_t i = 0; i < this->structured_extra_fields_.size(); ++i) string += this->structured_extra_fields_[i].ToVcfString() + "\n";
	for (uint32_t i = 0; i < this->filter_fields_.size(); ++i) string += this->filter_fields_[i].ToVcfString(idx++) + "\n";
	for (uint32_t i = 0; i < this->info_fields_.size(); ++i)   string += this->info_fields_[i].ToVcfString(idx++) + "\n";
	for (uint32_t i = 0; i < this->format_fields_.size(); ++i) string += this->format_fields_[i].ToVcfString(idx++) + "\n";
	for (uint32_t i = 0; i < this->extra_fields_.size(); ++i)  string += this->extra_fields_[i].ToVcfString() + "\n";
	return(string);
}

std::ostream& operator<<(std::ostream& stream, const yon_vnt_hdr_t& header) {
	utility::SerializeString(header.fileformat_string_, stream);
	utility::SerializeString(header.literals_, stream);

	size_t l_helper = header.samples_.size();
	utility::SerializePrimitive(l_helper, stream);
	for (uint32_t i = 0; i < header.samples_.size(); ++i) utility::SerializeString(header.samples_[i], stream);

	l_helper = header.contigs_.size();
	utility::SerializePrimitive(l_helper, stream);
	for (uint32_t i = 0; i < header.contigs_.size(); ++i) stream << header.contigs_[i];

	l_helper = header.info_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for (uint32_t i = 0; i < header.info_fields_.size(); ++i) stream << header.info_fields_[i];

	l_helper = header.format_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for (uint32_t i = 0; i < header.format_fields_.size(); ++i) stream << header.format_fields_[i];

	l_helper = header.filter_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for (uint32_t i = 0; i < header.filter_fields_.size(); ++i) stream << header.filter_fields_[i];

	l_helper = header.structured_extra_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for (uint32_t i = 0; i < header.structured_extra_fields_.size(); ++i) stream << header.structured_extra_fields_[i];

	l_helper = header.extra_fields_.size();
	utility::SerializePrimitive(l_helper, stream);
	for (uint32_t i = 0; i < header.extra_fields_.size(); ++i) stream << header.extra_fields_[i];

	return(stream);
}

std::istream& operator>>(std::istream& stream, yon_vnt_hdr_t& header) {
	utility::DeserializeString(header.fileformat_string_, stream);
	utility::DeserializeString(header.literals_, stream);

	size_t l_helper;
	utility::DeserializePrimitive(l_helper, stream);
	header.samples_.resize(l_helper);
	for (uint32_t i = 0; i < header.samples_.size(); ++i) utility::DeserializeString(header.samples_[i], stream);

	utility::DeserializePrimitive(l_helper, stream);
	header.contigs_.resize(l_helper);
	for (uint32_t i = 0; i < header.contigs_.size(); ++i) stream >> header.contigs_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.info_fields_.resize(l_helper);
	for (uint32_t i = 0; i < header.info_fields_.size(); ++i) stream >> header.info_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.format_fields_.resize(l_helper);
	for (uint32_t i = 0; i < header.format_fields_.size(); ++i) stream >> header.format_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.filter_fields_.resize(l_helper);
	for (uint32_t i = 0; i < header.filter_fields_.size(); ++i) stream >> header.filter_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.structured_extra_fields_.resize(l_helper);
	for (uint32_t i = 0; i < header.structured_extra_fields_.size(); ++i) stream >> header.structured_extra_fields_[i];

	utility::DeserializePrimitive(l_helper, stream);
	header.extra_fields_.resize(l_helper);
	for (uint32_t i = 0; i < header.extra_fields_.size(); ++i) stream >> header.extra_fields_[i];

	header.BuildMaps();
	header.BuildReverseMaps();

	return stream;
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const yon_vnt_hdr_t& header) {
	SerializeString(header.fileformat_string_, buffer);
	SerializeString(header.literals_, buffer);

	uint32_t l_helper = header.samples_.size();
	SerializePrimitive(l_helper, buffer);
	for (uint32_t i = 0; i < header.samples_.size(); ++i) SerializeString(header.samples_[i], buffer);

	l_helper = header.contigs_.size();
	SerializePrimitive(l_helper, buffer);
	for (uint32_t i = 0; i < header.contigs_.size(); ++i) buffer << header.contigs_[i];

	l_helper = header.info_fields_.size();
	SerializePrimitive(l_helper, buffer);
	for (uint32_t i = 0; i < header.info_fields_.size(); ++i) buffer << header.info_fields_[i];

	l_helper = header.format_fields_.size();
	SerializePrimitive(l_helper, buffer);
	for (uint32_t i = 0; i < header.format_fields_.size(); ++i) buffer << header.format_fields_[i];

	l_helper = header.filter_fields_.size();
	SerializePrimitive(l_helper, buffer);
	for (uint32_t i = 0; i < header.filter_fields_.size(); ++i) buffer << header.filter_fields_[i];

	l_helper = header.structured_extra_fields_.size();
	SerializePrimitive(l_helper, buffer);
	for (uint32_t i = 0; i < header.structured_extra_fields_.size(); ++i) buffer << header.structured_extra_fields_[i];

	l_helper = header.extra_fields_.size();
	SerializePrimitive(l_helper, buffer);
	for (uint32_t i = 0; i < header.extra_fields_.size(); ++i) buffer << header.extra_fields_[i];

	return(buffer);
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, yon_vnt_hdr_t& header) {
	DeserializeString(header.fileformat_string_, buffer);
	DeserializeString(header.literals_, buffer);

	uint32_t l_helper;
	DeserializePrimitive(l_helper, buffer);
	header.samples_.resize(l_helper);
	for (uint32_t i = 0; i < header.samples_.size(); ++i) DeserializeString(header.samples_[i], buffer);

	DeserializePrimitive(l_helper, buffer);
	header.contigs_.resize(l_helper);
	for (uint32_t i = 0; i < header.contigs_.size(); ++i)       buffer >> header.contigs_[i];

	DeserializePrimitive(l_helper, buffer);
	header.info_fields_.resize(l_helper);
	for (uint32_t i = 0; i < header.info_fields_.size(); ++i)   buffer >> header.info_fields_[i];

	DeserializePrimitive(l_helper, buffer);
	header.format_fields_.resize(l_helper);
	for (uint32_t i = 0; i < header.format_fields_.size(); ++i) buffer >> header.format_fields_[i];

	DeserializePrimitive(l_helper, buffer);
	header.filter_fields_.resize(l_helper);
	for (uint32_t i = 0; i < header.filter_fields_.size(); ++i) buffer >> header.filter_fields_[i];

	DeserializePrimitive(l_helper, buffer);
	header.structured_extra_fields_.resize(l_helper);
	for (uint32_t i = 0; i < header.structured_extra_fields_.size(); ++i) buffer >> header.structured_extra_fields_[i];

	DeserializePrimitive(l_helper, buffer);
	header.extra_fields_.resize(l_helper);
	for (uint32_t i = 0; i < header.extra_fields_.size(); ++i) buffer >> header.extra_fields_[i];

	header.BuildMaps();
	header.BuildReverseMaps();

	return buffer;
}

yon_ftr_t::yon_ftr_t() :
	offset_end_of_data(0),
	n_blocks(0),
	n_variants(0),
	controller(0)
{
	utility::HexToBytes(TACHYON_FILE_EOF, &this->EOF_marker[0]);
}

yon_ftr_t::yon_ftr_t(const char* const data) :
	offset_end_of_data(*reinterpret_cast<const uint64_t* const>(data)),
	n_blocks(*reinterpret_cast<const uint64_t* const>(&data[sizeof(uint64_t)])),
	n_variants(*reinterpret_cast<const uint64_t* const>(&data[sizeof(uint64_t)*2])),
	controller(*reinterpret_cast<const uint16_t* const>(&data[sizeof(uint64_t)*3]))
{
	memcpy(&this->EOF_marker[0], &data[sizeof(uint64_t)*3+sizeof(uint16_t)], TACHYON_FILE_EOF_LENGTH);
}

yon_ftr_t::yon_ftr_t(const self_type& other) :
	offset_end_of_data(other.offset_end_of_data),
	n_blocks(other.n_blocks),
	n_variants(other.n_variants),
	controller(other.controller)
{
	memcpy(&this->EOF_marker[0], &other.EOF_marker[0], TACHYON_FILE_EOF_LENGTH);
}

bool yon_ftr_t::Validate(void) const {
	if (this->offset_end_of_data == 0) return false;
	if (this->n_blocks  == 0)          return false;
	if (this->n_variants == 0)         return false;

	// Check EOF marker
	uint8_t reference[TACHYON_FILE_EOF_LENGTH];
	utility::HexToBytes(TACHYON_FILE_EOF, &reference[0]);

	if (strncmp(reinterpret_cast<const char* const>(&this->EOF_marker[0]), reinterpret_cast<const char* const>(&reference[0]), TACHYON_FILE_EOF_LENGTH) != 0) return false;
	return true;
}

std::ostream& operator<<(std::ostream& stream, const yon_ftr_t& yon_ftr_t) {
	utility::SerializePrimitive(yon_ftr_t.offset_end_of_data, stream);
	utility::SerializePrimitive(yon_ftr_t.n_blocks,   stream);
	utility::SerializePrimitive(yon_ftr_t.n_variants, stream);
	utility::SerializePrimitive(yon_ftr_t.controller, stream);
	stream.write(reinterpret_cast<const char*>(&yon_ftr_t.EOF_marker[0]), TACHYON_FILE_EOF_LENGTH);
	return(stream);
}

std::istream& operator>>(std::istream& stream, yon_ftr_t& yon_ftr_t) {
	utility::DeserializePrimitive(yon_ftr_t.offset_end_of_data, stream);
	utility::DeserializePrimitive(yon_ftr_t.n_blocks,   stream);
	utility::DeserializePrimitive(yon_ftr_t.n_variants, stream);
	utility::DeserializePrimitive(yon_ftr_t.controller, stream);
	stream.read(reinterpret_cast<char*>(&yon_ftr_t.EOF_marker[0]), TACHYON_FILE_EOF_LENGTH);
	return(stream);
}

}
