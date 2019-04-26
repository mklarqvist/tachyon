#include "tachyon.h"
#include "variant_record.h"

namespace tachyon{

/****************************
*  Allele description
****************************/
yon_allele::yon_allele() : l_allele(0), allele(nullptr) {}
yon_allele::yon_allele(const std::string& string) : l_allele(string.size()), allele(new char[l_allele]) { memcpy(allele, string.data(), l_allele); }
yon_allele& yon_allele::operator=(const std::string& value) { delete [] allele; l_allele = value.size(); allele = new char[l_allele]; memcpy(allele, value.data(), l_allele); return(*this); }
yon_allele& yon_allele::operator=(const char* value) { delete [] allele; l_allele = strlen(value); allele = new char[l_allele]; memcpy(allele, value, l_allele); return(*this); }
yon_allele& yon_allele::operator=(const char& value) { delete [] allele; l_allele = 1; allele = new char[l_allele]; allele[0] = value; return(*this); }

yon_allele::yon_allele(const yon_allele& other) : l_allele(other.l_allele), allele(new char[l_allele]) { memcpy(allele, other.allele, l_allele); }
yon_allele::yon_allele(yon_allele&& other) noexcept : l_allele(other.l_allele), allele(nullptr) { std::swap(allele, other.allele); other.l_allele = 0; }

yon_allele& yon_allele::operator=(const yon_allele& other) {
	delete [] allele;
	l_allele = other.l_allele;
	allele = new char[l_allele];
	memcpy(allele, other.allele, l_allele);
	return(*this);
}

yon_allele& yon_allele::operator=(yon_allele&& other) noexcept {
	if (this == &other) {
		// precautions against self-moves
		return *this;
	}
	l_allele = other.l_allele;
	std::swap(allele, other.allele);
	other.l_allele = 0;
	return(*this);
}

yon_allele::~yon_allele() { delete [] allele; }

yon_allele& yon_allele::ParseFromBuffer(const char* const in) {
	delete [] this->allele;
	this->l_allele = *reinterpret_cast<const uint16_t* const>(in);
	this->allele = new char[this->l_allele];
	memcpy(this->allele, &in[sizeof(uint16_t)], this->l_allele);
	return(*this);
}

/****************************
*  Core record
****************************/
yon1_vnt_t::yon1_vnt_t() :
	is_dirty(false), is_loaded_gt(false), n_base_ploidy(0), n_alleles(0),
	n_flt(0), n_fmt(0), n_info(0), info_pid(-1), flt_pid(-1), fmt_pid(-1),
	m_fmt(0), m_info(0), m_allele(0),
	qual(NAN), rid(0), pos(0),
	alleles(nullptr), gt(nullptr), gt_sum(nullptr), gt_sum_occ(nullptr),
	info(nullptr), fmt(nullptr)
{}

yon1_vnt_t::yon1_vnt_t(const yon1_vnt_t& other) :
	is_dirty(other.is_dirty), is_loaded_gt(other.is_loaded_gt), controller(other.controller),
	n_base_ploidy(other.n_base_ploidy), n_alleles(other.n_alleles),
	n_flt(other.n_flt), n_fmt(other.n_fmt), n_info(other.n_info),
	info_pid(other.info_pid), flt_pid(other.flt_pid), fmt_pid(other.fmt_pid),
	m_fmt(other.m_fmt), m_info(other.m_info), m_allele(other.m_allele),
	qual(other.qual), rid(other.rid), pos(other.pos), name(other.name),
	info_hdr(other.info_hdr), fmt_hdr(other.fmt_hdr), flt_hdr(other.flt_hdr),
	info_map(other.info_map), fmt_map(other.fmt_map), flt_map(other.flt_map),
	alleles(new yon_allele[n_alleles]),
	gt(nullptr), gt_sum(nullptr), gt_sum_occ(nullptr),
	info(nullptr), fmt(nullptr)
{
	if (n_info) info = new PrimitiveContainerInterface*[m_info];
	if (n_fmt)  fmt  = new PrimitiveGroupContainerInterface*[m_fmt];
	for (int i = 0; i < n_alleles; ++i) alleles[i] = other.alleles[i];
	for (int i = 0; i < n_info; ++i) info[i] = other.info[i]->Clone();
	for (int i = 0; i < n_fmt; ++i)  fmt[i]  = other.fmt[i]->Clone();
	if (other.gt != nullptr) gt = new yon_gt(*other.gt);
	if (other.gt_sum != nullptr) gt_sum = new yon_gt_summary(*other.gt_sum);
	if (other.gt_sum_occ != nullptr) {
		gt_sum_occ = new yon_gt_summary[gt->n_o];
		for (int i = 0; i < gt->n_o; ++i) gt_sum_occ[i] = other.gt_sum_occ[i];
	}
}

yon1_vnt_t& yon1_vnt_t::operator=(const yon1_vnt_t& other) {
	// Delete existing data without consideration.
	delete [] alleles;
	delete gt; delete gt_sum;
	delete [] gt_sum_occ;
	for (int i = 0; i < n_info; ++i) delete info[i];
	for (int i = 0; i < n_fmt;  ++i) delete fmt[i];
	delete [] info; delete [] fmt;

	// Copy data.
	is_dirty = other.is_dirty; is_loaded_gt = other.is_loaded_gt; controller = other.controller;
	n_base_ploidy = other.n_base_ploidy; n_alleles = other.n_alleles;
	n_flt = other.n_flt; n_fmt = other.n_fmt; n_info = other.n_info;
	info_pid = other.info_pid; flt_pid = other.flt_pid; fmt_pid = other.fmt_pid;
	m_fmt = other.m_fmt; m_info = other.m_info; m_allele = other.m_allele;
	qual = other.qual; rid = other.rid; pos = other.pos; name = other.name;
	info_hdr = other.info_hdr; fmt_hdr = other.fmt_hdr; flt_hdr = other.flt_hdr;
	info_map = other.info_map; fmt_map = other.fmt_map; flt_map = other.flt_map;
	alleles = new yon_allele[n_alleles];
	gt = nullptr; gt_sum = nullptr; gt_sum_occ = nullptr,
	info = nullptr; fmt = nullptr;

	if (n_info) info = new PrimitiveContainerInterface*[m_info];
	if (n_fmt)  fmt  = new PrimitiveGroupContainerInterface*[m_fmt];
	for (int i = 0; i < n_alleles; ++i) alleles[i] = other.alleles[i];
	for (int i = 0; i < n_info; ++i) info[i] = other.info[i]->Clone();
	for (int i = 0; i < n_fmt; ++i)  fmt[i]  = other.fmt[i]->Clone();
	if (other.gt != nullptr) gt = new yon_gt(*other.gt);
	if (other.gt_sum != nullptr) gt_sum = new yon_gt_summary(*other.gt_sum);
	if (other.gt_sum_occ != nullptr) {
		gt_sum_occ = new yon_gt_summary[gt->n_o];
		for (int i = 0; i < gt->n_o; ++i) gt_sum_occ[i] = other.gt_sum_occ[i];
	}

	return(*this);
}

yon1_vnt_t::yon1_vnt_t(yon1_vnt_t&& other) noexcept :
	is_dirty(other.is_dirty), is_loaded_gt(other.is_loaded_gt), controller(other.controller),
	n_base_ploidy(other.n_base_ploidy), n_alleles(other.n_alleles),
	n_flt(other.n_flt), n_fmt(other.n_fmt), n_info(other.n_info),
	info_pid(other.info_pid), flt_pid(other.flt_pid), fmt_pid(other.fmt_pid),
	m_fmt(other.m_fmt), m_info(other.m_info), m_allele(other.m_allele),
	qual(other.qual), rid(other.rid), pos(other.pos), name(std::move(other.name)),
	info_hdr(other.info_hdr), fmt_hdr(other.fmt_hdr), flt_hdr(other.flt_hdr),
	info_map(other.info_map), fmt_map(other.fmt_map), flt_map(other.flt_map),
	alleles(nullptr),
	gt(nullptr), gt_sum(nullptr), gt_sum_occ(nullptr),
	info(nullptr), fmt(nullptr)
{
	std::swap(alleles, other.alleles);
	std::swap(gt, other.gt);
	std::swap(gt_sum, other.gt_sum);
	std::swap(gt_sum_occ, other.gt_sum_occ);
	if (n_info) {
		info = new PrimitiveContainerInterface*[n_info];
		for (int i = 0; i < n_info; ++i) info[i] = other.info[i]->Move();
	}

	if (n_fmt) {
		fmt = new PrimitiveGroupContainerInterface*[n_fmt];
		for (int i = 0; i < n_fmt; ++i) fmt[i] = other.fmt[i]->Move();
	}
}
//yon1_vnt_t& operator=(yon1_vnt_t&& other) noexcept = delete; // temporarily deleted

yon1_vnt_t::~yon1_vnt_t() {
	delete [] alleles;
	delete gt; delete gt_sum;
	delete [] gt_sum_occ;
	for (int i = 0; i < n_info; ++i) delete info[i];
	for (int i = 0; i < n_fmt;  ++i) delete fmt[i];
	delete [] info; delete [] fmt;
}

bool yon1_vnt_t::UpdateBase(const bcf1_t* record) {
	n_alleles = record->n_allele;
	qual = record->qual;
	rid  = record->rid;
	pos  = record->pos;
	name = std::string(record->d.id);
	delete [] alleles;
	alleles = new yon_allele[n_alleles];

	// Fix for the special case when ALT is not encoded
	if (this->n_alleles == 1) {
		delete [] alleles;
		this->n_alleles = 2;
		this->alleles = new yon_allele[n_alleles];
		this->alleles[0] = record->d.allele[0];
		this->alleles[1].allele    = new char[1];
		this->alleles[1].allele[0] = '.';
		this->alleles[1].l_allele  = 1;
	} else {
		for (uint32_t i = 0; i < this->n_alleles; ++i) {
			this->alleles[i] = record->d.allele[i];
		}
	}

	if (this->n_alleles == 2) this->controller.biallelic = true;
	this->controller.simple_snv = (this->alleles[0].length() == 1 && this->alleles[1].length() == 1);

	return true;
}

std::string yon1_vnt_t::GetAlleleString(void) const {
	std::string ret = this->alleles[0].ToString();
	for (uint32_t i = 1; i < this->n_alleles; ++i)
		ret += "," + this->alleles[i].ToString();
	return(ret);
}

bool yon1_vnt_t::EvaluateSummary(bool lazy_evaluate) {
	assert(this->gt != nullptr);
	assert(this->gt->rcds != nullptr);

	if (this->gt_sum != nullptr)
		return true;

	this->gt_sum = new yon_gt_summary(this->gt->m, this->gt->n_allele);
	*this->gt_sum += *this->gt;
	if (lazy_evaluate) this->gt_sum->LazyEvaluate();
	return true;
}

bool yon1_vnt_t::EvaluateOccSummary(bool lazy_evaluate) {
	assert(this->gt != nullptr);
	assert(this->gt->rcds != nullptr);

	if (this->gt_sum_occ != nullptr)
		return false;

	if (this->gt->n_o == 0)
		return false;

	this->gt_sum_occ = new yon_gt_summary[this->gt->n_o];
	for (int i = 0; i < this->gt->n_o; ++i) {
		this->gt_sum_occ[i].Setup(this->gt->m, this->gt->n_allele);
		this->gt_sum_occ[i].Add(*this->gt, this->gt->n_i_occ[i], this->gt->d_occ[i]);
		if (lazy_evaluate) this->gt_sum_occ[i].LazyEvaluate();
	}

	return true;
}


bool yon1_vnt_t::EvaluateOcc(yon_occ& occ) {
	if (occ.row_names.size() == 0)
		return false;

	if (this->gt == nullptr)
		return false;

	this->gt->n_o     = occ.occ.size();
	this->gt->n_i_occ = new uint32_t[this->gt->n_o];
	this->gt->d_occ   = new yon_gt_rcd*[this->gt->n_o];

	uint32_t* cum_sums = new uint32_t[this->gt->n_o]; // Total cumulative genotypes observed.
	uint32_t* cum_sums_hit = new uint32_t[this->gt->n_o]; // Number of non-zero runs observed.
	memset(cum_sums, 0, sizeof(uint32_t)*this->gt->n_o);
	memset(cum_sums_hit, 0, sizeof(uint32_t)*this->gt->n_o);

	for (uint32_t i = 0; i < this->gt->n_o; ++i) {
		this->gt->d_occ[i]   = new yon_gt_rcd[std::min(this->gt->n_i, occ.cum_sums[i])];
		this->gt->n_i_occ[i] = 0;
	}

	// Iterate over available gt rcds.
	for (uint32_t j = 0; j < this->gt->n_i; ++j) {
		// Iterate over available groupings in transposed occ.
		for (uint32_t i = 0; i < this->gt->n_o; ++i) {
			const uint32_t to   = occ.vocc[cum_sums[i] + this->gt->rcds[j].run_length][i];
			const uint32_t from = occ.vocc[cum_sums[i]][i];
			if (to - from != 0) {
				// Allocate memory for alleles.
				this->gt->d_occ[i][this->gt->n_i_occ[i]].allele = new uint8_t[this->gt->m];

				// Copy allelic data from recerence gt rcd.
				for (uint32_t k = 0; k < this->gt->m; ++k) {
					this->gt->d_occ[i][this->gt->n_i_occ[i]].allele[k] = this->gt->rcds[j].allele[k];
				}

				// Set run-length representation.
				this->gt->d_occ[i][this->gt->n_i_occ[i]].run_length = to - from;
				assert(this->gt->n_i_occ[i] < this->gt->n_i);
				++this->gt->n_i_occ[i];
				cum_sums_hit[i] += to - from;
			}
			cum_sums[i] += this->gt->rcds[j].run_length;
		}
	}
	// Debug assertions.
	for (int i = 0; i < this->gt->n_o; ++i) {
		assert(cum_sums[i] == this->gt->n_s);
		assert(cum_sums_hit[i] == occ.cum_sums[i]);
	}

	delete [] cum_sums;
	delete [] cum_sums_hit;

	return(true);
}

bool yon1_vnt_t::UsePackedRefAlt(void) const {
	if (this->controller.biallelic == false || this->controller.diploid == false)
		return false;

	if (std::regex_match(std::string(this->alleles[0].allele, this->alleles[0].l_allele), YON_REGEX_PACKED_ALLELES) &&
	   std::regex_match(std::string(this->alleles[1].allele, this->alleles[1].l_allele), YON_REGEX_PACKED_ALLELES)) {
		return true;
	}
	return false;
}

uint8_t yon1_vnt_t::PackRefAltByte(void) const {
	assert(this->UsePackedRefAlt());
	uint8_t ref_alt = 0; // start out with empty

	if (this->alleles[0].l_allele == 9 && strncmp(this->alleles[0].allele, "<NON_REF>", 9) == 0) {
		ref_alt ^= YON_ALLELE_NON_REF << 4;
	} else {
		switch(this->alleles[0].allele[0]) {
		case 'A': ref_alt ^= YON_ALLELE_A << 4; break;
		case 'T': ref_alt ^= YON_ALLELE_T << 4; break;
		case 'G': ref_alt ^= YON_ALLELE_G << 4; break;
		case 'C': ref_alt ^= YON_ALLELE_C << 4; break;
		case 'N': ref_alt ^= YON_ALLELE_N << 4; break;
		case '.': ref_alt ^= YON_ALLELE_MISS << 4; break;
		default:
			std::cerr << utility::timestamp("ERROR") << "Illegal SNV reference..." << std::endl;
			std::cerr << std::string(this->alleles[0].allele , this->alleles[0].l_allele) << std::endl;
			std::cerr << std::string(this->alleles[1].allele , this->alleles[1].l_allele) << std::endl;
			exit(1);
		}
	}

	if (this->alleles[1].l_allele == 9 && strncmp(this->alleles[1].allele, "<NON_REF>", 9) == 0) {
		ref_alt ^= YON_ALLELE_NON_REF << 0;
	} else {
		switch(this->alleles[1].allele[0]) {
		case 'A': ref_alt ^= YON_ALLELE_A << 0; break;
		case 'T': ref_alt ^= YON_ALLELE_T << 0; break;
		case 'G': ref_alt ^= YON_ALLELE_G << 0; break;
		case 'C': ref_alt ^= YON_ALLELE_C << 0; break;
		case 'N': ref_alt ^= YON_ALLELE_N << 0; break;
		case '.': ref_alt ^= YON_ALLELE_MISS << 0; break;
		default:
			std::cerr << utility::timestamp("ERROR") << "Illegal SNV alt..." << std::endl;
			exit(1);
		}
	}
	return(ref_alt);
}

bcf1_t* yon1_vnt_t::UpdateHtslibVcfRecord(bcf1_t* rec, bcf_hdr_t* hdr) const {
	rec->rid = this->rid;
	rec->pos = this->pos;
	bcf_update_id(hdr, rec, this->name.data());
	bcf_update_alleles_str(hdr, rec, this->GetAlleleString().data());
	if (std::isnan(this->qual)) bcf_float_set_missing(rec->qual);
	else rec->qual = this->qual;

	return(rec);
}

void yon1_vnt_t::OutputHtslibVcfInfo(bcf1_t* rec, bcf_hdr_t* hdr) {
	for (uint32_t j = 0; j < n_info; ++j) {
		if (info_hdr[j]->yon_type == YON_VCF_HEADER_FLAG) {
			bcf_update_info_flag(hdr, rec, info_hdr[j]->id.data(), NULL, 1);
		} else {
			info[j]->UpdateHtslibVcfRecordInfo(rec, hdr, info_hdr[j]->id);
		}
	}
}

void yon1_vnt_t::OutputHtslibVcfFormat(bcf1_t* rec,
						   bcf_hdr_t* hdr,
						   const bool display_genotypes,
						   yon_gt_rcd* external_exp) const
{
	if (n_fmt) {
		// Case when the only available FORMAT field is the GT field.
		if (n_fmt == 1 && is_loaded_gt &&
		   controller.gt_available &&
		   display_genotypes)
		{
			gt->ExpandExternal(external_exp);
			gt->d_exp = external_exp;
			gt->UpdateHtslibGenotypes(rec, hdr);
			gt->d_exp = nullptr;
		}
		// Case when there are > 1 Vcf Format fields and the GT field
		// is available.
		else if (n_fmt > 1 && is_loaded_gt &&
				controller.gt_available &&
				display_genotypes)
		{
			gt->ExpandExternal(external_exp);
			gt->d_exp = external_exp;

			gt->UpdateHtslibGenotypes(rec, hdr);
			for (uint32_t g = 1; g < n_fmt; ++g) {
				if (fmt_hdr[g]->yon_type == YON_VCF_HEADER_FLOAT)
					fmt[g]->UpdateHtslibVcfRecordFormatFloat(rec, hdr, fmt_hdr[g]->id);
				else if (fmt_hdr[g]->yon_type == YON_VCF_HEADER_INTEGER)
					fmt[g]->UpdateHtslibVcfRecordFormatInt32(rec, hdr, fmt_hdr[g]->id);
				else if (fmt_hdr[g]->yon_type == YON_VCF_HEADER_STRING || fmt_hdr[g]->yon_type == YON_VCF_HEADER_CHARACTER)
					fmt[g]->UpdateHtslibVcfRecordFormatString(rec, hdr, fmt_hdr[g]->id);
			}
			gt->d_exp = nullptr;
		}
		// All other cases.
		else {
			for (uint32_t g = 0; g < n_fmt; ++g) {
				if (fmt_hdr[g]->yon_type == YON_VCF_HEADER_FLOAT)
					fmt[g]->UpdateHtslibVcfRecordFormatFloat(rec, hdr, fmt_hdr[g]->id);
				else if (fmt_hdr[g]->yon_type == YON_VCF_HEADER_INTEGER)
					fmt[g]->UpdateHtslibVcfRecordFormatInt32(rec, hdr, fmt_hdr[g]->id);
				else if (fmt_hdr[g]->yon_type == YON_VCF_HEADER_STRING || fmt_hdr[g]->yon_type == YON_VCF_HEADER_CHARACTER)
					fmt[g]->UpdateHtslibVcfRecordFormatString(rec, hdr, fmt_hdr[g]->id);
			}
		}

	}
}

void yon1_vnt_t::OutputHtslibVcfFilter(bcf1_t* rec, bcf_hdr_t* hdr) const {
	if (n_flt) {
		for (uint32_t k = 0; k < n_flt; ++k) {
			int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, flt_hdr[k]->id.data());
			bcf_update_filter(hdr, rec, &tmpi, 1);
		}
	}
}

void yon1_vnt_t::ToVcfString(const yon_vnt_hdr_t& header,
		   yon_buffer_t& buffer,
		   uint32_t display,
		   yon_gt_rcd* external_rcd) const
{
	buffer.Add(header.contigs_[rid].name.data(), header.contigs_[rid].name.size());
	buffer += '\t';
	buffer.AddReadble(pos+1);
	buffer += '\t';
	if (name.size()) buffer += name;
	else buffer += '.';
	buffer += '\t';
	buffer.Add(alleles[0].allele, alleles[0].l_allele);
	buffer += '\t';
	buffer.Add(alleles[1].allele, alleles[1].l_allele);
	for (int i = 2; i < n_alleles; ++i) {
		buffer += ',';
		buffer.Add(alleles[i].allele, alleles[i].l_allele);
	}
	buffer += '\t';
	if (std::isnan(qual)) buffer += '.';
	else buffer.AddReadble(qual);
	buffer += '\t';

	if (flt_hdr.size() == 0) buffer += '.';
	else {
		buffer.Add(flt_hdr[0]->id.data(), flt_hdr[0]->id.size());
		for (int j = 1; j < flt_hdr.size(); ++j) {
			buffer += ';';
			buffer.Add(flt_hdr[j]->id.data(), flt_hdr[j]->id.size());
		}
	}
	buffer += '\t';

	if (n_info) {
		for (int j = 0; j < n_info; ++j) {
			if (j != 0) buffer += ';';
			if (info_hdr[j]->yon_type == YON_VCF_HEADER_FLAG) {
				buffer += info_hdr[j]->id;
			} else {
				buffer += info_hdr[j]->id;
				buffer += '=';
				info[j]->ToVcfString(buffer);
			}
		}
	} else {
		buffer += '.';
	}

	if (n_fmt) {
		buffer += '\t';
		buffer += fmt_hdr[0]->id;
		for (int j = 1; j < n_fmt; ++j) {
			buffer += ':';
			buffer += fmt_hdr[j]->id;
		}
		buffer += '\t';

		if (n_fmt == 1 &&
		   controller.gt_available &&
		   (display & YON_BLK_BV_GT))
		{
			gt->ExpandExternal(external_rcd);
			gt->d_exp = external_rcd;

			// Iterate over samples and print FORMAT:GT value in Vcf format.
			gt->d_exp[0].PrintVcf(buffer, gt->m);
			for (uint32_t s = 1; s < header.GetNumberSamples(); ++s) {
				buffer += '\t';
				gt->d_exp[s].PrintVcf(buffer, gt->m);
			}

			gt->d_exp = nullptr;
		}
		// Case when there are > 1 Vcf Format fields and the GT field
		// is available.
		else if (n_fmt > 1 && is_loaded_gt &&
				controller.gt_available &&
				(display & YON_BLK_BV_GT))
		{
			gt->ExpandExternal(external_rcd);
			gt->d_exp = external_rcd;

			gt->d_exp[0].PrintVcf(buffer, gt->m);
			for (uint32_t g = 1; g < n_fmt; ++g) {
				buffer += ':';
				fmt[g]->ToVcfString(buffer, 0);
			}
			for (uint32_t s = 1; s < header.GetNumberSamples(); ++s) {
				buffer += '\t';
				gt->d_exp[s].PrintVcf(buffer, gt->m);
				for (uint32_t g = 1; g < n_fmt; ++g) {
					buffer += ':';
					fmt[g]->ToVcfString(buffer, s);
				}
			}

			gt->d_exp = nullptr;
		}
		// All other cases.
		else {
			fmt[0]->ToVcfString(buffer, 0);
			for (uint32_t g = 1; g < n_fmt; ++g) {
				buffer += ':';
				fmt[g]->ToVcfString(buffer, 0);
			}

			for (uint32_t s = 1; s < header.GetNumberSamples(); ++s) {
				buffer += '\t';
				fmt[0]->ToVcfString(buffer, s);
				for (uint32_t g = 1; g < n_fmt; ++g) {
					buffer += ':';
					fmt[g]->ToVcfString(buffer, s);
				}
			}
		}
	}

	buffer += '\n';
}

bool yon1_vnt_t::AddInfoFlag(const std::string& tag, yon_vnt_hdr_t& header) {
	const YonInfo* info_tag = header.GetInfo(tag);
	assert(info_tag != nullptr);

	int32_t offset = GetInfoOffset(tag);
	if (offset >= 0) {
		delete info[offset];
		info[offset] = new PrimitiveContainer<int8_t>();
	} else {
		info[n_info++] = new PrimitiveContainer<int8_t>();
		info_hdr.push_back(info_tag);
	}
	return true;
}

bool yon1_vnt_t::AddGenotypeStatistics(yon_vnt_hdr_t& header, const bool replace_existing) {
	if (gt_sum == nullptr) {
		std::cerr << utility::timestamp("LOG") << "Failed to update record. Genotype statistics have not been processed!" << std::endl;
		return false;
	}

	if (gt_sum->d == nullptr) {
		std::cerr << utility::timestamp("LOG") << "Failed to update record. Genotype statistics have not been lazy evaluated!" << std::endl;
		return false;
	}

	if (n_info + 10 > m_info) {
		m_info = std::max(n_info + 10, m_info + 10);
		PrimitiveContainerInterface** old = info;
		info = new PrimitiveContainerInterface*[m_info];
		for (int i = 0; i < n_info; ++i) info[i] = old[i];
		delete [] old;
	}

	this->AddInfo("NM", header, gt_sum->d->nm);
	this->AddInfo("NPM", header, gt_sum->d->npm);
	this->AddInfo("AN", header, gt_sum->d->an);
	this->AddInfo("HWE_P", header, gt_sum->d->hwe_p);

	if (gt_sum->d->n_ac_af > 2) {
		this->AddInfo("AC", header, &gt_sum->d->ac[2], gt_sum->d->n_ac_af - 2);
		this->AddInfo("AF", header, &gt_sum->d->af[2], gt_sum->d->n_ac_af - 2);
	}

	if (this->n_alleles > 2) this->AddInfoFlag("MULTI_ALLELIC", header);

	// Special case for AC_P
	PrimitiveContainer<uint32_t>* ac_p = new PrimitiveContainer<uint32_t>();
	ac_p->resize(n_base_ploidy * (gt_sum->d->n_ac_af - 2));

	int j = 0;
	for (uint32_t p = 0; p < n_base_ploidy; ++p) {
		for (uint32_t i = 2; i < gt_sum->d->n_ac_af; ++i, ++j) {
			ac_p->at(j) = gt_sum->alleles_strand[p][i];
			++(*ac_p);
		}
	}
	this->AddInfo("AC_P", header, ac_p);

	if (n_base_ploidy == 2) {
		this->AddInfo("F_PIC", header, gt_sum->d->f_pic);
		this->AddInfo("HET", header, gt_sum->d->heterozygosity);
	}

	assert(n_info == info_hdr.size());
	return true;
}

bool yon1_vnt_t::AddGenotypeStatisticsOcc(yon_vnt_hdr_t& header, std::vector<std::string>& names) {
	if (names.size() == 0)
		return true;

	if (n_info + names.size()*10 + 3 > m_info) {
		m_info += names.size()*10 + 3;
		PrimitiveContainerInterface** old = info;
		info = new PrimitiveContainerInterface*[m_info];
		for (int i = 0; i < n_info; ++i) info[i] = old[i];
		delete [] old;
	}

	if (EvaluatePopGenOcc(header) == false) {
		//std::cerr << utility::timestamp("ERROR") << "Failed to evaluate popgen for occ..." << std::endl;
		return false;
	}

	for (int i = 0; i < names.size(); ++i) {
		this->AddInfo(names[i]+"_NM",    header, gt_sum_occ[i].d->nm);
		this->AddInfo(names[i]+"_NPM",   header, gt_sum_occ[i].d->npm);
		this->AddInfo(names[i]+"_AN",    header, gt_sum_occ[i].d->an);
		this->AddInfo(names[i]+"_HWE_P", header, gt_sum_occ[i].d->hwe_p);

		if (gt_sum_occ[i].d->n_ac_af > 2) {
			this->AddInfo(names[i]+"_AC", header, &gt_sum_occ[i].d->ac[2], gt_sum_occ[i].d->n_ac_af - 2);
			this->AddInfo(names[i]+"_AF", header, &gt_sum_occ[i].d->af[2], gt_sum_occ[i].d->n_ac_af - 2);
		}

		if (this->n_alleles > 2) this->AddInfoFlag(names[i]+"_MULTI_ALLELIC", header);

		// Special case for AC_P
		PrimitiveContainer<uint32_t>* ac_p = new PrimitiveContainer<uint32_t>();
		ac_p->resize(n_base_ploidy * (gt_sum_occ[i].d->n_ac_af - 2));

		int j = 0;
		for (uint32_t p = 0; p < n_base_ploidy; ++p) {
			for (uint32_t k = 2; k < gt_sum_occ[i].d->n_ac_af; ++k, ++j) {
				ac_p->at(j) = gt_sum_occ->alleles_strand[p][k];
				++(*ac_p);
			}
		}
		this->AddInfo(names[i]+"_AC_P", header, ac_p);

		if (n_base_ploidy == 2) {
			this->AddInfo(names[i]+"_F_PIC", header, gt_sum_occ[i].d->f_pic);
			this->AddInfo(names[i]+"_HET",   header, gt_sum_occ[i].d->heterozygosity);
		}
	}

	assert(n_info == info_hdr.size());
	return true;
}

bool yon1_vnt_t::EvaluatePopGenOcc(yon_vnt_hdr_t& header) {
	if (this->gt == nullptr) return false;
	if (this->gt->n_o == 0) return false;
	if (gt->n_allele != 2 || gt->m != 2) return false;

	// Number of genotypes across populations.
	uint32_t n_gt_all = 0;
	uint32_t n_gt_ref = 0, n_gt_alt, n_gt_het = 0;

	float* obs = new float[this->gt->n_o];
	float* exp = new float[this->gt->n_o];
	uint32_t* n_s = new uint32_t[this->gt->n_o];

	for (int i = 0; i < this->gt->n_o; ++i) {
		const uint32_t het  = this->gt_sum_occ[i].gt[2][3].n_cnt + this->gt_sum_occ[i].gt[3][2].n_cnt;
		const uint32_t n_gt = het + this->gt_sum_occ[i].gt[2][2].n_cnt + this->gt_sum_occ[i].gt[3][3].n_cnt;

		if (n_gt == 0) {
			obs[i] = 0; exp[i] = 0; n_s[i] = 0;
			continue;
		}
		n_gt_all += n_gt;
		n_gt_ref += this->gt_sum_occ[i].gt[2][2].n_cnt;
		n_gt_alt += this->gt_sum_occ[i].gt[3][3].n_cnt;
		n_gt_het += het;

		const float p  = (2*(float)this->gt_sum_occ[i].gt[2][2].n_cnt + het) / (2*n_gt);
		const float q  = (2*(float)this->gt_sum_occ[i].gt[3][3].n_cnt + het) / (2*n_gt);
		const float pq2 = 2*p*q;


		// F = (Hexp - Hobs) / Hexp
		const float Hobs = (float)het/n_gt; // observed heterozygosities in populations
		const float Hexp = (float)pq2;      // expected heterozygosities in populations

		/*
		std::cerr << "occ" << i << "\t" <<
		this->gt_sum_occ[i].gt[2][2].n_cnt << "\t" <<
		this->gt_sum_occ[i].gt[2][3].n_cnt << "\t" <<
		this->gt_sum_occ[i].gt[3][2].n_cnt << "\t" <<
		this->gt_sum_occ[i].gt[3][3].n_cnt << "\t" <<
		p << "\t" << (1-p) << "\t" <<
		Hexp << "\t" << Hobs << "\t" << "F=" <<
		(Hexp - Hobs) / Hexp << "\t" << Hobs*n_gt << "/" << Hexp*n_gt << "/" << n_gt << std::endl;
		*/

		obs[i] = Hobs;
		exp[i] = Hexp;
		n_s[i] = n_gt;
	}

	if (n_gt_all == 0)
		return false;

	float exp_sum = 0, obs_sum = 0;
	for (int i = 0; i < this->gt->n_o; ++i) {
		exp_sum += exp[i] * n_s[i];
		obs_sum += obs[i] * n_s[i];
	}

	const float h_i = obs_sum / n_gt_all;
	const float h_s = exp_sum / n_gt_all;
	// p-bar (frequency of allele A) over all populations
	// q-bar
	const float p_bar = (2*(float)n_gt_ref + n_gt_het) / (2*n_gt_all);
	const float q_bar = (2*(float)n_gt_alt + n_gt_het) / (2*n_gt_all);
	const float h_t   =  2*p_bar*q_bar;

	//std::cerr << "pbar=" << p_bar << " qbar=" << q_bar << " with " << n_gt_all << std::endl;
	//std::cerr << "sums: " << exp_sum << ", " << obs_sum << std::endl;
	//std::cerr << "hi=" << h_i << ", hs=" << h_s << ", ht=" << h_t << std::endl;

	// Fstatistics
	// FIS [= (HS - HI)/HS]
	//std::cerr << "FIS=" << (h_s - h_i) / h_s << std::endl;
	// FST [= (HT - HS)/HT]
	//std::cerr << "FST=" << (h_t - h_s) / h_t << std::endl;
	// FIT [= (HT - HI)/HT]
	//std::cerr << "FIT=" << (h_t - h_i) / h_t << std::endl;
	float fis = 0, fst = 0, fit = 0;

	if (h_s != 0) fis = (h_s - h_i) / h_s;
	if (h_t != 0) fst = (h_t - h_s) / h_t;
	if (h_t != 0) fit = (h_t - h_i) / h_t;

	this->AddInfo("FIS", header, fis);
	this->AddInfo("FST", header, fst);
	this->AddInfo("FIT", header, fit);

	delete [] exp; delete [] obs; delete [] n_s;

	return true;
}

PrimitiveContainerInterface* yon1_vnt_t::GetInfo(const std::string& name) const {
	std::unordered_map<std::string, uint32_t>::const_iterator it = info_map.find(name);
	if (it != info_map.end()) return(info[it->second]);
	return(nullptr);
}

PrimitiveGroupContainerInterface* yon1_vnt_t::GetFmt(const std::string& name) const {
	std::unordered_map<std::string, uint32_t>::const_iterator it = info_map.find(name);
	if (it != info_map.end()) return(fmt[it->second]);
	return(nullptr);
}

const VcfFilter* yon1_vnt_t::GetFlt(const std::string& name) const {
	std::unordered_map<std::string, uint32_t>::const_iterator it = flt_map.find(name);
	if (it != flt_map.end()) return(flt_hdr[it->second]);
	return(nullptr);
}

int32_t yon1_vnt_t::GetInfoOffset(const std::string& name) const {
	std::unordered_map<std::string, uint32_t>::const_iterator it = info_map.find(name);
	if (it != info_map.end()) return(it->second);
	return(-1);
}

int32_t yon1_vnt_t::GetFormatOffset(const std::string& name) const {
	std::unordered_map<std::string, uint32_t>::const_iterator it = info_map.find(name);
	if (it != info_map.end()) return(it->second);
	return(-1);
}

int32_t yon1_vnt_t::GetFilterOffset(const std::string& name) const {
	std::unordered_map<std::string, uint32_t>::const_iterator it = flt_map.find(name);
	if (it != flt_map.end()) return(it->second);
	return(-1);
}

/****************************
*  Ts/Tv statistics
****************************/
yon_stats_sample::yon_stats_sample(void) :
	n_ins(0), n_del(0), n_singleton(0), n_ts(0), n_tv(0), ts_tv_ratio(0), ins_del_dist(new uint64_t[512])
{
	for (uint32_t i = 0; i < 9; ++i) {
		this->base_conv[i] = new uint64_t[9];
		memset(&this->base_conv[i][0], 0, sizeof(uint64_t)*9);
	}
	memset(ins_del_dist, 0, 512*sizeof(uint64_t));
}

yon_stats_sample::~yon_stats_sample(void) {
	for (uint32_t i = 0; i < 9; ++i) {
		delete [] this->base_conv[i];
	}
	delete [] ins_del_dist;
}

yon_stats_sample& yon_stats_sample::operator+=(const yon_stats_sample& other) {
	this->n_ins += other.n_ins;
	this->n_del += other.n_del;
	this->n_singleton += other.n_singleton;
	this->n_ts += other.n_ts;
	this->n_tv += other.n_tv;
	for (int i = 0; i < 9; ++i) {
		for (uint32_t j = 0; j < 9; ++j) {
			this->base_conv[i][j] += other.base_conv[i][j];
		}
	}

	for (int i = 0; i < 512; ++i)
		this->ins_del_dist[i] += other.ins_del_dist[i];

	return(*this);
}

bool yon_stats_sample::LazyEvalute(void) {
	// Transversions: A->C, C->A, T->G, G->T, A->T, T->A, C->G, G->C
	// Transitions:   A->G, G->A, C->T, T->C
	this->n_tv = this->base_conv[YON_GT_TSTV_A][YON_GT_TSTV_C] + this->base_conv[YON_GT_TSTV_C][YON_GT_TSTV_A] +
				 this->base_conv[YON_GT_TSTV_T][YON_GT_TSTV_G] + this->base_conv[YON_GT_TSTV_G][YON_GT_TSTV_T] +
				 this->base_conv[YON_GT_TSTV_A][YON_GT_TSTV_T] + this->base_conv[YON_GT_TSTV_T][YON_GT_TSTV_A] +
				 this->base_conv[YON_GT_TSTV_C][YON_GT_TSTV_G] + this->base_conv[YON_GT_TSTV_G][YON_GT_TSTV_C];
	this->n_ts = this->base_conv[YON_GT_TSTV_A][YON_GT_TSTV_G] + this->base_conv[YON_GT_TSTV_G][YON_GT_TSTV_A] +
				 this->base_conv[YON_GT_TSTV_T][YON_GT_TSTV_C] + this->base_conv[YON_GT_TSTV_C][YON_GT_TSTV_T];

	if (this->n_tv == 0) this->ts_tv_ratio = 0;
	else this->ts_tv_ratio = ((double)this->n_ts / this->n_tv);

	this->n_ins = this->base_conv[YON_GT_TSTV_A][YON_GT_TSTV_INS] + this->base_conv[YON_GT_TSTV_T][YON_GT_TSTV_INS] +
				  this->base_conv[YON_GT_TSTV_C][YON_GT_TSTV_INS] + this->base_conv[YON_GT_TSTV_C][YON_GT_TSTV_INS] +
				  this->base_conv[YON_GT_TSTV_UNKNOWN][YON_GT_TSTV_INS];
	this->n_del = this->base_conv[YON_GT_TSTV_UNKNOWN][YON_GT_TSTV_DEL];

	return true;
}

yon_buffer_t& yon_stats_sample::ToJsonString(yon_buffer_t& buffer, const std::string& sample_name) const {
	buffer +=  "\"" + sample_name + "\":{";
	buffer +=  "\"n_ins\":" + std::to_string(this->n_ins);
	buffer += ",\"n_del\":" + std::to_string(this->n_del);
	buffer += ",\"n_singleton\":" + std::to_string(this->n_singleton);
	buffer += ",\"n_ts\":"  + std::to_string(this->n_ts);
	buffer += ",\"n_tv\":"  + std::to_string(this->n_tv);
	buffer += ",\"ts_tv\":" + std::to_string(this->ts_tv_ratio);
	buffer += ",\"conv\":[";
	for (uint32_t i = 0; i < 9; ++i) {
		if (i != 0) buffer += ',';
		buffer += '[';
		buffer.AddReadble((uint64_t)this->base_conv[i][0]);
		for (uint32_t j = 1; j < 9; ++j) {
			buffer += ',';
			buffer.AddReadble((uint64_t)this->base_conv[i][j]);
		}
		buffer += ']';
	}

	buffer += ",\"in_dist\":[";
	//buffer.AddReadble(this->ins_del_dist[0]);
	for (int i = 1; i < 512; i+=2) {
		if (i != 1) buffer += ',';
		buffer.AddReadble(this->ins_del_dist[i]);
	}
	buffer += "]";
	buffer += ",\"del_dist\":[";
	//buffer.AddReadble(this->ins_del_dist[0]);
	for (int i = 2; i < 512; i+=2) {
		if (i != 2) buffer += ',';
		buffer.AddReadble(this->ins_del_dist[i]);
	}
	buffer += ']';

	buffer += ']';
	buffer += '}';
	return(buffer);
}

void yon_stats_sample::reset(void) {
	n_ins = 0, n_del = 0, n_singleton = 0;
	n_ts = 0, n_tv = 0; ts_tv_ratio = 0;

	for (int i = 0; i < 9; ++i)
		memset(&this->base_conv[i][0], 0, sizeof(uint64_t)*9);

	for (int i = 0; i < 512; ++i) this->ins_del_dist[i] = 0;
}

yon_stats_tstv::yon_stats_tstv() : n_s(0), n_rcds(0), n_snp(0), n_mnp(0), n_ins(0), n_del(0), n_other(0), n_no_alts(0), n_singleton(0),
				   n_biallelic(0), n_multi_allele(0), n_multi_allele_snp(0), sample(nullptr), alt_count(new uint64_t[32]) { memset(alt_count, 0, sizeof(uint64_t)*32); }
yon_stats_tstv::yon_stats_tstv(const uint32_t n_samples) : n_s(n_samples), n_rcds(0), n_snp(0), n_mnp(0), n_ins(0), n_del(0), n_other(0), n_no_alts(0), n_singleton(0),
										   n_biallelic(0), n_multi_allele(0), n_multi_allele_snp(0), sample(new yon_stats_sample[n_samples]), alt_count(new uint64_t[32]) { memset(alt_count, 0, sizeof(uint64_t)*32); }
yon_stats_tstv::~yon_stats_tstv(void) { delete [] this->sample; delete [] this->alt_count; }

yon_stats_tstv& yon_stats_tstv::operator+=(const yon_stats_tstv& other) {
	assert(other.n_s == this->n_s);
	this->n_rcds  += other.n_rcds;
	this->n_snp   += other.n_snp;
	this->n_mnp   += other.n_mnp;
	this->n_ins   += other.n_ins;
	this->n_del   += other.n_del;
	this->n_other += other.n_other;
	this->n_no_alts   += other.n_no_alts;
	this->n_singleton += other.n_singleton;
	this->n_biallelic += other.n_biallelic;
	this->n_multi_allele += other.n_multi_allele;
	this->n_multi_allele_snp += other.n_multi_allele_snp;

	for (int i = 0; i < this->n_s; ++i) this->sample[i] += other.sample[i];
	for (int i = 0; i < 32; ++i) this->alt_count[i] = other.alt_count[i];

	return(*this);
}

yon_stats_tstv& yon_stats_tstv::Add(const yon_stats_tstv& other, const yon_gt_ppa* ppa) {
	if (ppa == nullptr)
		return(*this += other);

	assert(ppa->n_s == this->n_s);
	assert(other.n_s == this->n_s);

	this->n_rcds += other.n_rcds;
	this->n_snp += other.n_snp;
	this->n_mnp += other.n_mnp;
	this->n_ins += other.n_ins;
	this->n_del += other.n_del;
	this->n_other += other.n_other;
	this->n_no_alts += other.n_no_alts;
	this->n_singleton += other.n_singleton;
	this->n_biallelic += other.n_biallelic;
	this->n_multi_allele += other.n_multi_allele;
	this->n_multi_allele_snp += other.n_multi_allele_snp;

	for (int i = 0; i < this->n_s; ++i) this->sample[(*ppa)[i]] += other.sample[i];
	for (int i = 0; i < 32; ++i) this->alt_count[i] = other.alt_count[i];

	return(*this);
}

void yon_stats_tstv::SetSize(const uint32_t n_samples) {
	delete [] sample;
	n_s = n_samples;
	sample = new yon_stats_sample[n_samples];
}

yon_buffer_t& yon_stats_tstv::ToJsonString(yon_buffer_t& buffer, const std::vector<std::string>& sample_names) const {
	buffer += '{';
	buffer +=  "\"VI\":{";
	buffer +=  "\"n_samples\":";
	buffer.AddReadble(this->n_s);
	buffer +=  ",\"n_records\":";
	buffer.AddReadble(this->n_rcds);
	buffer +=  ",\"n_biallelic\":";
	buffer.AddReadble(this->n_biallelic);
	buffer +=  ",\"n_del\":";
	buffer.AddReadble(this->n_del);
	buffer +=  ",\"n_ins\":";
	buffer.AddReadble(this->n_ins);
	buffer +=  ",\"n_snp\":";
	buffer.AddReadble(this->n_snp);
	buffer +=  ",\"n_mnp\":";
	buffer.AddReadble(this->n_mnp);
	buffer +=  ",\"n_other\":";
	buffer.AddReadble(this->n_other);
	buffer +=  ",\"n_no_alts\":";
	buffer.AddReadble(this->n_no_alts);
	buffer +=  ",\"n_multi_allele\":";
	buffer.AddReadble(this->n_multi_allele);
	buffer +=  ",\"n_multi_allele_snp\":";
	buffer.AddReadble(this->n_multi_allele_snp);
	buffer +=  ",\"n_singleton\":";
	buffer.AddReadble(this->n_singleton);
	buffer +=  ",\"n_alts\":[";
	for (int i = 0; i < 32; ++i) {
		if (i != 0) buffer += ',';
		buffer.AddReadble(this->alt_count[i]);
	}
	buffer += "]},\n";
	buffer += "\"PSI\":{\n";
	for (uint32_t i = 0; i < this->n_s; ++i) {
		if (i != 0) buffer += ",\n";
		this->sample[i].ToJsonString(buffer, sample_names[i]);
	}
	buffer += "\n}\n}\n";

	return(buffer);
}

bool yon_stats_tstv::GetEncodings(const yon1_vnt_t& rcd, yon_stats_tstv_obj& helper)
{
	// Update number of records.
	++this->n_rcds;

	if (rcd.n_alleles == 1) {
		std::cerr << utility::timestamp("LOG") << "Cannot have site with no ALT alleles described..." << std::endl;
		return false;
	}

	// Update count for target variant line type.
	// Case: variant site is multi-allelic.
	if (rcd.n_alleles > 2) {
		bool is_snp = true;
		for (uint32_t i = 0; i < rcd.gt->n_allele; ++i) {
			if (rcd.alleles[i].l_allele != 1) {
				is_snp = false;
				break;
			}
		}

		if (is_snp) ++this->n_multi_allele_snp;
		else ++this->n_multi_allele;
	}
	// Case: variant site is bi-allelic.
	else if (rcd.n_alleles == 2) {
		++this->n_biallelic;
	}

	memset(helper.non_ref_encodings, 1, sizeof(uint8_t)*(rcd.gt->n_allele + 2));
	memset(helper.non_ref_encodings, 0, sizeof(uint8_t)*3);

	// For SNV to SNV or insertion. It is not possible to have a deletion
	// if the reference value is represented as a SNV.
	// If the reference allele is of a single character wide then
	// we assume the reference site is a single base. This implicitly
	// means there cannot be any encodings for deletions.
	if (rcd.alleles[0].size() == 1) {
		// Encode alleles.
		helper.allele_encodings[0] = YON_GT_TSTV_MISS;
		helper.allele_encodings[1] = YON_GT_TSTV_EOV;
		memset(helper.b_size, 0, sizeof(int32_t)*2);
		bool has_insert = false;

		// Iterate over available alleles.
		for (uint32_t i = 2; i < rcd.gt->n_allele + 2; ++i) {
			if (rcd.alleles[i - 2].l_allele == 1) {
				helper.allele_encodings[i] = YON_STATS_TSTV_LOOKUP[rcd.alleles[i - 2].allele[0]];
				helper.b_size[i] = 0;
			} else {
			if (rcd.alleles[i - 2].l_allele > 1 &&
				   std::regex_match(rcd.alleles[i - 2].ToString(), YON_REGEX_CANONICAL_BASES))
				{
					//std::cerr << "is insertion: " << rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele << std::endl;
					helper.allele_encodings[i] = YON_GT_TSTV_INS;
					helper.b_size[i] = std::min(255, rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele);
					has_insert = true;
				} else {
					helper.allele_encodings[i] = YON_GT_TSTV_UNKNOWN;
					helper.b_size[i] = 0;
				}
			}
		}

		// Ascertain that the reference allele is valid for this purpose.
		if (helper.allele_encodings[YON_GT_RCD_REF] > 3) {
			std::cerr << utility::timestamp("LOG") << "Bad reference allele: " << rcd.alleles[0].ToString() << std::endl;
			return false;
		}

		this->n_ins += has_insert;
	}
	// For insertion/deletion to SNV/insertion/deletion.
	else {
		// Encode alleles.
		helper.allele_encodings[0] = YON_GT_TSTV_MISS;
		helper.allele_encodings[1] = YON_GT_TSTV_EOV;
		helper.allele_encodings[2] = YON_GT_TSTV_UNKNOWN;
		const uint16_t& ref_length = rcd.alleles[0].l_allele;
		bool has_insert = false;
		bool has_deletion = false;

		memset(helper.b_size, 0, sizeof(int32_t)*3);

		// Iterate over available alleles.
		for (uint32_t i = 3; i < rcd.gt->n_allele + 2; ++i) {
			// If the target allele is a simple SNV.
			if (rcd.alleles[i - 2].l_allele == 1) {
				if (std::regex_match(rcd.alleles[i - 2].ToString(), YON_REGEX_CANONICAL_BASES)) {
					//std::cerr << "target is deletion: " << i - 2 << "/" << rcd.n_alleles << "; " << rcd.alleles[0].ToString() << "->" << rcd.alleles[i - 2].ToString() << " size: " << rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele << std::endl;
					helper.allele_encodings[i] = YON_GT_TSTV_DEL;
					helper.b_size[i] = std::max(-255, (int)rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele);
					has_deletion = true;
				} else {
					helper.allele_encodings[i] = YON_GT_TSTV_UNKNOWN;
					helper.b_size[i] = 0;
				}
			}
			// If the target allele length is shorter than the reference
			// allele length and is comprised of only canonical bases then
			// classify this allele as a deletion.
			else if (rcd.alleles[i - 2].l_allele < ref_length &&
					std::regex_match(rcd.alleles[i - 2].ToString(), YON_REGEX_CANONICAL_BASES))
			{
				//std::cerr << "target is deletion: " << i - 2 << "/" << rcd.n_alleles << "; " << rcd.alleles[0].ToString() << "->" << rcd.alleles[i - 2].ToString() << " size: " << rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele << std::endl;
				helper.allele_encodings[i] = YON_GT_TSTV_DEL;
				helper.b_size[i] = std::max(-255, (int)rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele);
				has_deletion = true;
			} else {
				if (rcd.alleles[i - 2].l_allele > ref_length &&
				   std::regex_match(rcd.alleles[i - 2].ToString(), YON_REGEX_CANONICAL_BASES))
				{
					//std::cerr << "is insertion: " << rcd.alleles[i - 2].ToString() << ": " << rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele << "/" << rcd.n_alleles << std::endl;
					helper.allele_encodings[i] = YON_GT_TSTV_INS;
					helper.b_size[i] = std::min(255, rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele);
					has_insert = true;
				} else {
					//std::cerr << "is same: " << rcd.alleles[i - 2].toString() << ": " << i - 2 << "/" << rcd.n_alleles << std::endl;
					helper.allele_encodings[i] = YON_GT_TSTV_UNKNOWN;
					helper.b_size[i] = 0;
				}
			}
		}

		this->n_ins += has_insert;
		this->n_del += has_deletion;
	}

	return true;
}

bool yon_stats_tstv::Update(const yon1_vnt_t& rcd, yon_gt_rcd** rcds) {
	if (rcd.is_loaded_gt == false)
		return false;

	if ((rcd.gt->eval_cont & YON_GT_UN_RCDS) == false) {
		bool eval = rcd.gt->Evaluate();
		if (eval == false) {
			std::cerr << "failed to evaluate" << std::endl;
			return false;
		}
	}

	assert(this->n_s == rcd.gt->n_s);
	yon_stats_tstv_obj helper(rcd.gt->n_allele + 2);
	if (this->GetEncodings(rcd, helper) == false)
		return false;

	// Update number of alt alleles.
	++this->alt_count[std::min(31, rcd.n_alleles - 1)];

	// Perform actual work.
	for (uint32_t i = 0; i < rcd.gt->n_s; ++i) {
		for (uint32_t j = 0; j < rcd.gt->m; ++j) {
			const uint8_t allele = YON_GT_RCD_ALLELE_UNPACK(rcds[i]->allele[j]);
			++this->sample[i].base_conv[helper.allele_encodings[YON_GT_RCD_REF]][helper.allele_encodings[allele]];
			++this->sample[i].ins_del_dist[(helper.b_size[allele] << 1) ^ (helper.b_size[allele] >> 31)];
			helper.n_non_ref += helper.non_ref_encodings[allele];
			helper.t_non_ref += helper.non_ref_encodings[allele] * i;
		}
	}

	if (helper.n_non_ref == 0) ++this->n_no_alts;
	else if (helper.n_non_ref == 1) {
		++this->n_singleton;
		++this->sample[helper.t_non_ref].n_singleton;
		assert(helper.t_non_ref < this->n_s);
	}

	return true;
}

bool yon_stats_tstv::Update(const yon1_vnt_t& rcd) {
	if (rcd.is_loaded_gt == false)
		return false;

	if ((rcd.gt->eval_cont & YON_GT_UN_RCDS) == false) {
		bool eval = rcd.gt->Evaluate();
		if (eval == false) {
			std::cerr << "failed to evaluate" << std::endl;
			return false;
		}
	}

	yon_stats_tstv_obj helper(rcd.gt->n_allele + 2);
	if (this->GetEncodings(rcd, helper) == false)
		return false;

	// Update number of alt alleles.
	++this->alt_count[std::min(31, rcd.n_alleles - 1)];

	if (rcd.gt->m == 2) this->UpdateDiploid(rcd.gt, helper);
	else this->UpdateNPloidy(rcd.gt, helper);

	if (helper.n_non_ref == 0) ++this->n_no_alts;
	else if (helper.n_non_ref == 1) {
		++this->n_singleton;
		++this->sample[helper.t_non_ref].n_singleton;
		assert(helper.t_non_ref < this->n_s);
	}

	return true;
}

void yon_stats_tstv::UpdateDiploid(const yon_gt* gt,
                                   yon_stats_tstv_obj& helper)
{
	uint32_t sample_offset = 0;
	for (uint32_t i = 0; i < gt->n_i; ++i) {
		if ((YON_GT_RCD_ALLELE_UNPACK(gt->rcds[i].allele[0]) == YON_GT_RCD_REF) && (YON_GT_RCD_ALLELE_UNPACK(gt->rcds[i].allele[1]) == YON_GT_RCD_REF)) {
			sample_offset += gt->rcds[i].run_length;
			continue;
		}

		for (uint32_t r = 0; r < gt->rcds[i].run_length; ++r, ++sample_offset) {
			for (uint32_t j = 0; j < gt->m; ++j) {
				const uint8_t allele = YON_GT_RCD_ALLELE_UNPACK(gt->rcds[i].allele[j]);
				++this->sample[sample_offset].base_conv[helper.allele_encodings[YON_GT_RCD_REF]][helper.allele_encodings[allele]];
				++this->sample[sample_offset].ins_del_dist[(helper.b_size[allele] << 1) ^ (helper.b_size[allele] >> 31)];
				helper.n_non_ref += helper.non_ref_encodings[allele];
				helper.t_non_ref += helper.non_ref_encodings[allele] * sample_offset;
			}
		}
	}
	assert(sample_offset == this->n_s);
}

void yon_stats_tstv::UpdateNPloidy(const yon_gt* gt,
                                   yon_stats_tstv_obj& helper)
{
	uint32_t sample_offset = 0;
	for (uint32_t i = 0; i < gt->n_i; ++i) {
		// If current run-length encoded object has all reference
		// template then continue.
		uint32_t n_refs = 0;
		for (uint32_t j = 0; j < gt->m; ++j)
			n_refs += (YON_GT_RCD_ALLELE_UNPACK(gt->rcds[i].allele[j]) == YON_GT_RCD_REF);

		if (n_refs == gt->m) {
			sample_offset += gt->rcds[i].run_length;
			continue;
		}

		// Iterate over samples in the current run-length encoded object.
		for (uint32_t r = 0; r < gt->rcds[i].run_length; ++r, ++sample_offset) {
			for (uint32_t j = 0; j < gt->m; ++j) {
				const uint8_t allele = YON_GT_RCD_ALLELE_UNPACK(gt->rcds[i].allele[j]);
				++this->sample[sample_offset].base_conv[helper.allele_encodings[YON_GT_RCD_REF]][helper.allele_encodings[allele]];
				++this->sample[sample_offset].ins_del_dist[(helper.b_size[allele] << 1) ^ (helper.b_size[allele] >> 31)];
				helper.n_non_ref += helper.non_ref_encodings[allele];
				helper.t_non_ref += helper.non_ref_encodings[allele] * sample_offset;
			}
		}
	}
	assert(sample_offset == this->n_s);
}

void yon_stats_tstv::reset(void) {
	n_rcds = 0;
	n_snp = 0, n_mnp = 0, n_ins = 0, n_del = 0, n_other = 0;
	n_no_alts = 0, n_biallelic = 0, n_multi_allele = 0, n_multi_allele_snp = 0, n_singleton = 0;
	for (int i = 0; i < this->n_s; ++i) this->sample[i].reset();
	for (int i = 0; i < 32; ++i) this->alt_count[i] = 0;
}

}
