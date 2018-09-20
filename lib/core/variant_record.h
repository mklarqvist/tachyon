#ifndef CORE_VARIANT_RECORD_H_
#define CORE_VARIANT_RECORD_H_

#include "containers/components/variant_block_footer.h"
#include "containers/primitive_container.h"
#include "containers/primitive_group_container.h"
#include "genotypes.h"

#include <cmath>

namespace tachyon{

/****************************
*  Core
****************************/
const char* const YON_REFALT_LOOKUP = "ATGC.XN";
#define YON_ALLELE_A       0
#define YON_ALLELE_T       1
#define YON_ALLELE_G       2
#define YON_ALLELE_C       3
#define YON_ALLELE_MISS    4
#define YON_ALLELE_NON_REF 5
#define YON_ALLELE_N       6

/****************************
*  TS/TV objects
****************************/
#define YON_GT_TSTV_A       0
#define YON_GT_TSTV_T       1
#define YON_GT_TSTV_G       2
#define YON_GT_TSTV_C       3
#define YON_GT_TSTV_UNKNOWN 4
#define YON_GT_TSTV_MISS    5
#define YON_GT_TSTV_EOV     6
#define YON_GT_TSTV_INS     7
#define YON_GT_TSTV_DEL     8

// ASCII positions for A,T,G,C are filled with their respective
// encodings and the background is filled with YON_GT_TSTV_UNKNOWN.
const uint8_t YON_STATS_TSTV_LOOKUP[256] =
{4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,YON_GT_TSTV_A,4,YON_GT_TSTV_C,
 4,4,4,YON_GT_TSTV_G,4,4,4,4,4,4,4,4,4,4,4,4,YON_GT_TSTV_T,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};


struct yon_allele {
public:
	yon_allele() : l_allele(0), allele(nullptr){}
	yon_allele(const std::string& string) : l_allele(string.size()), allele(new char[l_allele]){ memcpy(allele, string.data(), l_allele); }
	yon_allele& operator=(const std::string& value){ delete [] allele; l_allele = value.size(); allele = new char[l_allele]; memcpy(allele, value.data(), l_allele); return(*this); }
	yon_allele& operator=(const char* value){ delete [] allele; l_allele = strlen(value); allele = new char[l_allele]; memcpy(allele, value, l_allele); return(*this); }
	yon_allele& operator=(const char& value){ delete [] allele; l_allele = 1; allele = new char[l_allele]; allele[0] = value; return(*this); }

	yon_allele(const yon_allele& other) : l_allele(other.l_allele), allele(new char[l_allele]){ memcpy(allele, other.allele, l_allele); }
	yon_allele(yon_allele&& other) noexcept : l_allele(other.l_allele), allele(nullptr){ std::swap(allele, other.allele); other.l_allele = 0; }

	yon_allele& operator=(const yon_allele& other){
		delete [] allele;
		l_allele = other.l_allele;
		allele = new char[l_allele];
		memcpy(allele, other.allele, l_allele);
		return(*this);
	}

	yon_allele& operator=(yon_allele&& other) noexcept {
		if(this == &other){
			// precautions against self-moves
			return *this;
		}
		l_allele = other.l_allele;
		std::swap(allele, other.allele);
		other.l_allele = 0;
		return(*this);
	}

	~yon_allele(){ delete [] allele; }

	yon_allele& ParseFromBuffer(const char* const in){
		delete [] this->allele;
		this->l_allele = *reinterpret_cast<const uint16_t* const>(in);
		this->allele = new char[this->l_allele];
		memcpy(this->allele, &in[sizeof(uint16_t)], this->l_allele);
		return(*this);
	}

	inline const uint16_t& size(void) const{ return(this->l_allele); }
	inline const uint16_t& length(void) const{ return(this->l_allele); }
	inline const std::string ToString(void) const{ return(std::string(this->allele, this->l_allele)); }

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const yon_allele& entry);

public:
	uint16_t l_allele;
	char*    allele;
};

/**<
 * Primary evaluated record of a Tachyon variant. Construction is done
 * outside of this definition. Evaluation into this object is relatively
 * inexpensive in isolation but quiet expensive when considered
 * across the an entire dataset as this structure needs be evaluated many
 * millions of times over many different containers.
 */
struct yon1_vnt_t {
public:
	yon1_vnt_t() :
		is_dirty(false), is_loaded_gt(false), n_base_ploidy(0), n_alleles(0),
		n_flt(0), n_fmt(0), n_info(0), info_pid(-1), flt_pid(-1), fmt_pid(-1),
		m_fmt(0), m_info(0), m_allele(0),
		qual(NAN), rid(0), pos(0),
		alleles(nullptr), gt(nullptr), gt_sum(nullptr), gt_sum_occ(nullptr),
		info(nullptr), fmt(nullptr)
	{}

	yon1_vnt_t(const yon1_vnt_t& other) :
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
		if(n_info) info = new containers::PrimitiveContainerInterface*[m_info];
		if(n_fmt)  fmt  = new containers::PrimitiveGroupContainerInterface*[m_fmt];
		for(int i = 0; i < n_alleles; ++i) alleles[i] = other.alleles[i];
		for(int i = 0; i < n_info; ++i) info[i] = other.info[i]->Clone();
		for(int i = 0; i < n_fmt; ++i)  fmt[i]  = other.fmt[i]->Clone();
		if(other.gt != nullptr) gt = new yon_gt(*other.gt);
		if(other.gt_sum != nullptr) gt_sum = new yon_gt_summary(*other.gt_sum);
		if(other.gt_sum_occ != nullptr){
			gt_sum_occ = new yon_gt_summary[gt->n_o];
			for(int i = 0; i < gt->n_o; ++i) gt_sum_occ[i] = other.gt_sum_occ[i];
		}
	}

	yon1_vnt_t& operator=(const yon1_vnt_t& other){
		// Delete existing data without consideration.
		delete [] alleles;
		delete gt, gt_sum;
		delete [] gt_sum_occ;
		for(int i = 0; i < n_info; ++i) delete info[i];
		for(int i = 0; i < n_fmt;  ++i) delete fmt[i];
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

		if(n_info) info = new containers::PrimitiveContainerInterface*[m_info];
		if(n_fmt)  fmt  = new containers::PrimitiveGroupContainerInterface*[m_fmt];
		for(int i = 0; i < n_alleles; ++i) alleles[i] = other.alleles[i];
		for(int i = 0; i < n_info; ++i) info[i] = other.info[i]->Clone();
		for(int i = 0; i < n_fmt; ++i)  fmt[i]  = other.fmt[i]->Clone();
		if(other.gt != nullptr) gt = new yon_gt(*other.gt);
		if(other.gt_sum != nullptr) gt_sum = new yon_gt_summary(*other.gt_sum);
		if(other.gt_sum_occ != nullptr){
			gt_sum_occ = new yon_gt_summary[gt->n_o];
			for(int i = 0; i < gt->n_o; ++i) gt_sum_occ[i] = other.gt_sum_occ[i];
		}

		return(*this);
	}

	yon1_vnt_t(yon1_vnt_t&& other) noexcept :
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
		if(n_info){
			info = new containers::PrimitiveContainerInterface*[n_info];
			for(int i = 0; i < n_info; ++i) info[i] = other.info[i]->Move();
		}

		if(n_fmt){
			fmt = new containers::PrimitiveGroupContainerInterface*[n_fmt];
			for(int i = 0; i < n_fmt; ++i) fmt[i] = other.fmt[i]->Move();
		}
	}
	//yon1_vnt_t& operator=(yon1_vnt_t&& other) noexcept = delete; // temporarily deleted

	~yon1_vnt_t(){
		delete [] alleles;
		delete gt, gt_sum;
		delete [] gt_sum_occ;
		for(int i = 0; i < n_info; ++i) delete info[i];
		for(int i = 0; i < n_fmt;  ++i) delete fmt[i];
		delete [] info; delete [] fmt;
	}

	bool EvaluateSummary(bool lazy_evaluate = true);
	bool EvaluateOcc(yon_occ& occ);
	bool EvaluateOccSummary(bool lazy_evaluate = true);

	bool UpdateBase(const bcf1_t* record){
		n_alleles = record->n_allele;
		qual = record->qual;
		rid = record->rid;
		pos = record->pos;
		name = std::string(record->d.id);
		delete [] alleles;
		alleles = new yon_allele[n_alleles];

		// Fix for the special case when ALT is not encoded
		if(this->n_alleles == 1){
			delete [] alleles;
			this->n_alleles = 2;
			this->alleles = new yon_allele[n_alleles];
			this->alleles[0] = record->d.allele[0];
			this->alleles[1].allele    = new char[1];
			this->alleles[1].allele[0] = '.';
			this->alleles[1].l_allele  = 1;
		} else {
			for(uint32_t i = 0; i < this->n_alleles; ++i){
				this->alleles[i] = record->d.allele[i];
			}
		}

		if(this->n_alleles == 2) this->controller.biallelic = true;
		this->controller.simple_snv = (this->alleles[0].length() == 1 && this->alleles[1].length() == 1);

		return true;
	}

	/**<
	 * Check if it is possible to pack the REF and ALT allele strings into
	 * a single uint8_t. The input allelic data has to be diploid, biallelic and
	 * match the regular expression pattern "^([ATGCN\\.]{1}){1}|(<NON_REF>){1}$"
	 * @return Returns TRUE if it is possible to pack the ref and alt alleles into a byte or FALSE otherwise.
	 */
	bool UsePackedRefAlt(void) const;

	/**<
	 * Bitpack biallelic, diploid REF and ALT data into a single uint8_t. Failure
	 * to check for validity beforehand with UsePackedRefAlt() may result in
	 * errors.
	 * @return Returns a bitpacked unsigned byte.
	 */
	uint8_t PackRefAltByte(void) const;

	/**<
	 * Returns the alleles as a concatenated string. For example the alleles
	 * A and T is returned as "A,T". This function is used primarily in the
	 * UpdateHtslibVcfRecord() function for exporting into a htslib bcf1_t record.
	 * @return Returns a concatenated string of alleles.
	 */
	std::string GetAlleleString(void) const{
		std::string ret = this->alleles[0].ToString();
		for(uint32_t i = 1; i < this->n_alleles; ++i)
			ret += "," + this->alleles[i].ToString();
		return(ret);
	}

	// Base setters for function pointer overloading.
	void SetController(const uint16_t value){ this->controller = value; }
	void SetBasePloidy(const uint8_t value){ this->n_base_ploidy = value; }
	void SetInfoPId(const int32_t value){ this->info_pid = value; }
	void SetFormatPId(const int32_t value){ this->fmt_pid = value; }
	void SetFilterPId(const int32_t value){ this->flt_pid = value; }
	void SetQuality(const float value){ this->qual = value; }
	void SetPosition(const int64_t value){ this->pos = value; }
	void SetChromosome(const uint32_t value){ this->rid = value; }
	void SetName(const std::string& value);

	/**<
	 * Updates a htslib bcf1_t record with data available in this meta record.
	 * This function is used when converting yon1_t records to bcf1_t records.
	 * @param rec Input bcf1_t record that has been allocated.
	 * @param hdr Input bcf hdr structure converted from tachyon header.
	 * @return Returns the input bcf1_t record pointer.
	 */
	bcf1_t* UpdateHtslibVcfRecord(bcf1_t* rec, bcf_hdr_t* hdr) const{
		rec->rid = this->rid;
		rec->pos = this->pos;
		bcf_update_id(hdr, rec, this->name.data());
		bcf_update_alleles_str(hdr, rec, this->GetAlleleString().data());
		if(std::isnan(this->qual)) bcf_float_set_missing(rec->qual);
		else rec->qual = this->qual;

		return(rec);
	}

	void OutputHtslibVcfInfo(bcf1_t* rec, bcf_hdr_t* hdr, DataBlockSettings& settings){
		for(uint32_t j = 0; j < n_info; ++j){
			if(info_hdr[j]->yon_type == YON_VCF_HEADER_FLAG){
				bcf_update_info_flag(hdr, rec, info_hdr[j]->id.data(), NULL, 1);
			} else {
				info[j]->UpdateHtslibVcfRecordInfo(rec, hdr, info_hdr[j]->id);
			}
		}
	}

	void OutputHtslibVcfFormat(bcf1_t* rec,
	                           bcf_hdr_t* hdr,
	                           DataBlockSettings& settings,
	                           yon_gt_rcd* external_exp) const
	{
		if(n_fmt){
			// Case when the only available FORMAT field is the GT field.
			if(n_fmt == 1 && is_loaded_gt &&
			   controller.gt_available &&
			   (settings.display_static & YON_BLK_BV_GT))
			{
				gt->ExpandExternal(external_exp);
				gt->d_exp = external_exp;
				gt->UpdateHtslibGenotypes(rec, hdr);
				gt->d_exp = nullptr;
			}
			// Case when there are > 1 Vcf Format fields and the GT field
			// is available.
			else if(n_fmt > 1 && is_loaded_gt &&
					controller.gt_available &&
					(settings.display_static & YON_BLK_BV_GT))
			{
				gt->ExpandExternal(external_exp);
				gt->d_exp = external_exp;

				gt->UpdateHtslibGenotypes(rec, hdr);
				for(uint32_t g = 1; g < n_fmt; ++g){
					if(fmt_hdr[g]->yon_type == YON_VCF_HEADER_FLOAT)
						fmt[g]->UpdateHtslibVcfRecordFormatFloat(rec, hdr, fmt_hdr[g]->id);
					else if(fmt_hdr[g]->yon_type == YON_VCF_HEADER_INTEGER)
						fmt[g]->UpdateHtslibVcfRecordFormatInt32(rec, hdr, fmt_hdr[g]->id);
					else if(fmt_hdr[g]->yon_type == YON_VCF_HEADER_STRING || fmt_hdr[g]->yon_type == YON_VCF_HEADER_CHARACTER)
						fmt[g]->UpdateHtslibVcfRecordFormatString(rec, hdr, fmt_hdr[g]->id);
				}
				gt->d_exp = nullptr;
			}
			// All other cases.
			else {
				for(uint32_t g = 0; g < n_fmt; ++g){
					if(fmt_hdr[g]->yon_type == YON_VCF_HEADER_FLOAT)
						fmt[g]->UpdateHtslibVcfRecordFormatFloat(rec, hdr, fmt_hdr[g]->id);
					else if(fmt_hdr[g]->yon_type == YON_VCF_HEADER_INTEGER)
						fmt[g]->UpdateHtslibVcfRecordFormatInt32(rec, hdr, fmt_hdr[g]->id);
					else if(fmt_hdr[g]->yon_type == YON_VCF_HEADER_STRING || fmt_hdr[g]->yon_type == YON_VCF_HEADER_CHARACTER)
						fmt[g]->UpdateHtslibVcfRecordFormatString(rec, hdr, fmt_hdr[g]->id);
				}
			}

		}
	}

	void OutputHtslibVcfFilter(bcf1_t* rec, bcf_hdr_t* hdr) const{
		if(n_flt){
			for(uint32_t k = 0; k < n_flt; ++k){
				int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, flt_hdr[k]->id.data());
				bcf_update_filter(hdr, rec, &tmpi, 1);
			}
		}
	}

	void ToVcfString(const VariantHeader& header,
	           io::BasicBuffer& buffer,
	           uint32_t display,
	           yon_gt_rcd* external_rcd = nullptr) const
	{
		buffer.Add(header.contigs_[rid].name.data(), header.contigs_[rid].name.size());
		buffer += '\t';
		buffer.AddReadble(pos+1);
		buffer += '\t';
		if(name.size()) buffer += name;
		else buffer += '.';
		buffer += '\t';
		buffer.Add(alleles[0].allele, alleles[0].l_allele);
		buffer += '\t';
		buffer.Add(alleles[1].allele, alleles[1].l_allele);
		for(int i = 2; i < n_alleles; ++i){
			buffer += ',';
			buffer.Add(alleles[i].allele, alleles[i].l_allele);
		}
		buffer += '\t';
		if(std::isnan(qual)) buffer += '.';
		else buffer.AddReadble(qual);
		buffer += '\t';

		if(flt_hdr.size() == 0) buffer += '.';
		else {
			buffer.Add(flt_hdr[0]->id.data(), flt_hdr[0]->id.size());
			for(int j = 1; j < flt_hdr.size(); ++j){
				buffer += ',';
				buffer.Add(flt_hdr[j]->id.data(), flt_hdr[j]->id.size());
			}
		}
		buffer += '\t';

		if(n_info){
			for(int j = 0; j < n_info; ++j){
				if(j != 0) buffer += ';';
				if(info_hdr[j]->yon_type == YON_VCF_HEADER_FLAG){
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

		if(n_fmt){
			buffer += '\t';
			buffer += fmt_hdr[0]->id;
			for(int j = 1; j < n_fmt; ++j){
				buffer += ':';
				buffer += fmt_hdr[j]->id;
			}
			buffer += '\t';

			if(n_fmt == 1 &&
			   controller.gt_available &&
			   (display & YON_BLK_BV_GT))
			{
				gt->ExpandExternal(external_rcd);
				gt->d_exp = external_rcd;

				// Iterate over samples and print FORMAT:GT value in Vcf format.
				gt->d_exp[0].PrintVcf(buffer, gt->m);
				for(uint32_t s = 1; s < header.GetNumberSamples(); ++s){
					buffer += '\t';
					gt->d_exp[s].PrintVcf(buffer, gt->m);
				}

				gt->d_exp = nullptr;
			}
			// Case when there are > 1 Vcf Format fields and the GT field
			// is available.
			else if(n_fmt > 1 && is_loaded_gt &&
					controller.gt_available &&
					(display & YON_BLK_BV_GT))
			{
				gt->ExpandExternal(external_rcd);
				gt->d_exp = external_rcd;

				gt->d_exp[0].PrintVcf(buffer, gt->m);
				for(uint32_t g = 1; g < n_fmt; ++g){
					buffer += ':';
					fmt[g]->ToVcfString(buffer, 0);
				}
				for(uint32_t s = 1; s < header.GetNumberSamples(); ++s){
					buffer += '\t';
					gt->d_exp[s].PrintVcf(buffer, gt->m);
					for(uint32_t g = 1; g < n_fmt; ++g){
						buffer += ':';
						fmt[g]->ToVcfString(buffer, s);
					}
				}

				gt->d_exp = nullptr;
			}
			// All other cases.
			else {
				fmt[0]->ToVcfString(buffer, 0);
				for(uint32_t g = 1; g < n_fmt; ++g){
					buffer += ':';
					fmt[g]->ToVcfString(buffer, 0);
				}

				for(uint32_t s = 1; s < header.GetNumberSamples(); ++s){
					buffer += '\t';
					fmt[0]->ToVcfString(buffer, s);
					for(uint32_t g = 1; g < n_fmt; ++g){
						buffer += ':';
						fmt[g]->ToVcfString(buffer, s);
					}
				}
			}
		}

		buffer += '\n';
	}

	bool AddInfoFlag(const std::string& tag, VariantHeader& header){
		const YonInfo* info_tag = header.GetInfo(tag);
		assert(info_tag != nullptr);

		int32_t offset = GetInfoOffset(tag);
		if(offset >= 0){
			delete info[offset];
			info[offset] = new containers::PrimitiveContainer<int8_t>();
		} else {
			info[n_info++] = new containers::PrimitiveContainer<int8_t>();
			info_hdr.push_back(info_tag);
		}
		return true;
	}

	template <class int_t>
	bool AddInfo(const std::string& tag, VariantHeader& header, int_t& data_point){
		const YonInfo* info_tag = header.GetInfo(tag);
		assert(info_tag != nullptr);

		int32_t offset = GetInfoOffset(tag);
		if(offset >= 0){
			delete info[offset];
			info[offset] = new containers::PrimitiveContainer<int_t>(data_point);
		} else {
			info[n_info++] = new containers::PrimitiveContainer<int_t>(data_point);
			info_hdr.push_back(info_tag);
		}
		return true;
	}

	template <class int_t>
	bool AddInfo(const std::string& tag, VariantHeader& header, int_t* data_point, const size_t n_points){
		const YonInfo* info_tag = header.GetInfo(tag);
		assert(info_tag != nullptr);

		int32_t offset = GetInfoOffset(tag);
		if(offset >= 0){
			delete info[offset];
			info[offset] = new containers::PrimitiveContainer<int_t>(data_point, n_points);
		} else {
			info[n_info++] = new containers::PrimitiveContainer<int_t>(data_point, n_points);
			info_hdr.push_back(info_tag);
		}
		return true;
	}

	template <class int_t>
	bool AddInfo(const std::string& tag, VariantHeader& header, containers::PrimitiveContainer<int_t>*& container){
		const YonInfo* info_tag = header.GetInfo(tag);
		assert(info_tag != nullptr);

		int32_t offset = GetInfoOffset(tag);
		if(offset >= 0){
			delete info[offset];
			info[offset] = container;
		} else {
			info[n_info++] = container;
			info_hdr.push_back(info_tag);
		}
		container = nullptr;
		return true;
	}

	bool AddFilter(const std::string& tag, VariantHeader& header){
		const io::VcfFilter* filter_tag = header.GetFilter(tag);
		assert(filter_tag != nullptr);

		int32_t offset = GetFilterOffset(tag);
		if(offset >= 0) return true;
		else {
			++n_flt;
			flt_hdr.push_back(filter_tag);
		}
		return true;
	}

	bool AddGenotypeStatistics(VariantHeader& header, const bool replace_existing = true){
		if(gt_sum == nullptr){
			std::cerr << utility::timestamp("LOG") << "Failed to update record. Genotype statistics have not been processed!" << std::endl;
			return false;
		}

		if(gt_sum->d == nullptr){
			std::cerr << utility::timestamp("LOG") << "Failed to update record. Genotype statistics have not been lazy evaluated!" << std::endl;
			return false;
		}

		if(n_info + 10 > m_info){
			m_info = std::max(n_info + 10, m_info + 10);
			containers::PrimitiveContainerInterface** old = info;
			info = new containers::PrimitiveContainerInterface*[m_info];
			for(int i = 0; i < n_info; ++i) info[i] = old[i];
			delete [] old;
		}

		this->AddInfo("NM", header, gt_sum->d->nm);
		this->AddInfo("NPM", header, gt_sum->d->npm);
		this->AddInfo("AN", header, gt_sum->d->an);
		this->AddInfo("HWE_P", header, gt_sum->d->hwe_p);

		if(gt_sum->d->n_ac_af > 2){
			this->AddInfo("AC", header, &gt_sum->d->ac[2], gt_sum->d->n_ac_af - 2);
			this->AddInfo("AF", header, &gt_sum->d->af[2], gt_sum->d->n_ac_af - 2);
		}

		if(this->n_alleles > 2) this->AddInfoFlag("MULTI_ALLELIC", header);

		// Special case for AC_P
		containers::PrimitiveContainer<uint32_t>* ac_p = new containers::PrimitiveContainer<uint32_t>();
		ac_p->resize(n_base_ploidy * (gt_sum->d->n_ac_af - 2));

		int j = 0;
		for(uint32_t p = 0; p < n_base_ploidy; ++p){
			for(uint32_t i = 2; i < gt_sum->d->n_ac_af; ++i, ++j){
				ac_p->at(j) = gt_sum->alleles_strand[p][i];
				++(*ac_p);
			}
		}
		this->AddInfo("AC_P", header, ac_p);

		if(n_base_ploidy == 2){
			this->AddInfo("F_PIC", header, gt_sum->d->f_pic);
			this->AddInfo("HET", header, gt_sum->d->heterozygosity);
		}

		assert(n_info == info_hdr.size());
		return true;
	}

	bool AddGenotypeStatisticsOcc(VariantHeader& header, std::vector<std::string>& names){
		if(n_info + names.size()*10 + 3 > m_info){
			m_info += names.size()*10 + 3;
			containers::PrimitiveContainerInterface** old = info;
			info = new containers::PrimitiveContainerInterface*[m_info];
			for(int i = 0; i < n_info; ++i) info[i] = old[i];
			delete [] old;
		}

		EvaluatePopGenOcc(header);

		for(int i = 0; i < names.size(); ++i){
			this->AddInfo(names[i]+"_NM", header, gt_sum_occ[i].d->nm);
			this->AddInfo(names[i]+"_NPM", header, gt_sum_occ[i].d->npm);
			this->AddInfo(names[i]+"_AN", header, gt_sum_occ[i].d->an);
			this->AddInfo(names[i]+"_HWE_P", header, gt_sum_occ[i].d->hwe_p);

			if(gt_sum_occ[i].d->n_ac_af > 2){
				this->AddInfo(names[i]+"_AC", header, &gt_sum_occ[i].d->ac[2], gt_sum_occ[i].d->n_ac_af - 2);
				this->AddInfo(names[i]+"_AF", header, &gt_sum_occ[i].d->af[2], gt_sum_occ[i].d->n_ac_af - 2);
			}

			if(this->n_alleles > 2) this->AddInfoFlag(names[i]+"_MULTI_ALLELIC", header);

			// Special case for AC_P
			containers::PrimitiveContainer<uint32_t>* ac_p = new containers::PrimitiveContainer<uint32_t>();
			ac_p->resize(n_base_ploidy * (gt_sum_occ[i].d->n_ac_af - 2));

			int j = 0;
			for(uint32_t p = 0; p < n_base_ploidy; ++p){
				for(uint32_t k = 2; k < gt_sum_occ[i].d->n_ac_af; ++k, ++j){
					ac_p->at(j) = gt_sum_occ->alleles_strand[p][k];
					++(*ac_p);
				}
			}
			this->AddInfo(names[i]+"_AC_P", header, ac_p);

			if(n_base_ploidy == 2){
				this->AddInfo(names[i]+"_F_PIC", header, gt_sum_occ[i].d->f_pic);
				this->AddInfo(names[i]+"_HET", header, gt_sum_occ[i].d->heterozygosity);
			}
		}

		assert(n_info == info_hdr.size());
		return true;
	}

	/**<
	 * Calculates various F-statistics over the provided groupings using
	 * the Occ table. Calculates one set of statistics over all the groupings
	 * and for each individual group. The calculations here are available only
	 * for diploid individuals at biallelic sites.
	 *
	 * Warning: this function returns incorrect results if the provided
	 *          grouping sets are overlapping. The incorrect results arise
	 *          from the fact that a single individual cannot simultaneously
	 *          exist in multiple distinct populations.
	 * @param header Src VariantHeader.
	 * @return       Returns TRUE upon success or FALSE otherwise.
	 */
	bool EvaluatePopGenOcc(VariantHeader& header){
		if(this->gt == nullptr) return false;
		if(this->gt->n_o == 0) return false;
		if(gt->n_allele != 2 || gt->m != 2) return false;

		// Number of genotypes across populations.
		const uint32_t n_gt_all =
				this->gt_sum->gt[2][2].n_cnt +
				this->gt_sum->gt[2][3].n_cnt +
				this->gt_sum->gt[3][2].n_cnt +
				this->gt_sum->gt[3][3].n_cnt;

		float* obs = new float[this->gt->n_o];
		float* exp = new float[this->gt->n_o];
		uint32_t* n_s = new uint32_t[this->gt->n_o];

		for(int i = 0; i < this->gt->n_o; ++i){
			const uint32_t het  = this->gt_sum_occ[i].gt[2][3].n_cnt + this->gt_sum_occ[i].gt[3][2].n_cnt;
			const uint32_t n_gt = het + this->gt_sum_occ[i].gt[2][2].n_cnt + this->gt_sum_occ[i].gt[3][3].n_cnt;

			const float p   = (2*(float)this->gt_sum_occ[i].gt[2][2].n_cnt + het) / (2*n_gt);
			const float q   = (2*(float)this->gt_sum_occ[i].gt[3][3].n_cnt + het) / (2*n_gt);;
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

		float exp_sum = 0, obs_sum = 0;
		for(int i = 0; i < this->gt->n_o; ++i){
			exp_sum += exp[i] * n_s[i];
			obs_sum += obs[i] * n_s[i];
		}

		const float h_i = obs_sum / n_gt_all;
		const float h_s = exp_sum / n_gt_all;
		// p-bar (frequency of allele A) over all populations
		// q-bar
		const float p_bar = (2*(float)this->gt_sum->gt[2][2].n_cnt + this->gt_sum->gt[2][3].n_cnt + this->gt_sum->gt[3][2].n_cnt) / (2*n_gt_all);
		const float q_bar = (2*(float)this->gt_sum->gt[3][3].n_cnt + this->gt_sum->gt[2][3].n_cnt + this->gt_sum->gt[3][2].n_cnt) / (2*n_gt_all);
		const float h_t   = 2*p_bar*q_bar;

		/*
		std::cerr << "pbar=" << p_bar << " qbar=" << q_bar << " with " << n_gt_all << std::endl;
		std::cerr << "sums: " << exp_sum << ", " << obs_sum << std::endl;
		std::cerr << "hi=" << h_i << ", hs=" << h_s << ", ht=" << h_t << std::endl;
		*/

		// Fstatistics
		// FIS [= (HS - HI)/HS]
		//std::cerr << "FIS=" << (h_s - h_i) / h_s << std::endl;
		// FST [= (HT - HS)/HT]
		//std::cerr << "FST=" << (h_t - h_s) / h_t << std::endl;
		// FIT [= (HT - HI)/HT]
		//std::cerr << "FIT=" << (h_t - h_i) / h_t << std::endl;
		float fis = (h_s - h_i) / h_s;
		float fst = (h_t - h_s) / h_t;
		float fit = (h_t - h_i) / h_t;

		this->AddInfo("FIS", header, fis);
		this->AddInfo("FST", header, fst);
		this->AddInfo("FIT", header, fit);



		delete [] exp; delete [] obs; delete [] n_s;

		return true;
	}

	containers::PrimitiveContainerInterface* GetInfo(const std::string& name) const {
		std::unordered_map<std::string, uint32_t>::const_iterator it = info_map.find(name);
		if(it != info_map.end()) return(info[it->second]);
		return(nullptr);
	}

	containers::PrimitiveGroupContainerInterface* GetFmt(const std::string& name) const {
		std::unordered_map<std::string, uint32_t>::const_iterator it = info_map.find(name);
		if(it != info_map.end()) return(fmt[it->second]);
		return(nullptr);
	}

	const io::VcfFilter* GetFlt(const std::string& name) const {
		std::unordered_map<std::string, uint32_t>::const_iterator it = flt_map.find(name);
		if(it != flt_map.end()) return(flt_hdr[it->second]);
		return(nullptr);
	}

	int32_t GetInfoOffset(const std::string& name) const {
		std::unordered_map<std::string, uint32_t>::const_iterator it = info_map.find(name);
		if(it != info_map.end()) return(it->second);
		return(-1);
	}

	int32_t GetFormatOffset(const std::string& name) const {
		std::unordered_map<std::string, uint32_t>::const_iterator it = info_map.find(name);
		if(it != info_map.end()) return(it->second);
		return(-1);
	}

	int32_t GetFilterOffset(const std::string& name) const {
		std::unordered_map<std::string, uint32_t>::const_iterator it = flt_map.find(name);
		if(it != flt_map.end()) return(it->second);
		return(-1);
	}

public:
	bool is_dirty; // if data has been modified
	bool is_loaded_gt;
	yon_vnt_cnt controller;
	uint8_t   n_base_ploidy;
	uint16_t  n_alleles;
	int32_t   n_flt, n_fmt, n_info;
	int32_t   info_pid, flt_pid, fmt_pid;
	int32_t   m_fmt, m_info, m_allele; // allocated sizes
	float     qual;
	uint32_t  rid;
	int64_t   pos;
	std::string name;

	std::vector<const YonInfo*> info_hdr;
	std::vector<const YonFormat*> fmt_hdr;
	std::vector<const io::VcfFilter*> flt_hdr;

	std::unordered_map<std::string, uint32_t> info_map;
	std::unordered_map<std::string, uint32_t> fmt_map;
	std::unordered_map<std::string, uint32_t> flt_map;

	yon_allele* alleles;
	yon_gt* gt;
	yon_gt_summary* gt_sum;
	yon_gt_summary* gt_sum_occ; // summary if occ has been invoked
	containers::PrimitiveContainerInterface** info;
	containers::PrimitiveGroupContainerInterface** fmt;
};

// Genotype helper
struct GenotypeSummary {
public:
	GenotypeSummary(void) :
		base_ploidy(0),
		phase_if_uniform(0),
		mixed_phasing(false),
		invariant(false),
		n_missing(0),
		n_vector_end(0)
	{}

	~GenotypeSummary() = default;

	/**<
	 * Gathers summary statistics for a vector of genotypes
	 * at a given site. Collects information regarding the
	 * number of missing genotypes and count of sentinel
	 * nodes, checks if the phasing is uniform and whether
	 * all the genotypes are identical.
	 * Todo: Use hashes to check for uniformity of genotypes.
	 * @param n_samples Total number of samples in the input vector.
	 *                  This is equivalent to the samples in the file.
	 * @param fmt       The target htslib format structure.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	template<class T>
	bool Evaluate(const size_t& n_samples, const bcf_fmt_t& fmt){
		if(fmt.p_len == 0) return true;
		assert(fmt.size/fmt.n == sizeof(T));

		// Set the base ploidy. This corresponds to the LARGEST
		// ploidy observed of ANY individual sample at the given
		// locus. If a genotype has a ploidy < base ploidy then
		// it is trailed with the sentinel node symbol to signal
		// that the remainder of the vector is NULL.
		this->base_ploidy = fmt.n;

		// Find first phase
		this->phase_if_uniform = 0;
		// Iterate over genotypes to find the first valid phase
		// continuing the search if the current value is the
		// sentinel node symbol.
		int j = fmt.n - 1;
		for(uint32_t i = 0; i < n_samples; ++i){
			if(io::VcfGenotype<T>::IsMissing(fmt.p[j]) == true
			   || io::VcfType<T>::IsVectorEnd(fmt.p[j]) == true)
				j += fmt.n;
			else {
				this->phase_if_uniform = fmt.p[j] & 1;
				break;
			}
		}

		// Iterate over genotypes to compute summary statistics
		// regarding missingness, number of special sentinel
		// symbols and assess uniformity of phasing.
		j = fmt.n - 1;
		for(uint32_t i = 0; i < n_samples; ++i){
			if(io::VcfGenotype<T>::IsMissing(fmt.p[j]) == false
			   && io::VcfType<int8_t>::IsVectorEnd(fmt.p[j]) == false
			   && (fmt.p[j] & 1) != this->phase_if_uniform)
			{
				this->mixed_phasing = true;
			}

			// Iterate over the number of chromosomes / individual
			for(int k = 0; k < fmt.n; ++k, ++j){
				this->n_missing    += io::VcfGenotype<T>::IsMissing(fmt.p[j]);
				this->n_vector_end += io::VcfType<T>::IsVectorEnd(fmt.p[j]);
			}
		}

		return true;
	}

	/**<
	 * Gathers summary statistics for a vector of genotypes
	 * at a given site. Collects information regarding the
	 * number of missing genotypes and count of sentinel
	 * nodes, checks if the phasing is uniform and whether
	 * all the genotypes are identical.
	 * Todo: Use hashes to check for uniformity of genotypes.
	 * @param rec       The src yon1_vnt_t record.
	 * @param n_samples Total number of samples in the input vector.
	 *                  This is equivalent to the samples in the file.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	bool Evaluate(const yon1_vnt_t& rec, const size_t& n_samples){
		if(rec.gt == nullptr)
			return false;

		if(rec.gt->rcds == nullptr)
			return false;

		// Set the base ploidy. This corresponds to the LARGEST
		// ploidy observed of ANY individual sample at the given
		// locus. If a genotype has a ploidy < base ploidy then
		// it is trailed with the sentinel node symbol to signal
		// that the remainder of the vector is NULL.
		this->base_ploidy = rec.gt->m;

		// Find first phase
		this->phase_if_uniform = 0;
		const uint8_t p_offset = std::max(0, (int)rec.gt->m - 1);
		// Iterate over genotypes to find the first valid phase
		// continuing the search if the current value is the
		// sentinel node symbol.
		for(uint32_t i = 0; i < rec.gt->n_i; ++i){
			if(YON_GT_RCD_ALLELE_UNPACK(rec.gt->rcds[i].allele[p_offset]) > 1){
				this->phase_if_uniform = (YON_GT_RCD_PHASE(rec.gt->rcds[i].allele[p_offset]));
				break;
			}
		}

		// Iterate over genotypes to compute summary statistics
		// regarding missingness, number of special sentinel
		for(uint32_t i = 0; i < rec.gt->n_i; ++i){
			if(YON_GT_RCD_ALLELE_UNPACK(rec.gt->rcds[i].allele[p_offset]) > 1
			   && (YON_GT_RCD_PHASE(rec.gt->rcds[i].allele[p_offset])) != this->phase_if_uniform)
			{
				this->mixed_phasing = true;
			}

			// Iterate over the number of chromosomes / individual
			for(int k = 0; k < rec.gt->m; ++k){
				this->n_missing    += YON_GT_RCD_ALLELE_UNPACK(rec.gt->rcds[i].allele[k]) == YON_GT_RCD_MISS ? rec.gt->rcds[i].run_length : 0;
				this->n_vector_end += YON_GT_RCD_ALLELE_UNPACK(rec.gt->rcds[i].allele[k]) == YON_GT_RCD_EOV ? rec.gt->rcds[i].run_length : 0;
			}
		}

		//std::cerr << utility::timestamp("DEBUG") << n_missing << "," << n_vector_end << "," << (int)base_ploidy << "," << phase_if_uniform << ":" << mixed_phasing << ":" << invariant << std::endl;

		return true;
	}

	inline bool isBaseDiploid(void) const{ return(this->base_ploidy == 2); }

public:
	uint8_t  base_ploidy;
	bool     phase_if_uniform;
	bool     mixed_phasing;
	bool     invariant;
	uint64_t n_missing;
	uint64_t n_vector_end;
};

struct yon_stats_sample {
public:
	yon_stats_sample(void);
	~yon_stats_sample(void);

	yon_stats_sample& operator+=(const yon_stats_sample& other);

	/**<
	 * Compute the number of transitions and transversions and then use
	 * those numbers to calculate the transition/transversion ratio. These
	 * values are stored internally in the struct.
	 * @return Returns TRUE.
	 */
	bool LazyEvalute(void);

	/**<
	 * Converts this struct into a partial (albeit valid) JSON string.
	 * Takes as argument an input buffer reference and a reference to
	 * the target sample name as described in the global header.
	 * @param buffer      Dst buffer.
	 * @param sample_name Src sample name.
	 * @return            Returns TRUE.
	 */
	io::BasicBuffer& ToJsonString(io::BasicBuffer& buffer, const std::string& sample_name) const;

	/**<
	 * Resets all internal counters. Useful when reusing
	 * objects without releasing memory.
	 */
	void reset(void);

public:
	uint64_t  n_ins, n_del, n_singleton;
	uint64_t  n_ts, n_tv;
	double    ts_tv_ratio;
	uint64_t* base_conv[9]; // {A,T,G,C,unknown,'.',EOV,ins,del}
	uint64_t* ins_del_dist; // indel distribution. values are zigzagged.
};

/**<
 * Supportive structure for computing variant-centric
 * summary statistics using the genotypes and site-specific
 * meta information.
 */
struct yon_stats_tstv {
public:
	// Supportive structure.
	struct yon_stats_tstv_obj {
		yon_stats_tstv_obj(): n_allele(0), n_non_ref(0), t_non_ref(0), allele_encodings(nullptr), non_ref_encodings(nullptr), b_size(nullptr){}
		yon_stats_tstv_obj(const uint32_t n_allele): n_allele(n_allele), n_non_ref(0), t_non_ref(0), allele_encodings(new uint8_t[n_allele]), non_ref_encodings(new uint8_t[n_allele]), b_size(new int32_t[n_allele]){ }
		~yon_stats_tstv_obj(){
			delete [] allele_encodings;
			delete [] non_ref_encodings;
			delete [] b_size;
		}

		uint32_t  n_allele;
		uint32_t  n_non_ref; // Number of alleles with non-ref alleles.
		uint32_t  t_non_ref; // Target id for non-ref allele.
		uint8_t*  allele_encodings;
		uint8_t*  non_ref_encodings;
		int32_t*  b_size;
	};

public:
	yon_stats_tstv();
	yon_stats_tstv(const uint32_t n_samples);
	~yon_stats_tstv(void);

	yon_stats_tstv& operator+=(const yon_stats_tstv& other);

	// Accessors
	inline yon_stats_sample& operator[](const uint32_t pos){ return(this->sample[pos]); }
	inline const yon_stats_sample& operator[](const uint32_t pos) const{ return(this->sample[pos]); }

	/**<
	 * Addition operator in the situation where the sample order
	 * in the other object is permuted according to a provided
	 * permutation array.
	 * @param other Src object to read data from.
	 * @param ppa   Src permutation array describing the relationship between the current object and the provided one
	 * @return      Returns an reference to the dst object.
	 */
	yon_stats_tstv& Add(const yon_stats_tstv& other, const yon_gt_ppa* ppa);

	/**<
	 * Allocates memory for a number of samples. All previous data
	 * will be deleted without consideration.
	 * @param n_samples
	 */
	void SetSize(const uint32_t n_samples);

	/**<
	 * Construct a valid JSON string from the internal structures. Requires
	 * as input a reference to the global header sample list.
	 * @param buffer       Dst data buffer.
	 * @param sample_names Src vector of sample names.
	 * @return             Returns a reference to the dst data buffer.
	 */
	io::BasicBuffer& ToJsonString(io::BasicBuffer& buffer, const std::vector<std::string>& sample_names) const;

	/**<
	 * Construct mappings for the alleles into relative offset. Function
	 * will allocate new memory for allele_encodings and non_ref_encodings.
	 * It is not legal to pass the same pointer to allele_encodings and
	 * non_ref_encodings.
	 * @param rcd               Src yon1_t structure that must have genotype data available.
	 * @param allele_encodings  Pointer to empty integer array.
	 * @param non_ref_encodings Pointer to empty integer array.
	 * @return                  Returns TRUE upon success or FALSE otherwise.
	 */
	bool GetEncodings(const yon1_vnt_t& rcd, yon_stats_tstv_obj& helper);

	/**<
	 * Update the current structure with the genotype data from the
	 * provided yon1_t record and a preallocated array of pointers
	 * to yon_gt_rcds. This function requires that the yon1_t record
	 * has genotype data available and that the yon_gt_rcds have been
	 * preprocessed. If genotype data not available or there and/or
	 * there are no meta data available then the function does not
	 * continue. This function is considerably slower than the
	 * Update(const yon1_t&) function that operates directly on the
	 * run-length encoded objects.
	 * @param rcd  Src yon1_t record.
	 * @param rcds Src yon_gt_rcd pointers.
	 * @return     Returns TRUE upon success or FALSE otherwise.
	 */
	bool Update(const yon1_vnt_t& rcd, yon_gt_rcd** rcds);

	/**<
	 * Update the current structure with the genotype data from the
	 * provided yon1_t record. This function requires that the yon1_t record
	 * has genotype data available. If genotype data not available or there and/or
	 * there are no meta data available then the function does not
	 * continue. This function does not guarantee that the ordering of the
	 * updates are correct in terms of global ordering. The updates will
	 * be completed according to the (possibly) local permutation array.
	 * Restoring global order is up to the user.
	 * @param rcd  Src yon1_t record.
	 * @param rcds Src yon_gt_rcd pointers.
	 * @return     Returns TRUE upon success or FALSE otherwise.
	 */
	bool Update(const yon1_vnt_t& rcd);
	void UpdateDiploid(const yon_gt* gt, yon_stats_tstv_obj& helper);
	void UpdateNPloidy(const yon_gt* gt, yon_stats_tstv_obj& helper);

	// Todo: expand if necessary.
	void Evaluate(void){
		for(uint32_t i = 0; i < this->n_s; ++i)
			this->sample[i].LazyEvalute();
	}

	/**<
	 * Resets all internal counters. Useful when reusing
	 * objects without releasing memory. Iteratively calls
	 * reset on the child objects.
	 */
	void reset(void);

public:
	uint32_t n_s;
	uint64_t n_rcds;
	uint64_t n_snp, n_mnp, n_ins, n_del, n_other; // determined during lazy evaluation
	uint64_t n_no_alts, n_biallelic, n_multi_allele, n_multi_allele_snp, n_singleton;
	yon_stats_sample* sample;
	uint64_t* alt_count; // number of alt distribution
};

}

#endif /* CORE_VARIANT_RECORD_H_ */
