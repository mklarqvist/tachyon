#ifndef CORE_VARIANT_RECORD_H_
#define CORE_VARIANT_RECORD_H_

#include "containers/components/variant_block_footer.h"
#include "containers/genotype_container_interface.h"
#include "containers/primitive_container.h"
#include "containers/primitive_group_container.h"
#include "occ.h"
#include "genotypes.h"

namespace tachyon{

const char* const REFALT_LOOKUP = "ATGC.XN";
#define YON_ALLELE_A       0
#define YON_ALLELE_T       1
#define YON_ALLELE_G       2
#define YON_ALLELE_C       3
#define YON_ALLELE_MISS    4
#define YON_ALLELE_NON_REF 5
#define YON_ALLELE_N       6


struct yon_allele {
public:
	yon_allele() : l_allele(0), allele(nullptr){}
	yon_allele(const std::string& string) : l_allele(string.size()), allele(new char[l_allele]){ memcpy(allele, string.data(), l_allele); }
	yon_allele& operator=(const std::string& value){ delete [] allele; l_allele = value.size(); allele = new char[l_allele]; memcpy(allele, value.data(), l_allele); return(*this); }
	yon_allele& operator=(const char* value){ delete [] allele; l_allele = strlen(value); allele = new char[l_allele]; memcpy(allele, value, l_allele); return(*this); }
	yon_allele& operator=(const char& value){ delete [] allele; l_allele = 1; allele = new char[l_allele]; allele[0] = value; return(*this); }

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
 * Primary lazy-evaluated record of a variant. Construction is done
 * outside of this definition. The struct is constructed from pointers
 * to the appropriate location in other lazy-evaluated objects such as
 * Info, Format, Genotype, and Meta containers. Evaluation into this object
 * is relatively inexpensive in isolation but quiet expensive when considered
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
		alleles(nullptr), gt(nullptr), gt_sum(nullptr), occ(nullptr),
		info(nullptr), fmt(nullptr) //gt_raw(nullptr),
		//info_ids(nullptr), fmt_ids(nullptr), flt_ids(nullptr)
	{}

	~yon1_vnt_t(){
		delete [] alleles;
		delete gt, gt_sum, occ;
		for(int i = 0; i < n_info; ++i) delete info[i];
		for(int i = 0; i < n_fmt;  ++i) delete fmt[i];
		delete [] info; delete [] fmt; //delete gt_raw;
		//delete info_ids, fmt_ids, flt_ids;
	}

	bool EvaluateSummary(bool lazy_evaluate = true);
	bool EvaluateOcc();

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
	void SetPosition(const uint64_t value){ this->pos = value; }
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
		if(n_info){
			for(uint32_t j = 0; j < n_info; ++j){
				if(info_hdr[j]->yon_type == YON_VCF_HEADER_FLAG){
					bcf_update_info_flag(hdr, rec, info_hdr[j]->id.data(), NULL, 1);
				} else {
					info[j]->UpdateHtslibVcfRecordInfo(rec, hdr, info_hdr[j]->id);
				}
			}

			if(settings.annotate_extra){
				this->EvaluateSummary(true);
				gt_sum->d->UpdateHtslibVcfRecord(rec, hdr);
			}
		} else {
			if(settings.annotate_extra){
				this->EvaluateSummary(true);
				gt_sum->d->UpdateHtslibVcfRecord(rec, hdr);
			}
		}
	}

	void OutputHtslibVcfFormat(bcf1_t* rec, bcf_hdr_t* hdr, DataBlockSettings& settings, yon_gt_rcd** external_exp) const{
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

	void Print(const VariantHeader& header,
	           io::BasicBuffer& buffer,
	           uint32_t display,
	           yon_gt_rcd** external_rcd = nullptr) const
	{
		buffer.Add(header.contigs_[rid].name.data(), header.contigs_[rid].name.size());
		buffer += '\t';
		buffer.AddReadble(pos+1);
		buffer += '\t';
		buffer += name;
		buffer += '\t';
		buffer.Add(alleles[0].allele, alleles[0].l_allele);
		buffer += '\t';
		buffer.Add(alleles[1].allele, alleles[1].l_allele);
		for(int i = 2; i < n_alleles; ++i){
			buffer += ',';
			buffer.Add(alleles[i].allele, alleles[i].l_allele);
		}
		buffer += '\t';
		buffer.AddReadble(qual);
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

		buffer += '\t';
		if(n_fmt){
			buffer += fmt_hdr[0]->id;
			for(int j = 1; j < n_fmt; ++j){
				buffer += ':';
				buffer += fmt_hdr[j]->id;
			}
			buffer += '\t';

			if(n_fmt == 1 && is_loaded_gt &&
			   controller.gt_available &&
			   (display & YON_BLK_BV_GT))
			{
				gt->ExpandExternal(external_rcd);
				gt->d_exp = external_rcd;

				// Iterate over samples and print FORMAT:GT value in Vcf format.
				gt->d_exp[0]->PrintVcf(buffer, gt->m);
				for(uint32_t s = 1; s < header.GetNumberSamples(); ++s){
					buffer += '\t';
					gt->d_exp[s]->PrintVcf(buffer, gt->m);
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

				gt->d_exp[0]->PrintVcf(buffer, gt->m);
				for(uint32_t g = 1; g < n_fmt; ++g){
					buffer += ':';
					fmt[g]->ToVcfString(buffer, 0);
				}
				for(uint32_t s = 1; s < header.GetNumberSamples(); ++s){
					buffer += '\t';
					gt->d_exp[s]->PrintVcf(buffer, gt->m);
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
		} else {
			buffer += '.';
		}

		buffer += '\n';
	}

	bool UpdateYonRecord(VariantHeader& header){
		std::cerr << "NM=" << header.info_fields_[header.info_fields_map_["NM"]].id << std::endl;
		info[n_info++] = new containers::PrimitiveContainer<uint32_t>(&gt_sum->d->nm, 1);
		info_hdr.push_back(&header.info_fields_[header.info_fields_map_["NM"]]);
		std::cerr << "NPM=" << header.info_fields_[header.info_fields_map_["NPM"]].id << std::endl;
		info[n_info++] = new containers::PrimitiveContainer<uint32_t>(&gt_sum->d->npm, 1);
		info_hdr.push_back(&header.info_fields_[header.info_fields_map_["NPM"]]);
		std::cerr << "AN=" << header.info_fields_[header.info_fields_map_["AN"]].id << std::endl;
		info[n_info++] = new containers::PrimitiveContainer<uint32_t>(&gt_sum->d->an, 1);
		info_hdr.push_back(&header.info_fields_[header.info_fields_map_["AN"]]);
		std::cerr << "HWE_P=" << header.info_fields_[header.info_fields_map_["HWE_P"]].id << std::endl;
		info[n_info++] = new containers::PrimitiveContainer<float>(&gt_sum->d->hwe_p, 1);
		info_hdr.push_back(&header.info_fields_[header.info_fields_map_["HWE_P"]]);
		if(gt_sum->d->n_ac_af > 2){
			std::cerr << "AC=" << header.info_fields_[header.info_fields_map_["AC"]].id << std::endl;
			info[n_info++] = new containers::PrimitiveContainer<uint32_t>(gt_sum->d->ac, gt_sum->d->n_ac_af - 2);
			info_hdr.push_back(&header.info_fields_[header.info_fields_map_["AC"]]);
			std::cerr << "AF=" << header.info_fields_[header.info_fields_map_["AF"]].id << std::endl;
			info[n_info++] = new containers::PrimitiveContainer<float>(gt_sum->d->af, gt_sum->d->n_ac_af - 2);
			info_hdr.push_back(&header.info_fields_[header.info_fields_map_["AF"]]);
		}
		return true;
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
	uint64_t  pos;
	std::string name;

	std::vector<const YonInfo*> info_hdr;
	std::vector<const YonFormat*> fmt_hdr;
	std::vector<const io::VcfFilter*> flt_hdr;

	yon_allele* alleles;
	yon_gt* gt;
	yon_gt_summary* gt_sum;
	yon_occ* occ;
	containers::PrimitiveContainerInterface** info;
	containers::PrimitiveGroupContainerInterface** fmt;
	//containers::GenotypeContainerInterface* gt_raw;

	//std::vector<int>* info_ids;
	//std::vector<int>* fmt_ids;
	//std::vector<int>* flt_ids;
};

}

#endif /* CORE_VARIANT_RECORD_H_ */
