#ifndef CORE_TS_TV_OBJECT_H_
#define CORE_TS_TV_OBJECT_H_

#include "meta_entry.h"
#include "variant_record.h"

namespace tachyon{

// 0,1,2,3,4 -> A,T,G,C,X
const uint8_t YON_STATS_TSTV_LOOKUP[256] =
{4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,0,4,3,
 4,4,4,2,4,4,4,4,4,4,4,4,4,4,4,4,1,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4};

struct yon_stats_sample {
public:
	yon_stats_sample(void) :
		n_ins(0), n_del(0), n_singleton(0), n_ts(0), n_tv(0), ts_tv_ratio(0)
	{
		for(uint32_t i = 0; i < 9; ++i){
			this->base_conv[i] = new uint64_t[9];
			memset(&this->base_conv[i][0], 0, sizeof(uint64_t)*9);
		}
	}

	~yon_stats_sample(void){
		for(uint32_t i = 0; i < 9; ++i){
			delete [] this->base_conv[i];
		}
	}

	yon_stats_sample& operator+=(const yon_stats_sample& other){
		this->n_ins += other.n_ins;
		this->n_del += other.n_del;
		this->n_singleton += other.n_singleton;
		this->n_ts += other.n_ts;
		this->n_tv += other.n_tv;
		for(int i = 0; i < 9; ++i){
			for(uint32_t j = 0; j < 9; ++j){
				this->base_conv[i][j] += other.base_conv[i][j];
			}
		}

		return(*this);
	}

	inline double GetTiTVRatio(void) const{
		// Prevent division by 0
		if(this->n_ts == 0) return 0;
		return((double)this->n_ts / this->n_tv);
	}

	bool LazyEvalute(void){
		// Transversions: A->C, C->A, T->G, G->T, A->T, T->A, C->G, G->C
		// Transitions:   A->G, G->A, C->T, T->C
		// Maps:          A->0, T->1, G->2, C->3
		this->n_tv = this->base_conv[0][3] + this->base_conv[3][0] +
		             this->base_conv[1][2] + this->base_conv[2][1] +
		             this->base_conv[0][1] + this->base_conv[1][0] +
		             this->base_conv[3][2] + this->base_conv[2][3];
		this->n_ts = this->base_conv[0][2] + this->base_conv[2][0] +
		             this->base_conv[1][3] + this->base_conv[3][1];

		if(this->n_ts == 0) this->ts_tv_ratio = 0;
		else this->ts_tv_ratio = ((double)this->n_ts / this->n_tv);

		this->n_ins = this->base_conv[0][7] + this->base_conv[1][7] +
		              this->base_conv[2][7] + this->base_conv[3][7] +
		              this->base_conv[4][7];
		this->n_del = this->base_conv[4][8];

		return true;
	}

	io::BasicBuffer& ToJsonString(io::BasicBuffer& buffer, const std::string& sample_name) const {
		buffer +=  "\"" + sample_name + "\":{";
		buffer +=  "\"n_ins\":" + std::to_string(this->n_ins);
		buffer += ",\"n_del\":" + std::to_string(this->n_del);
		buffer += ",\"n_singleton\":" + std::to_string(this->n_singleton);
		buffer += ",\"n_ts\":"  + std::to_string(this->n_ts);
		buffer += ",\"n_tv\":"  + std::to_string(this->n_tv);
		buffer += ",\"ts_tv\":" + std::to_string(this->ts_tv_ratio);
		buffer += ",\"conv\":[";
		for(uint32_t i = 0; i < 9; ++i){
			if(i != 0) buffer += ',';
			buffer += '[';
			buffer.AddReadble((uint64_t)this->base_conv[i][0]);
			for(uint32_t j = 1; j < 9; ++j){
				buffer += ',';
				buffer.AddReadble((uint64_t)this->base_conv[i][j]);
			}
			buffer += ']';
		}
		buffer += ']';
		buffer += '}';
		return(buffer);
	}

	void reset(void){
		n_ins = 0, n_del = 0, n_singleton = 0;
		n_ts = 0, n_tv = 0;
		ts_tv_ratio = 0;
		for(int i = 0; i < 9; ++i){
			memset(&this->base_conv[i][0], 0, sizeof(uint64_t)*9);
		}
	}

public:
	uint64_t  n_ins, n_del, n_singleton;
	uint64_t  n_ts, n_tv;
	double    ts_tv_ratio;
	uint64_t* base_conv[9]; // {A,T,G,C,unknown,'.',EOV,ins,del}
};

struct yon_stats_tstv {
public:
	yon_stats_tstv() : n_s(0), n_rcds(0), n_snp(0), n_mnp(0), n_ins(0), n_del(0), n_other(0), n_no_alts(0), n_singleton(0),
	                   n_biallelic(0), n_multi_allele(0), n_multi_allele_snp(0), sample(nullptr){}
	yon_stats_tstv(const uint32_t n_samples) : n_s(n_samples), n_rcds(0), n_snp(0), n_mnp(0), n_ins(0), n_del(0), n_other(0), n_no_alts(0), n_singleton(0),
                                               n_biallelic(0), n_multi_allele(0), n_multi_allele_snp(0), sample(new yon_stats_sample[n_samples]){}
	~yon_stats_tstv(void){ delete [] this->sample; }

	yon_stats_tstv& operator+=(const yon_stats_tstv& other){
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

		for(int i = 0; i < this->n_s; ++i){
			this->sample[i] += other.sample[i];
		}

		return(*this);
	}

	yon_stats_tstv& Add(const yon_stats_tstv& other, const yon_gt_ppa& ppa){
		assert(ppa.n_s == this->n_s);
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

		for(int i = 0; i < this->n_s; ++i){
			this->sample[i] += other.sample[ppa[i]];
		}

		return(*this);
	}

	void SetSize(const uint32_t n_samples){
		delete [] sample;
		n_s = n_samples;
		sample = new yon_stats_sample[n_samples];
	}

	io::BasicBuffer& ToJsonString(io::BasicBuffer& buffer, const std::vector<std::string>& sample_names) const {
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
		buffer += "},\n";
		buffer += "\"PSI\":{\n";
		for(uint32_t i = 0; i < this->n_s; ++i){
			this->sample[i].LazyEvalute();
			if(i != 0) buffer += ",\n";
			this->sample[i].ToJsonString(buffer, sample_names[i]);
		}
		buffer += "\n}\n}\n";
		return(buffer);
	}

	inline yon_stats_sample& operator[](const uint32_t pos){ return(this->sample[pos]); }
	inline const yon_stats_sample& operator[](const uint32_t pos) const{ return(this->sample[pos]); }

	bool GetEncodings(const yon1_t& rcd,
	                  __restrict__ uint8_t*& allele_encodings,
	                  __restrict__ uint8_t*& non_ref_encodings)
	{
		// Update number of records.
		++this->n_rcds;

		// Update count for target variant line type.
		// Case: variant site is multi-allelic.
		if(rcd.meta->n_alleles > 2){
			bool is_snp = true;
			for(uint32_t i = 0; i < rcd.gt->n_allele; ++i){
				if(rcd.meta->alleles[i].l_allele != 1){
					is_snp = false;
					break;
				}
			}

			if(is_snp) ++this->n_multi_allele_snp;
			else ++this->n_multi_allele;
		}
		// Case: variant site is biallelic.
		else if(rcd.meta->n_alleles == 2){
			++this->n_biallelic;
		}

		// For SNV to SNV or insertion. It is not possible to have a deletion
		// if the reference value is represented as a SNV.
		allele_encodings  = new uint8_t[rcd.gt->n_allele + 2];
		non_ref_encodings = new uint8_t[rcd.gt->n_allele + 2];

		// If the reference allele is of a single character wide then
		// we assume the reference site is a single base. This implicitly
		// means there cannot be any encodings for deletions.
		if(rcd.meta->alleles[0].size() == 1){
			// Encode alleles.
			memset(non_ref_encodings, 1, sizeof(uint8_t)*(rcd.gt->n_allele + 2));
			allele_encodings[0] = 5; non_ref_encodings[0] = 0;
			allele_encodings[1] = 6; non_ref_encodings[1] = 0;
			non_ref_encodings[2] = 0;
			for(uint32_t i = 2; i < rcd.gt->n_allele + 2; ++i){
				if(rcd.meta->alleles[i - 2].l_allele == 1){
					allele_encodings[i] = YON_STATS_TSTV_LOOKUP[rcd.meta->alleles[i - 2].allele[0]];
				} else {
					if(rcd.meta->alleles[i - 2].l_allele > 1 &&
					   std::regex_match(rcd.meta->alleles[i - 2].toString(), constants::YON_REGEX_CANONICAL_BASES))
					{
						//std::cerr << "is insertion: " << rcd.meta->alleles[i - 2].toString() << std::endl;
						allele_encodings[i] = 7;
					} else allele_encodings[i] = 4;
				}
			}

			if(allele_encodings[2] > 3){
				std::cerr << "bad reference allele: " << rcd.meta->alleles[0].toString() << std::endl;
				delete [] allele_encodings;
				delete [] non_ref_encodings;
				return false;
			}

			//if(n_non_ref <= 1)
			//	std::cerr << rcd.meta->position+1 << ": " << n_non_ref << std::endl;
		}
		// For insertion/deletion to SNV/insertion/deletion.
		else {
			//std::cerr << "ref allele is not 1: " << rcd.meta->alleles[0].toString() << std::endl;

			// Encode alleles.
			memset(non_ref_encodings, 1, sizeof(uint8_t)*(rcd.gt->n_allele + 2));
			allele_encodings[0] = 5;
			allele_encodings[1] = 6;
			memset(non_ref_encodings, 0, sizeof(uint8_t)*3);

			const uint16_t& ref_length = rcd.meta->alleles[0].l_allele;
			allele_encodings[2] = 4;

			// Iterate over available alleles.
			for(uint32_t i = 3; i < rcd.gt->n_allele + 2; ++i){
				// If the target allele is a simple SNV.
				if(rcd.meta->alleles[i - 2].l_allele == 1){
					//std::cerr << "target is deletion: " << i - 2 << "/" << rcd.meta->n_alleles << "; " << rcd.meta->alleles[i - 2].toString() << std::endl;
					allele_encodings[i] = 8;
				}
				// If the target allele length is shorter than the reference
				// allele length and is comprised of only canonical bases then
				// classify this allele as a deletion.
				else if(rcd.meta->alleles[i - 2].l_allele < ref_length &&
						std::regex_match(rcd.meta->alleles[i - 2].toString(), constants::YON_REGEX_CANONICAL_BASES))
				{
					//std::cerr << "is canonical deletion " << i - 2 << "/" << rcd.meta->n_alleles << std::endl;
					allele_encodings[i] = 8;
				} else {
					if(rcd.meta->alleles[i - 2].l_allele > ref_length &&
					   std::regex_match(rcd.meta->alleles[i - 2].toString(), constants::YON_REGEX_CANONICAL_BASES))
					{
						//std::cerr << "is insertion: " << rcd.meta->alleles[i - 2].toString() << ": " << i - 2 << "/" << rcd.meta->n_alleles << std::endl;
						allele_encodings[i] = 7;
					} else {
						//std::cerr << "is same: " << rcd.meta->alleles[i - 2].toString() << ": " << i - 2 << "/" << rcd.meta->n_alleles << std::endl;
						allele_encodings[i] = 4;
					}
				}
			}
		}

		return true;
	}

	bool Update(const yon1_t& rcd, yon_gt_rcd** rcds) {
		if(rcd.is_loaded_gt == false || rcd.is_loaded_meta == false)
			return false;

		if((rcd.gt->eval_cont & YON_GT_UN_RCDS) == false){
			bool eval = rcd.gt->Evaluate();
			if(eval == false){
				std::cerr << "failed to evaluate" << std::endl;
				return false;
			}
		}

		assert(this->n_s == rcd.gt->n_s);
		uint8_t *allele_encodings, *non_ref_encodings;
		if(this->GetEncodings(rcd, allele_encodings, non_ref_encodings) == false)
			return false;

		// Perform actual work.
		uint32_t n_non_ref = 0;
		uint32_t t_non_ref = 0;
		for(uint32_t i = 0; i < rcd.gt->n_s; ++i){
			for(uint32_t j = 0; j < rcd.gt->m; ++j){
				++this->sample[i].base_conv[allele_encodings[2]][allele_encodings[(rcds[i]->allele[j] >> 1)]];
				n_non_ref += non_ref_encodings[rcds[i]->allele[j] >> 1];
				t_non_ref += non_ref_encodings[rcds[i]->allele[j] >> 1] * i;
			}
		}

		if(n_non_ref == 0) ++this->n_no_alts;
		else if(n_non_ref == 1){
			++this->n_singleton;
			++this->sample[t_non_ref].n_singleton;
			assert(t_non_ref < this->n_s);
		}

		// Cleanup.
		delete [] allele_encodings;
		delete [] non_ref_encodings;

		return true;
	}

	bool Update(const yon1_t& rcd){
		if(rcd.is_loaded_gt == false || rcd.is_loaded_meta == false)
			return false;

		if((rcd.gt->eval_cont & YON_GT_UN_RCDS) == false){
			bool eval = rcd.gt->Evaluate();
			if(eval == false){
				std::cerr << "failed to evaluate" << std::endl;
				return false;
			}
		}

		uint8_t *allele_encodings, *non_ref_encodings;
		if(this->GetEncodings(rcd, allele_encodings, non_ref_encodings) == false)
			return false;

		uint32_t n_non_ref = 0;
		uint32_t t_non_ref = 0;
		if(rcd.gt->m == 2) this->UpdateDiploid(rcd.gt, allele_encodings, non_ref_encodings, n_non_ref, t_non_ref);
		else this->UpdateNPloidy(rcd.gt, allele_encodings, non_ref_encodings, n_non_ref, t_non_ref);

		if(n_non_ref == 0) ++this->n_no_alts;
		else if(n_non_ref == 1){
			++this->n_singleton;
			++this->sample[t_non_ref].n_singleton;
			assert(t_non_ref < this->n_s);
		}

		// Cleanup.
		delete [] allele_encodings;
		delete [] non_ref_encodings;

		return true;
	}

	void UpdateDiploid(const yon_gt* gt,
	                   __restrict__ const uint8_t* allele_encodings,
	                   __restrict__ const uint8_t* non_ref_encodings,
	                   uint32_t& n_non_ref,
	                   uint32_t& t_non_ref)
	{
		uint32_t sample_offset = 0;
		for(uint32_t i = 0; i < gt->n_i; ++i){
			if(((gt->rcds[i].allele[0] >> 1) == YON_GT_RCD_REF) && ((gt->rcds[i].allele[1] >> 1) == YON_GT_RCD_REF)){
				sample_offset += gt->rcds[i].run_length;
				continue;
			}

			for(uint32_t r = 0; r < gt->rcds[i].run_length; ++r, ++sample_offset){
				for(uint32_t j = 0; j < gt->m; ++j){
					++this->sample[sample_offset].base_conv[allele_encodings[2]][allele_encodings[(gt->rcds[i].allele[j] >> 1)]];
					n_non_ref += non_ref_encodings[gt->rcds[i].allele[j] >> 1];
					t_non_ref += non_ref_encodings[gt->rcds[i].allele[j] >> 1] * i;
				}
			}
		}
		assert(sample_offset == this->n_s);
	}

	void UpdateNPloidy(const yon_gt* gt,
	                   __restrict__ const uint8_t* allele_encodings,
	                   __restrict__ const uint8_t* non_ref_encodings,
	                   uint32_t& n_non_ref,
	                   uint32_t& t_non_ref)
	{
		uint32_t sample_offset = 0;
		for(uint32_t i = 0; i < gt->n_i; ++i){
			for(uint32_t r = 0; r < gt->rcds[i].run_length; ++r, ++sample_offset){
				for(uint32_t j = 0; j < gt->m; ++j){
					++this->sample[sample_offset].base_conv[allele_encodings[2]][allele_encodings[(gt->rcds[i].allele[j] >> 1)]];
					n_non_ref += non_ref_encodings[gt->rcds[i].allele[j] >> 1];
					t_non_ref += non_ref_encodings[gt->rcds[i].allele[j] >> 1] * i;
				}
			}
		}
		assert(sample_offset == this->n_s);
	}

	void reset(void){
		n_rcds = 0;
		n_snp = 0;
		n_snp = 0, n_mnp = 0, n_ins = 0, n_del = 0, n_other = 0;
		n_no_alts = 0, n_biallelic = 0, n_multi_allele = 0, n_multi_allele_snp = 0, n_singleton = 0;
		for(int i = 0; i < this->n_s; ++i)
			this->sample[i].reset();
	}

public:
	uint32_t n_s;
	uint64_t n_rcds;
	uint64_t n_snp, n_mnp, n_ins, n_del, n_other; // determined during lazy evaluation
	uint64_t n_no_alts, n_biallelic, n_multi_allele, n_multi_allele_snp, n_singleton;
	yon_stats_sample* sample;
};

}



#endif /* CORE_TS_TV_OBJECT_H_ */
