#ifndef CORE_TS_TV_OBJECT_H_
#define CORE_TS_TV_OBJECT_H_

#include "meta_entry.h"
#include "variant_record.h"
#include "support/type_definitions.h"

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
	yon_stats_sample(void) :
		n_ins(0), n_del(0), n_singleton(0), n_ts(0), n_tv(0), ts_tv_ratio(0)
	{
		for(U32 i = 0; i < 9; ++i){
			this->base_conv[i] = new uint64_t[9];
			memset(&this->base_conv[i][0], 0, sizeof(uint64_t)*9);
		}
	}

	~yon_stats_sample(void){
		for(U32 i = 0; i < 9; ++i){
			delete [] this->base_conv[i];
		}
	}

	inline double GetTiTVRatio(void) const{
		// Prevent division by 0
		if(this->n_ts == 0) return 0;
		return((double)this->n_ts / this->n_tv);
	}

	bool LazyEvalute(void){
		// Transversions: A->C, C->A, T->G, G->T, A->T, T->A, C->G, G->C
		// Transitions:   A->G, G->A, C->T, T->C
		// Maps: A->0, T->1, G->2, C->3
		this->n_tv = this->base_conv[0][3] + this->base_conv[3][0] +
		             this->base_conv[1][2] + this->base_conv[2][1] +
		             this->base_conv[0][1] + this->base_conv[1][0] +
		             this->base_conv[3][2] + this->base_conv[2][3];
		this->n_ts = this->base_conv[0][2] + this->base_conv[2][0] +
		             this->base_conv[1][3] + this->base_conv[3][1];

		if(this->n_ts == 0) this->ts_tv_ratio = 0;
		else this->ts_tv_ratio = ((double)this->n_ts / this->n_tv);

		return true;
	}

	uint64_t  n_ins, n_del, n_singleton;
	uint64_t  n_ts, n_tv;
	double ts_tv_ratio;
	uint64_t* base_conv[9]; // {A,T,G,C,unknown,'.',EOV,ins,del}
};

struct yon_stats_tstv {
public:
	yon_stats_tstv() : n_s(0), n_rcds(0), n_snp(0), n_mnp(0), n_ins(0), n_del(0), n_other(0), n_no_alts(0), n_singleton(0),
	                   n_biallelic(0), n_multi_allele(0), n_multi_allele_snp(0), sample(nullptr){}
	yon_stats_tstv(const uint32_t n_samples) : n_s(n_samples), n_rcds(0), n_snp(0), n_mnp(0), n_ins(0), n_del(0), n_other(0), n_no_alts(0), n_singleton(0),
                                               n_biallelic(0), n_multi_allele(0), n_multi_allele_snp(0), sample(new yon_stats_sample[n_samples]){}
	~yon_stats_tstv(void){ delete [] this->sample; }

	inline yon_stats_sample& operator[](const uint32_t pos){ return(this->sample[pos]); }
	inline const yon_stats_sample& operator[](const uint32_t pos) const{ return(this->sample[pos]); }

	bool Update(const yon1_t& rcd, yon_gt_rcd** rcds) {
		if(rcd.is_loaded_gt == false || rcd.is_loaded_meta == false)
			return false;

		if((rcd.gt->eval_cont & YON_GT_UN_RCDS) == false){
			std::cerr << "eval because not done" << std::endl;
			bool eval = rcd.gt->Evaluate();
			if(eval == false){
				std::cerr << "failed to evaluate" << std::endl;
				return false;
			}
		}

		assert(this->n_s == rcd.gt->n_s);

		++this->n_rcds; // Update number of records.

		// Update count for target variant line type.
		if(rcd.meta->n_alleles > 2){
			bool is_snp = true;
			for(U32 i = 0; i < rcd.gt->n_allele; ++i){
				if(rcd.meta->alleles[i].l_allele != 1){
					is_snp = false;
					break;
				}
			}

			if(is_snp) ++this->n_multi_allele_snp;
			else ++this->n_multi_allele;
		} else if(rcd.meta->n_alleles == 2){
			++this->n_biallelic;
		}

		// For SNP->SNP,insertion
		if(rcd.meta->alleles[0].size() == 1){
			// Encode alleles.
			uint8_t* allele_encodings = new uint8_t[rcd.gt->n_allele + 2];
			uint8_t* non_ref_encodings = new uint8_t[rcd.gt->n_allele + 2];
			memset(non_ref_encodings, 1, sizeof(uint8_t)*(rcd.gt->n_allele + 2));
			allele_encodings[0] = 5; non_ref_encodings[0] = 0;
			allele_encodings[1] = 6; non_ref_encodings[1] = 0;
			non_ref_encodings[2] = 0;
			for(U32 i = 2; i < rcd.gt->n_allele + 2; ++i){
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

			uint32_t n_non_ref = 0;
			for(U32 i = 0; i < rcd.gt->n_s; ++i){
				for(U32 j = 0; j < rcd.gt->m; ++j){
					++this->sample[i].base_conv[allele_encodings[2]][allele_encodings[(rcds[i]->allele[j] >> 1)]];
					n_non_ref += non_ref_encodings[rcds[i]->allele[j] >> 1];
				}
			}

			//if(n_non_ref <= 1)
			//	std::cerr << rcd.meta->position+1 << ": " << n_non_ref << std::endl;

			if(n_non_ref == 0) ++this->n_no_alts;
			else if(n_non_ref == 1) ++this->n_singleton;

			delete [] allele_encodings;
			delete [] non_ref_encodings;

		} else {
			//std::cerr << "ref allele is not 1: " << rcd.meta->alleles[0].toString() << std::endl;
		}


		return true;
	}

	// todo update
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

		++this->n_rcds;
		if(rcd.meta->n_alleles > 2){
			bool is_snp = true;
			for(U32 i = 0; i < rcd.gt->n_allele; ++i){
				if(rcd.meta->alleles[i].l_allele != 1){
					is_snp = false;
					break;
				}
			}

			if(is_snp) ++this->n_multi_allele_snp;
			else ++this->n_multi_allele;
		} else if(rcd.meta->n_alleles == 2){
			++this->n_biallelic;
		} else if(rcd.meta->n_alleles == 1){
			++this->n_no_alts;
		}

		assert(this->n_s == rcd.gt->n_s);

		if(rcd.meta->alleles[0].size() == 1){
			// Encode alleles.
			uint8_t* allele_encodings = new uint8_t[rcd.gt->n_allele + 2];
			allele_encodings[0] = 5;
			allele_encodings[1] = 6;
			for(U32 i = 2; i < rcd.gt->n_allele + 2; ++i){
				if(rcd.meta->alleles[i - 2].l_allele == 1){
					allele_encodings[i] = YON_STATS_TSTV_LOOKUP[rcd.meta->alleles[i - 2].allele[0]];
				} else {
					if(rcd.meta->alleles[i - 2].l_allele > 1 &&
					   std::regex_match(rcd.meta->alleles[i - 2].toString(), constants::YON_REGEX_CANONICAL_BASES))
					{
						std::cerr << "is insertion" << std::endl;
					}
					allele_encodings[i] = 4;
				}
			}

			if(allele_encodings[2] > 3){
				std::cerr << "bad reference allele: " << rcd.meta->alleles[0].toString() << std::endl;
				return false;
			}

			uint32_t sample_offset = 0;
			for(U32 i = 0; i < rcd.gt->n_i; ++i){
				for(U32 r = 0; r < rcd.gt->rcds[i].run_length; ++r, ++sample_offset){
					for(U32 j = 0; j < rcd.gt->m; ++j){
						++this->sample[sample_offset].base_conv[allele_encodings[2]][allele_encodings[(rcd.gt->rcds[i].allele[j] >> 1)]];
					}
				}
			}
			assert(sample_offset == this->n_s);
		} else {
			//std::cerr << "ref allele is not 1: " << rcd.meta->alleles[0].toString() << std::endl;
		}
		return true;
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
