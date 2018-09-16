#include "support/magic_constants.h"
#include "ts_tv_object.h"

namespace tachyon {

yon_stats_sample::yon_stats_sample(void) :
	n_ins(0), n_del(0), n_singleton(0), n_ts(0), n_tv(0), ts_tv_ratio(0), ins_del_dist(new uint64_t[512])
{
	for(uint32_t i = 0; i < 9; ++i){
		this->base_conv[i] = new uint64_t[9];
		memset(&this->base_conv[i][0], 0, sizeof(uint64_t)*9);
	}
	memset(ins_del_dist,0,512*sizeof(uint64_t));
}

yon_stats_sample::~yon_stats_sample(void){
	for(uint32_t i = 0; i < 9; ++i){
		delete [] this->base_conv[i];
	}
	delete [] ins_del_dist;
}

yon_stats_sample& yon_stats_sample::operator+=(const yon_stats_sample& other){
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

	for(int i = 0; i < 512; ++i)
		this->ins_del_dist[i] += other.ins_del_dist[i];

	return(*this);
}

bool yon_stats_sample::LazyEvalute(void){
	// Transversions: A->C, C->A, T->G, G->T, A->T, T->A, C->G, G->C
	// Transitions:   A->G, G->A, C->T, T->C
	this->n_tv = this->base_conv[YON_GT_TSTV_A][YON_GT_TSTV_C] + this->base_conv[YON_GT_TSTV_C][YON_GT_TSTV_A] +
				 this->base_conv[YON_GT_TSTV_T][YON_GT_TSTV_G] + this->base_conv[YON_GT_TSTV_G][YON_GT_TSTV_T] +
				 this->base_conv[YON_GT_TSTV_A][YON_GT_TSTV_T] + this->base_conv[YON_GT_TSTV_T][YON_GT_TSTV_A] +
				 this->base_conv[YON_GT_TSTV_C][YON_GT_TSTV_G] + this->base_conv[YON_GT_TSTV_G][YON_GT_TSTV_C];
	this->n_ts = this->base_conv[YON_GT_TSTV_A][YON_GT_TSTV_G] + this->base_conv[YON_GT_TSTV_G][YON_GT_TSTV_A] +
				 this->base_conv[YON_GT_TSTV_T][YON_GT_TSTV_C] + this->base_conv[YON_GT_TSTV_C][YON_GT_TSTV_T];

	if(this->n_tv == 0) this->ts_tv_ratio = 0;
	else this->ts_tv_ratio = ((double)this->n_ts / this->n_tv);

	this->n_ins = this->base_conv[YON_GT_TSTV_A][YON_GT_TSTV_INS] + this->base_conv[YON_GT_TSTV_T][YON_GT_TSTV_INS] +
				  this->base_conv[YON_GT_TSTV_C][YON_GT_TSTV_INS] + this->base_conv[YON_GT_TSTV_C][YON_GT_TSTV_INS] +
				  this->base_conv[YON_GT_TSTV_UNKNOWN][YON_GT_TSTV_INS];
	this->n_del = this->base_conv[YON_GT_TSTV_UNKNOWN][YON_GT_TSTV_DEL];

	return true;
}

io::BasicBuffer& yon_stats_sample::ToJsonString(io::BasicBuffer& buffer, const std::string& sample_name) const {
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

	buffer += ",\"in_dist\":[";
	//buffer.AddReadble(this->ins_del_dist[0]);
	for(int i = 1; i < 512; i+=2){
		if(i != 1) buffer += ',';
		buffer.AddReadble(this->ins_del_dist[i]);
	}
	buffer += "],";
	buffer += ",\"del_dist\":[";
	//buffer.AddReadble(this->ins_del_dist[0]);
	for(int i = 2; i < 512; i+=2){
		if(i != 2) buffer += ',';
		buffer.AddReadble(this->ins_del_dist[i]);
	}
	buffer += ']';

	buffer += ']';
	buffer += '}';
	return(buffer);
}

void yon_stats_sample::reset(void){
	n_ins = 0, n_del = 0, n_singleton = 0;
	n_ts = 0, n_tv = 0;
	ts_tv_ratio = 0;

	for(int i = 0; i < 9; ++i){
		memset(&this->base_conv[i][0], 0, sizeof(uint64_t)*9);
	}

	for(int i = 0; i < 512; ++i)
		this->ins_del_dist[i] = 0;
}


yon_stats_tstv::yon_stats_tstv() : n_s(0), n_rcds(0), n_snp(0), n_mnp(0), n_ins(0), n_del(0), n_other(0), n_no_alts(0), n_singleton(0),
				   n_biallelic(0), n_multi_allele(0), n_multi_allele_snp(0), sample(nullptr), alt_count(new uint64_t[32]){ memset(alt_count, 0, sizeof(uint64_t)*32); }
yon_stats_tstv::yon_stats_tstv(const uint32_t n_samples) : n_s(n_samples), n_rcds(0), n_snp(0), n_mnp(0), n_ins(0), n_del(0), n_other(0), n_no_alts(0), n_singleton(0),
										   n_biallelic(0), n_multi_allele(0), n_multi_allele_snp(0), sample(new yon_stats_sample[n_samples]), alt_count(new uint64_t[32]){ memset(alt_count, 0, sizeof(uint64_t)*32); }
yon_stats_tstv::~yon_stats_tstv(void){ delete [] this->sample; delete [] this->alt_count; }

yon_stats_tstv& yon_stats_tstv::operator+=(const yon_stats_tstv& other){
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
	for(int i = 0; i < 32; ++i) this->alt_count[i] = other.alt_count[i];

	return(*this);
}

yon_stats_tstv& yon_stats_tstv::Add(const yon_stats_tstv& other, const yon_gt_ppa& ppa){
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
		this->sample[ppa[i]] += other.sample[i];
	}
	for(int i = 0; i < 32; ++i) this->alt_count[i] = other.alt_count[i];

	return(*this);
}

void yon_stats_tstv::SetSize(const uint32_t n_samples){
	delete [] sample;
	n_s = n_samples;
	sample = new yon_stats_sample[n_samples];
}

io::BasicBuffer& yon_stats_tstv::ToJsonString(io::BasicBuffer& buffer, const std::vector<std::string>& sample_names) const {
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
	for(int i = 0; i < 32; ++i){
		if(i != 0) buffer += ',';
		buffer.AddReadble(this->alt_count[i]);
	}
	buffer += "]},\n";
	buffer += "\"PSI\":{\n";
	for(uint32_t i = 0; i < this->n_s; ++i){
		if(i != 0) buffer += ",\n";
		this->sample[i].ToJsonString(buffer, sample_names[i]);
	}
	buffer += "\n}\n}\n";
	return(buffer);
}

bool yon_stats_tstv::GetEncodings(const yon1_vnt_t& rcd, yon_stats_tstv_obj& helper)
{
	// Update number of records.
	++this->n_rcds;

	if(rcd.n_alleles == 1){
		std::cerr << utility::timestamp("LOG") << "Cannot have site with no ALT alleles described..." << std::endl;
		return false;
	}

	// Update count for target variant line type.
	// Case: variant site is multi-allelic.
	if(rcd.n_alleles > 2){
		bool is_snp = true;
		for(uint32_t i = 0; i < rcd.gt->n_allele; ++i){
			if(rcd.alleles[i].l_allele != 1){
				is_snp = false;
				break;
			}
		}

		if(is_snp) ++this->n_multi_allele_snp;
		else ++this->n_multi_allele;
	}
	// Case: variant site is biallelic.
	else if(rcd.n_alleles == 2){
		++this->n_biallelic;
	}

	// For SNV to SNV or insertion. It is not possible to have a deletion
	// if the reference value is represented as a SNV.
	// If the reference allele is of a single character wide then
	// we assume the reference site is a single base. This implicitly
	// means there cannot be any encodings for deletions.
	if(rcd.alleles[0].size() == 1){
		// Encode alleles.
		memset(helper.non_ref_encodings, 1, sizeof(uint8_t)*(rcd.gt->n_allele + 2));
		helper.allele_encodings[0]  = YON_GT_TSTV_MISS; helper.non_ref_encodings[0] = 0;
		helper.allele_encodings[1]  = YON_GT_TSTV_EOV;  helper.non_ref_encodings[1] = 0;
		helper.non_ref_encodings[2] = 0;
		memset(helper.b_size, 0, sizeof(int32_t)*2);

		// Iterate over available alleles.
		for(uint32_t i = 2; i < rcd.gt->n_allele + 2; ++i){
			if(rcd.alleles[i - 2].l_allele == 1){
				helper.allele_encodings[i] = YON_STATS_TSTV_LOOKUP[rcd.alleles[i - 2].allele[0]];
				helper.b_size[i] = 0;
			} else {
			if(rcd.alleles[i - 2].l_allele > 1 &&
				   std::regex_match(rcd.alleles[i - 2].ToString(), constants::YON_REGEX_CANONICAL_BASES))
				{
					//std::cerr << "is insertion: " << rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele << std::endl;
					//std::cerr << "is insertion: " << rcd.alleles[i - 2].toString() << std::endl;
					helper.allele_encodings[i] = YON_GT_TSTV_INS;
					helper.b_size[i] = std::min(255, rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele);
				} else {
					helper.allele_encodings[i] = YON_GT_TSTV_UNKNOWN;
					helper.b_size[i] = 0;
				}
			}
		}

		// Ascertain that the reference allele is valid for this purpose.
		if(helper.allele_encodings[YON_GT_RCD_REF] > 3){
			std::cerr << utility::timestamp("LOG") << "Bad reference allele: " << rcd.alleles[0].ToString() << std::endl;
			return false;
		}
	}
	// For insertion/deletion to SNV/insertion/deletion.
	else {
		// Encode alleles.
		memset(helper.non_ref_encodings, 1, sizeof(uint8_t)*(rcd.gt->n_allele + 2));
		helper.allele_encodings[0] = YON_GT_TSTV_MISS;
		helper.allele_encodings[1] = YON_GT_TSTV_EOV;
		memset(helper.non_ref_encodings, 0, sizeof(uint8_t)*3);

		const uint16_t& ref_length = rcd.alleles[0].l_allele;
		helper.allele_encodings[2] = YON_GT_TSTV_UNKNOWN;

		memset(helper.b_size, 0, sizeof(int32_t)*3);

		// Iterate over available alleles.
		for(uint32_t i = 3; i < rcd.gt->n_allele + 2; ++i){
			// If the target allele is a simple SNV.
			if(rcd.alleles[i - 2].l_allele == 1){
				if(std::regex_match(rcd.alleles[i - 2].ToString(), constants::YON_REGEX_CANONICAL_BASES)){
					//std::cerr << "target is deletion: " << i - 2 << "/" << rcd.n_alleles << "; " << rcd.alleles[0].ToString() << "->" << rcd.alleles[i - 2].ToString() << " size: " << rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele << std::endl;
					helper.allele_encodings[i] = YON_GT_TSTV_DEL;
					helper.b_size[i] = std::max(-255, (int)rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele);
				} else {
					helper.allele_encodings[i] = YON_GT_TSTV_UNKNOWN;
					helper.b_size[i] = 0;
				}
			}
			// If the target allele length is shorter than the reference
			// allele length and is comprised of only canonical bases then
			// classify this allele as a deletion.
			else if(rcd.alleles[i - 2].l_allele < ref_length &&
					std::regex_match(rcd.alleles[i - 2].ToString(), constants::YON_REGEX_CANONICAL_BASES))
			{
				//std::cerr << "target is deletion: " << i - 2 << "/" << rcd.n_alleles << "; " << rcd.alleles[0].ToString() << "->" << rcd.alleles[i - 2].ToString() << " size: " << rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele << std::endl;
				helper.allele_encodings[i] = YON_GT_TSTV_DEL;
				helper.b_size[i] = std::max(-255, (int)rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele);
			} else {
				if(rcd.alleles[i - 2].l_allele > ref_length &&
				   std::regex_match(rcd.alleles[i - 2].ToString(), constants::YON_REGEX_CANONICAL_BASES))
				{
					//std::cerr << "is insertion: " << rcd.alleles[i - 2].ToString() << ": " << rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele << "/" << rcd.n_alleles << std::endl;
					helper.allele_encodings[i] = YON_GT_TSTV_INS;
					helper.b_size[i] = std::min(255, rcd.alleles[i - 2].l_allele - rcd.alleles[0].l_allele);
				} else {
					//std::cerr << "is same: " << rcd.alleles[i - 2].toString() << ": " << i - 2 << "/" << rcd.n_alleles << std::endl;
					helper.allele_encodings[i] = YON_GT_TSTV_UNKNOWN;
					helper.b_size[i] = 0;
				}
			}
		}
	}

	return true;
}

bool yon_stats_tstv::Update(const yon1_vnt_t& rcd, yon_gt_rcd** rcds) {
	if(rcd.is_loaded_gt == false)
		return false;

	if((rcd.gt->eval_cont & YON_GT_UN_RCDS) == false){
		bool eval = rcd.gt->Evaluate();
		if(eval == false){
			std::cerr << "failed to evaluate" << std::endl;
			return false;
		}
	}

	assert(this->n_s == rcd.gt->n_s);
	yon_stats_tstv_obj helper(rcd.gt->n_allele + 2);
	if(this->GetEncodings(rcd, helper) == false)
		return false;

	// Update number of alt alleles.
	++this->alt_count[std::min(31,rcd.n_alleles - 1)];

	// Perform actual work.
	uint32_t n_non_ref = 0;
	uint32_t t_non_ref = 0;
	for(uint32_t i = 0; i < rcd.gt->n_s; ++i){
		for(uint32_t j = 0; j < rcd.gt->m; ++j){
			const uint8_t allele = YON_GT_ALLELE_UNPACK(rcds[i]->allele[j]);
			++this->sample[i].base_conv[helper.allele_encodings[YON_GT_RCD_REF]][helper.allele_encodings[allele]];
			++this->sample[i].ins_del_dist[(helper.b_size[allele] << 1) ^ (helper.b_size[allele] >> 31)];
			n_non_ref += helper.non_ref_encodings[allele];
			t_non_ref += helper.non_ref_encodings[allele] * i;
		}
	}

	if(n_non_ref == 0) ++this->n_no_alts;
	else if(n_non_ref == 1){
		++this->n_singleton;
		++this->sample[t_non_ref].n_singleton;
		assert(t_non_ref < this->n_s);
	}

	return true;
}

bool yon_stats_tstv::Update(const yon1_vnt_t& rcd){
	if(rcd.is_loaded_gt == false)
		return false;

	if((rcd.gt->eval_cont & YON_GT_UN_RCDS) == false){
		bool eval = rcd.gt->Evaluate();
		if(eval == false){
			std::cerr << "failed to evaluate" << std::endl;
			return false;
		}
	}

	yon_stats_tstv_obj helper(rcd.gt->n_allele + 2);
	if(this->GetEncodings(rcd, helper) == false)
		return false;

	// Update number of alt alleles.
	++this->alt_count[std::min(31, rcd.n_alleles - 1)];

	if(rcd.gt->m == 2) this->UpdateDiploid(rcd.gt, helper);
	else this->UpdateNPloidy(rcd.gt, helper);

	if(helper.n_non_ref == 0) ++this->n_no_alts;
	else if(helper.n_non_ref == 1){
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
	for(uint32_t i = 0; i < gt->n_i; ++i){
		if((YON_GT_ALLELE_UNPACK(gt->rcds[i].allele[0]) == YON_GT_RCD_REF) && (YON_GT_ALLELE_UNPACK(gt->rcds[i].allele[1]) == YON_GT_RCD_REF)){
			sample_offset += gt->rcds[i].run_length;
			continue;
		}

		for(uint32_t r = 0; r < gt->rcds[i].run_length; ++r, ++sample_offset){
			for(uint32_t j = 0; j < gt->m; ++j){
				const uint8_t allele = YON_GT_ALLELE_UNPACK(gt->rcds[i].allele[j]);
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
	for(uint32_t i = 0; i < gt->n_i; ++i){
		// If current run-length encoded object has all reference
		// template then continue.
		uint32_t n_refs = 0;
		for(uint32_t j = 0; j < gt->m; ++j)
			n_refs += (YON_GT_ALLELE_UNPACK(gt->rcds[i].allele[j]) == YON_GT_RCD_REF);

		if(n_refs == gt->m){
			sample_offset += gt->rcds[i].run_length;
			continue;
		}

		// Iterate over samples in the current run-length encoded object.
		for(uint32_t r = 0; r < gt->rcds[i].run_length; ++r, ++sample_offset){
			for(uint32_t j = 0; j < gt->m; ++j){
				const uint8_t allele = YON_GT_ALLELE_UNPACK(gt->rcds[i].allele[j]);
				++this->sample[sample_offset].base_conv[helper.allele_encodings[YON_GT_RCD_REF]][helper.allele_encodings[allele]];
				++this->sample[sample_offset].ins_del_dist[(helper.b_size[allele] << 1) ^ (helper.b_size[allele] >> 31)];
				helper.n_non_ref += helper.non_ref_encodings[allele];
				helper.t_non_ref += helper.non_ref_encodings[allele] * sample_offset;
			}
		}
	}
	assert(sample_offset == this->n_s);
}

void yon_stats_tstv::reset(void){
	n_rcds = 0;
	n_snp = 0, n_mnp = 0, n_ins = 0, n_del = 0, n_other = 0;
	n_no_alts = 0, n_biallelic = 0, n_multi_allele = 0, n_multi_allele_snp = 0, n_singleton = 0;
	for(int i = 0; i < this->n_s; ++i) this->sample[i].reset();
	for(int i = 0; i < 32; ++i) this->alt_count[i] = 0;
}

}
