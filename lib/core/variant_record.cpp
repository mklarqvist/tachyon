#include "support/magic_constants.h"
#include "variant_record.h"

namespace tachyon{

bool yon1_vnt_t::EvaluateSummary(bool lazy_evaluate){
	assert(this->gt != nullptr);
	assert(this->gt->rcds != nullptr);

	if(this->gt_sum != nullptr)
		return true;

	this->gt_sum = new yon_gt_summary(this->gt->m, this->gt->n_allele);
	*this->gt_sum += *this->gt;
	if(lazy_evaluate) this->gt_sum->LazyEvaluate();
	return true;
}

bool yon1_vnt_t::EvaluateOccSummary(bool lazy_evaluate){
	assert(this->gt != nullptr);
	assert(this->gt->rcds != nullptr);

	if(this->gt_sum_occ != nullptr)
		return false;

	if(this->gt->n_o == 0)
		return false;

	this->gt_sum_occ = new yon_gt_summary[this->gt->n_o];
	for(int i = 0; i < this->gt->n_o; ++i){
		this->gt_sum_occ[i].Setup(this->gt->m, this->gt->n_allele);
		this->gt_sum_occ[i].Add(*this->gt, this->gt->n_i_occ[i], this->gt->d_occ[i]);
		if(lazy_evaluate) this->gt_sum_occ[i].LazyEvaluate();
	}

	return true;
}


bool yon1_vnt_t::EvaluateOcc(yon_occ& occ){
	if(occ.row_names.size() == 0)
		return false;

	if(this->gt == nullptr)
		return false;

	this->gt->n_o     = occ.occ.size();
	this->gt->n_i_occ = new uint32_t[this->gt->n_o];
	this->gt->d_occ   = new yon_gt_rcd*[this->gt->n_o];

	uint32_t* cum_sums = new uint32_t[this->gt->n_o]; // Total cumulative genotypes observed.
	uint32_t* cum_sums_hit = new uint32_t[this->gt->n_o]; // Number of non-zero runs observed.
	memset(cum_sums, 0, sizeof(uint32_t)*this->gt->n_o);
	memset(cum_sums_hit, 0, sizeof(uint32_t)*this->gt->n_o);

	for(uint32_t i = 0; i < this->gt->n_o; ++i){
		this->gt->d_occ[i]   = new yon_gt_rcd[std::min(this->gt->n_i, occ.cum_sums[i])];
		this->gt->n_i_occ[i] = 0;
	}

	// Iterate over available gt rcds.
	for(uint32_t j = 0; j < this->gt->n_i; ++j){
		// Iterate over available groupings in transposed occ.
		for(uint32_t i = 0; i < this->gt->n_o; ++i){
			const uint32_t to   = occ.vocc[cum_sums[i] + this->gt->rcds[j].run_length][i];
			const uint32_t from = occ.vocc[cum_sums[i]][i];
			if(to - from != 0){
				// Allocate memory for alleles.
				this->gt->d_occ[i][this->gt->n_i_occ[i]].allele = new uint8_t[this->gt->m];

				// Copy allelic data from recerence gt rcd.
				for(uint32_t k = 0; k < this->gt->m; ++k){
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
	for(int i = 0; i < this->gt->n_o; ++i){
		assert(cum_sums[i] == this->gt->n_s);
		assert(cum_sums_hit[i] == occ.cum_sums[i]);
	}

	delete [] cum_sums;
	delete [] cum_sums_hit;

	return(true);
}

bool yon1_vnt_t::UsePackedRefAlt(void) const{
	if(this->controller.biallelic == false || this->controller.diploid == false)
		return false;

	if(std::regex_match(std::string(this->alleles[0].allele, this->alleles[0].l_allele), constants::YON_REGEX_PACKED_ALLELES) &&
	   std::regex_match(std::string(this->alleles[1].allele, this->alleles[1].l_allele), constants::YON_REGEX_PACKED_ALLELES)){
		return true;
	}
	return false;
}

uint8_t yon1_vnt_t::PackRefAltByte(void) const{
	assert(this->UsePackedRefAlt());
	uint8_t ref_alt = 0; // start out with empty

	if(this->alleles[0].l_allele == 9 && strncmp(this->alleles[0].allele, "<NON_REF>", 9) == 0){
		ref_alt ^= YON_ALLELE_NON_REF << 4;
	} else {
		switch(this->alleles[0].allele[0]){
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

	if(this->alleles[1].l_allele == 9 && strncmp(this->alleles[1].allele, "<NON_REF>", 9) == 0){
		ref_alt ^= YON_ALLELE_NON_REF << 0;
	} else {
		switch(this->alleles[1].allele[0]){
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

}
