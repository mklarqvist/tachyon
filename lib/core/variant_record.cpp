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

bool yon1_vnt_t::EvaluateOcc(){
	assert(this->gt != nullptr);
	assert(occ != nullptr);

	this->gt->n_o   = occ->occ.size();
	this->gt->n_occ = new uint32_t[this->gt->n_o];
	this->gt->d_occ = new yon_gt_rcd*[this->gt->n_o];

	for(uint32_t i = 0; i < this->gt->n_o; ++i){
		this->gt->d_occ[i] = new yon_gt_rcd[this->gt->n_i];

		uint32_t cum_sum  = 0; // Total cumulative genotypes observed.
		uint32_t cum_sum_hit = 0; // Number of non-zero runs observed.
		uint32_t n_offset = 0; // Virtual offset in destination array.

		// Iterate over available gt rcds.
		for(uint32_t j = 0; j < this->gt->n_i; ++j){
			const uint32_t to   = this->occ->occ[i][cum_sum + this->gt->rcds[j].run_length];
			const uint32_t from = this->occ->occ[i][cum_sum];
			if(to - from != 0){
				// Allocate memory for alleles.
				this->gt->d_occ[i][n_offset].allele = new uint8_t[this->gt->m];

				// Copy allelic data from recerence gt rcd.
				for(uint32_t k = 0; k < this->gt->m; ++k){
					this->gt->d_occ[i][n_offset].allele[k] = this->gt->rcds[j].allele[k];
				}

				// Set run-length representation.
				this->gt->d_occ[i][n_offset].run_length = to - from;
				assert(n_offset < this->gt->n_i);
				++n_offset;
				cum_sum_hit += to - from;
			}
			cum_sum += this->gt->rcds[j].run_length;
		}
		assert(cum_sum == this->gt->n_s);
		assert(cum_sum_hit == this->occ->cum_sums[i]);
		this->gt->n_occ[i] = n_offset;
	}
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
