#include "variant_record.h"

namespace tachyon{

yon1_t::yon1_t(void) :
	is_dirty(false),
	is_loaded_meta(false),
	is_loaded_gt(false),
	n_format(0), n_info(0), n_filter(0),
	id_block(0),
	meta(nullptr),
	gt(nullptr),
	gt_sum(nullptr),
	occ(nullptr),
	info(nullptr),
	fmt(nullptr),
	info_containers(nullptr),
	format_containers(nullptr),
	gt_i(nullptr),
	info_ids(nullptr),
	format_ids(nullptr),
	filter_ids(nullptr),
	parent_container(nullptr)
{
}

yon1_t::~yon1_t(void){
	delete [] this->info;
	delete [] this->fmt;
	delete [] this->info_containers;
	delete [] this->format_containers;
	delete this->gt;
	delete this->gt_sum;
	// Do not delete occ. It is always borrowed!
}

bool yon1_t::EvaluateSummary(bool lazy_evaluate){
	assert(this->gt != nullptr);
	assert(this->gt->rcds != nullptr);

	if(this->gt_sum != nullptr)
		return true;

	this->gt_sum = new yon_gt_summary(this->gt->m, this->gt->n_allele);
	*this->gt_sum += *this->gt;
	if(lazy_evaluate) this->gt_sum->LazyEvaluate();
	return true;
}

bool yon1_t::EvaluateOcc(){
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

}
