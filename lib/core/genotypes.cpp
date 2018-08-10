#include "genotypes.h"

namespace tachyon{

yon_gt::~yon_gt(){
	delete [] d_bcf;
	delete [] d_bcf_ppa,
		delete [] rcds;
	delete [] d_exp;
	delete itree;
	delete gt_sum;
}

bool yon_gt::EvaluateSummary(bool lazy_evaluate){
	assert(this->rcds != nullptr);
	this->gt_sum = new yon_gt_summary(this->m, this->n_allele);
	*this->gt_sum += *this;
	if(lazy_evaluate) this->gt_sum->LazyEvaluate();
	return true;
}

}
