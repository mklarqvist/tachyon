#include "genotypes.h"

namespace tachyon{

yon_gt::~yon_gt(){
	delete [] d_bcf;
	delete [] d_bcf_ppa,
		delete [] rcds;
	delete [] d_exp;
	delete itree;
}

bool yon_gt_summary::AddGenotypeLayer(yon_gt_summary_obj* target, uint8_t depth){
	if(depth == this->n_ploidy) return false;
	if(target->children == nullptr){
		target->children = new yon_gt_summary_obj[this->n_alleles];
	}

	for(U32 i = 0; i < this->n_alleles; ++i)
		this->AddGenotypeLayer(&target->children[i], depth+1);

	return false;
}

yon_gt_summary& yon_gt_summary::operator+=(const yon_gt& gt){
	assert(gt.rcds != nullptr);

	// Iterate over available genotype records.
	for(U32 i = 0; i < gt.n_i; ++i){
		// Target root node.
		yon_gt_summary_obj* target = &this->gt[gt.rcds[i].allele[0] >> 1];
		target->n_cnt += gt.rcds[i].run_length;

		// Iterate over alleles given the base ploidy.
		for(U32 j = 0; j < gt.m; ++j){
			assert((gt.rcds[i].allele[j] >> 1) < this->n_alleles);
			// Add allelic counts.
			this->alleles[gt.rcds[i].allele[j] >> 1] += gt.rcds[i].run_length;
			// Add strand-specific (ploidy-aware) alleleic counts.
			this->alleles_strand[j][gt.rcds[i].allele[j] >> 1] += gt.rcds[i].run_length;
		}

		// Update remainder.
		for(U32 j = 1; j < gt.m; ++j){
			target = &target->children[gt.rcds[i].allele[j] >> 1];
			target->n_cnt += gt.rcds[i].run_length;
		}
	}
	return(*this);
}

std::vector<uint64_t> yon_gt_summary::GetAlleleCounts(void) const{
	std::vector<uint64_t> c_allele(this->n_alleles, 0);
	for(U32 i = 0; i < this->n_alleles; ++i)
		c_allele[i] = this->alleles[i];

	return(c_allele);
}

std::vector< std::pair<uint64_t,double> > yon_gt_summary::GetAlleleCountFrequency(void) const{
	std::vector< std::pair<uint64_t,double> > c_allele(this->n_alleles);
	uint64_t n_total = 0;
	for(U32 i = 0; i < 2; ++i){
		c_allele[i].first = this->alleles[i];
	}

	for(U32 i = 2; i < this->n_alleles; ++i){
		c_allele[i].first = this->alleles[i];
		n_total += this->alleles[i];
	}

	if(n_total != 0){
		for(U32 i = 0; i < this->n_alleles; ++i){
			c_allele[i].second = (double)c_allele[i].first / n_total;
		}
	}

	return(c_allele);
}

std::vector< std::vector<uint64_t> > yon_gt_summary::GetAlleleStrandCounts(void) const{
	std::vector< std::vector<uint64_t> > c_allele(this->n_ploidy, std::vector<uint64_t>(this->n_alleles));
	for(U32 i = 0; i < this->n_ploidy; ++i){
		for(U32 j = 0; j < this->n_alleles; ++j){
			c_allele[i][j] = this->alleles_strand[i][j];
		}
	}

	return(c_allele);
}

bool yon_gt_summary::GetGenotype(std::vector<uint64_t>& data,
                                 yon_gt_summary_obj* target,
                                 uint8_t depth) const
{
	if(depth + 1 == this->n_ploidy){
		assert(target->children != nullptr);
		for(U32 i = 0; i < this->n_alleles; ++i){
			if(target->children[i].n_cnt != 0){
				data.push_back(target->children[i].n_cnt);
			}
		}
		return false;
	}

	for(U32 i = 0; i < this->n_alleles; ++i){
		this->GetGenotype(data, &target->children[i], depth + 1);
	}

	return(true);
}

// Todo: unfinished.
std::vector<yon_gt_rcd> yon_gt_summary::GetGenotypeCounts(bool drop_empty) const{
	std::vector<yon_gt_rcd> genotypes;

	// Traverse the trie and store records when hitting the leafs.
	std::vector<uint64_t> d;

	// If the target ploidy is greater than one (haploid) we collect
	// the genotypes by traversing the trie. Otherwise the genotypes
	// are the allele counts at the roots.
	if(this->n_ploidy > 1){
		for(U32 i = 0; i < this->n_alleles; ++i){
			//std::cerr << "outer " << i << "/" << (int)this->n_alleles << std::endl;
			this->GetGenotype(d, &this->gt[i], 1);
		}
	} else {
		for(U32 i = 0; i < this->n_alleles; ++i){
			if(this->gt[i].n_cnt != 0){
				d.push_back(this->gt[i].n_cnt);
			}
		}
	}
	uint64_t n_total = 0;
	for(U32 i = 0; i < d.size(); ++i)
		n_total += d[i];
	std::cerr << "collected: " << d.size() << " total = " << n_total << std::endl;
	return(genotypes);
}

std::vector<double> yon_gt_summary::GetStrandBiasAlleles(const bool phred_scale) const{
	if(this->n_ploidy != 2 || this->n_alleles + 2 < 2)
		return std::vector<double>();

	std::vector<double> strand_bias_p_values;
	double fisher_left_p, fisher_right_p, fisher_twosided_p;

	uint64_t n_cnt_fwd = 0;
	uint64_t n_cnt_rev = 0;
	for(U32 i = 2; i < this->n_alleles; ++i){
		n_cnt_fwd += this->alleles_strand[0][i];
		n_cnt_rev += this->alleles_strand[1][i];
	}

	kt_fisher_exact(
	this->alleles_strand[0][2], // A: Allele on forward strand
	this->alleles_strand[1][2], // B: Allele on reverse strand
	n_cnt_fwd - this->alleles_strand[0][2], // C: Not allele on forward strand
	n_cnt_rev - this->alleles_strand[1][2], // D: Not allele on reverse strand
	&fisher_left_p, &fisher_right_p, &fisher_twosided_p);

	if(phred_scale) strand_bias_p_values.push_back(std::abs(-10 * log10(fisher_twosided_p)));
	else strand_bias_p_values.push_back(fisher_twosided_p);

	// If n_alleles = 2 then they are identical because of symmetry
	if(this->n_alleles - 2 > 2){
		for(U32 p = 3; p < this->n_alleles; ++p){
			kt_fisher_exact(
			this->alleles_strand[0][p], // A: Allele on forward strand
			this->alleles_strand[1][p], // B: Allele on reverse strand
			n_cnt_fwd - this->alleles_strand[0][p], // C: Not allele on forward strand
			n_cnt_rev - this->alleles_strand[1][p], // D: Not allele on reverse strand
			&fisher_left_p, &fisher_right_p, &fisher_twosided_p);

			if(phred_scale) strand_bias_p_values.push_back(std::abs(-10 * log10(fisher_twosided_p)));
			else strand_bias_p_values.push_back(fisher_twosided_p);
		}
	}
	return(strand_bias_p_values);
}

double yon_gt_summary::CalculateHardyWeinberg(void) const{
	if(this->n_ploidy != 2 || this->n_alleles - 2 != 2) return -1;

	U64 obs_hets = this->gt[2][3].n_cnt + this->gt[3][2].n_cnt; // alts
	U64 obs_hom1 = this->gt[2][2].n_cnt; // hom ref
	U64 obs_hom2 = this->gt[3][3].n_cnt; // hom alt

	U64 obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	U64 obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

	int64_t rare_copies = 2 * obs_homr + obs_hets;
	int64_t genotypes   = obs_hets + obs_homc + obs_homr;

	double* het_probs = new double[rare_copies + 1];

	int64_t i;
	for (i = 0; i <= rare_copies; ++i)
		het_probs[i] = 0.0;

	/* start at midpoint */
	int64_t mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

	/* check to ensure that midpoint and rare alleles have same parity */
	if ((rare_copies & 1) ^ (mid & 1))
		++mid;

	int64_t curr_hets = mid;
	int64_t curr_homr = (rare_copies - mid) / 2;
	int64_t curr_homc = genotypes - curr_hets - curr_homr;

	het_probs[mid] = 1.0;
	double sum = het_probs[mid];
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2){
		het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
						   / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
		sum += het_probs[curr_hets - 2];

		/* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
		++curr_homr;
		++curr_homc;
	}

	curr_hets = mid;
	curr_homr = (rare_copies - mid) / 2;
	curr_homc = genotypes - curr_hets - curr_homr;
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2){
		het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
						/((curr_hets + 2.0) * (curr_hets + 1.0));
		sum += het_probs[curr_hets + 2];

		/* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
		--curr_homr;
		--curr_homc;
	}

	for (i = 0; i <= rare_copies; i++)
		het_probs[i] /= sum;

	double p_hwe = 0.0;
	/*  p-value calculation for p_hwe  */
	for (i = 0; i <= rare_copies; i++){
		if (het_probs[i] > het_probs[obs_hets])
			continue;

		p_hwe += het_probs[i];
	}

	p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

	delete [] het_probs;

	return(p_hwe);
}

bool yon_gt_summary::LazyEvaluate(void){
	delete this->d;
	this->d = new yon_gt_summary_rcd;

	// Allele count and frequency
	this->d->n_ploidy = this->n_ploidy;
	this->d->n_ac_af  = this->n_alleles;
	this->d->ac       = new uint64_t[this->n_alleles];
	this->d->af       = new double[this->n_alleles];

	uint64_t n_total = 0;
	for(U32 i = 0; i < 2; ++i) this->d->ac[i] = this->alleles[i];
	for(U32 i = 2; i < this->n_alleles; ++i){
		this->d->ac[i] = this->alleles[i];
		n_total += this->alleles[i];
	}

	if(n_total != 0){
		for(U32 i = 0; i < this->n_alleles; ++i)
			this->d->af[i] = (double)this->d->ac[i] / n_total;

	} else {
		for(U32 i = 0; i < this->n_alleles; ++i)
			this->d->af[i] = 0;
	}

	this->d->npm = this->d->ac[0];
	this->d->nm  = n_total;
	this->d->an  = n_total + this->d->ac[0];
	this->d->hwe_p = this->CalculateHardyWeinberg();
	this->d->ac_p = this->alleles_strand;

	// Strand-specific bias and inbreeding coefficient (F-statistic)
	if(this->n_ploidy == 2){
		this->d->heterozygosity = ((double)this->gt[2][3].n_cnt + this->gt[3][2].n_cnt) /
			(this->gt[2][2].n_cnt + this->gt[2][3].n_cnt + this->gt[3][2].n_cnt + this->gt[3][3].n_cnt);

		uint8_t n_fs_used = (this->n_alleles - 2 == 2 ? 1 : this->n_alleles - 2);
		this->d->n_fs = n_fs_used;
		this->d->fs_a = new double[n_fs_used];
		double fisher_left_p, fisher_right_p, fisher_twosided_p;
		uint64_t n_cnt_fwd = 0;
		uint64_t n_cnt_rev = 0;
		for(U32 i = 2; i < this->n_alleles; ++i){
			n_cnt_fwd += this->alleles_strand[0][i];
			n_cnt_rev += this->alleles_strand[1][i];
		}

		kt_fisher_exact(
		this->alleles_strand[0][2], // A: Allele on forward strand
		this->alleles_strand[1][2], // B: Allele on reverse strand
		n_cnt_fwd - this->alleles_strand[0][2], // C: Not allele on forward strand
		n_cnt_rev - this->alleles_strand[1][2], // D: Not allele on reverse strand
		&fisher_left_p, &fisher_right_p, &fisher_twosided_p);

		this->d->fs_a[0] = std::abs(-10 * log10(fisher_twosided_p));

		// If n_alleles = 2 then they are identical because of symmetry
		uint8_t pos = 1;
		if(this->n_alleles - 2 > 2){
			for(U32 p = 3; p < this->n_alleles; ++p){
				kt_fisher_exact(
				this->alleles_strand[0][p], // A: Allele on forward strand
				this->alleles_strand[1][p], // B: Allele on reverse strand
				n_cnt_fwd - this->alleles_strand[0][p], // C: Not allele on forward strand
				n_cnt_rev - this->alleles_strand[1][p], // D: Not allele on reverse strand
				&fisher_left_p, &fisher_right_p, &fisher_twosided_p);

				this->d->fs_a[pos++] = std::abs(-10 * log10(fisher_twosided_p));
			}
		}

		// Total number of genotypes is the sum of the root
		// nodes excluding special missing and sentinel node
		// (0 and 1).
		uint64_t n_genotypes = this->gt[2][2].n_cnt + this->gt[2][3].n_cnt + this->gt[3][2].n_cnt + this->gt[3][3].n_cnt;

		// Allele frequency of A
		const double p = ((double)2*this->gt[2][2].n_cnt + this->gt[2][3].n_cnt + this->gt[3][2].n_cnt) / (2*n_genotypes);
		// Genotype frequency of heterozyotes
		const double pg = ((double)this->gt[2][3].n_cnt + this->gt[3][2].n_cnt) / n_genotypes;
		// Expected heterozygosity
		const double exp = 2*p*(1-p);
		// Population inbreeding coefficient: F
		const double f_pic = exp > 0 ? (exp-pg)/exp : 0;
		this->d->f_pic = f_pic;
	}
}

}
