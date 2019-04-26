#include <fstream>

#include "genotypes.h"
#include "support/fisher_math.h"

namespace tachyon{

yon_gt_ppa::yon_gt_ppa(void) : n_s(0), ordering(nullptr) {}
yon_gt_ppa::yon_gt_ppa(const uint32_t n_samples) : n_s(n_samples), ordering(new uint32_t[n_samples]) { this->reset(); }
yon_gt_ppa::~yon_gt_ppa(void) { delete [] this->ordering; }

yon_gt_ppa::yon_gt_ppa(const yon_gt_ppa& other) :
	n_s(other.n_s),
	ordering(new uint32_t[other.n_s])
{
	memcpy(this->ordering, other.ordering, sizeof(uint32_t)*this->n_s);
}

yon_gt_ppa::yon_gt_ppa(yon_gt_ppa&& other) noexcept :
	n_s(other.n_s),
	ordering(other.ordering)
{
	other.ordering = nullptr;
}

yon_gt_ppa& yon_gt_ppa::operator=(const yon_gt_ppa& other) {
	delete [] ordering;
	n_s = other.n_s;
	ordering = new uint32_t[n_s];
	memcpy(this->ordering, other.ordering, sizeof(uint32_t)*this->n_s);
	return(*this);
}

yon_gt_ppa& yon_gt_ppa::operator=(yon_gt_ppa&& other) noexcept{
	if (this == &other) {
		// precautions against self-moves
		return *this;
	}
	delete [] ordering; ordering = nullptr;
	n_s = other.n_s;
	std::swap(ordering, other.ordering);
	return(*this);
}

void yon_gt_ppa::Allocate(const uint32_t n_s) {
	delete [] this->ordering;
	this->n_s = n_s;
	this->ordering = new uint32_t[n_s];
	this->reset();
}

void yon_gt_ppa::reset(void) {
	for (uint32_t i = 0; i < this->n_s; ++i)
		this->ordering[i] = i;
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, yon_gt_ppa& ppa) {
	DeserializePrimitive(ppa.n_s, buffer);
	ppa.ordering = new uint32_t[ppa.n_s];
	for (uint32_t i = 0; i < ppa.n_s; ++i)
		DeserializePrimitive(ppa.ordering[i], buffer);

	return(buffer);
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const yon_gt_ppa& ppa) {
	SerializePrimitive(ppa.n_s, buffer);
	for (uint32_t i = 0; i < ppa.n_s; ++i)
		SerializePrimitive(ppa.ordering[i], buffer);

	return(buffer);
}

yon_radix_gt::yon_radix_gt() :
	n_ploidy(0),
	n_allocated(4),
	id(0),
	alleles(new uint16_t[this->n_allocated])
{
	memset(this->alleles, 0, sizeof(uint16_t)*this->n_allocated);
}

yon_radix_gt::~yon_radix_gt() { delete [] this->alleles; }

yon_radix_gt::yon_radix_gt(const yon_radix_gt& other) :
	n_ploidy(other.n_ploidy),
	n_allocated(other.n_allocated),
	id(other.id),
	alleles(new uint16_t[this->n_allocated])
{
	memcpy(this->alleles, other.alleles, sizeof(uint16_t)*this->n_allocated);
}

yon_radix_gt::yon_radix_gt(yon_radix_gt&& other) :
	n_ploidy(other.n_ploidy),
	n_allocated(other.n_allocated),
	id(other.id),
	alleles(other.alleles)
{
	other.alleles = nullptr;
}

yon_radix_gt& yon_radix_gt::operator=(const yon_radix_gt& other) // copy assignment
{
	this->id = other.id;
	this->n_ploidy = other.n_ploidy;
	this->n_allocated = other.n_allocated;
	delete [] this->alleles;
	this->alleles = new uint16_t[this->n_allocated];
	memcpy(this->alleles, other.alleles, sizeof(uint16_t)*this->n_allocated);
	return *this;
}

yon_radix_gt& yon_radix_gt::operator=(yon_radix_gt&& other) // move assignment
{
	if (this!=&other) // prevent self-move
	{
		this->id = other.id;
		this->n_ploidy = other.n_ploidy;
		this->n_allocated = other.n_allocated;
		delete [] this->alleles;
		this->alleles = other.alleles;
		other.alleles = nullptr;
	}
	return *this;
}

bool yon_radix_gt::operator<(const yon_radix_gt& other) const {
	// Do not compare incremental sample identification
	// numbers as that is not the desired outcome of
	// the sort.
	if (this->n_ploidy < other.n_ploidy) return true;
	if (other.n_ploidy < this->n_ploidy) return false;

	for (uint32_t i = 0; i < this->n_ploidy; ++i) {
		if (this->alleles[i] < other.alleles[i])
			return true;
	}
	return false;
}

bool yon_radix_gt::operator==(const yon_radix_gt& other) const {
	// Do not compare incremental sample identification
	// numbers as that is not the desired outcome of
	// the comparison.
	if (this->n_ploidy != other.n_ploidy)
		return false;

	for (uint32_t i = 0; i < this->n_ploidy; ++i) {
		if (this->alleles[i] != other.alleles[i])
			return false;
	}
	return true;
}

std::ostream& operator<<(std::ostream& stream, const yon_radix_gt& genotype) {
	stream << genotype.id << ":";
	if (genotype.n_ploidy) {
		stream << genotype.alleles[0];
		for (uint32_t i = 1; i < genotype.n_ploidy; ++i) {
			stream << "," << genotype.alleles[i];
		}
	}
	return(stream);
}

uint64_t yon_radix_gt::GetPackedInteger(const uint8_t& shift_size) const {
	uint64_t packed = 0;
	for (uint32_t i = 0; i < this->n_ploidy; ++i) {
		packed <<= shift_size;
		assert(((this->alleles[i] << shift_size) >> shift_size) == this->alleles[i]);
		packed |= (this->alleles[i] & ((1 << shift_size)) - 1);
	}
	return packed;
}

void yon_radix_gt::resize(const uint8_t new_ploidy) {
	uint16_t* temp = new uint16_t[new_ploidy];
	memcpy(temp, this->alleles, this->n_allocated * sizeof(uint16_t));
	delete [] this->alleles;
	this->alleles = temp;
	this->n_allocated = new_ploidy;
}

yon_gt_rcd::yon_gt_rcd() : run_length(0), allele(nullptr) {}
yon_gt_rcd::~yon_gt_rcd() { delete [] this->allele; }
yon_gt_rcd::yon_gt_rcd(yon_gt_rcd&& other) :
		run_length(other.run_length),
		allele(other.allele)
{
	other.allele = nullptr;
}

yon_gt_rcd& yon_gt_rcd::operator=(yon_gt_rcd&& other) {
	if (this == &other) return(*this);
	delete[] allele; allele = nullptr;
	std::swap(allele, other.allele);
	run_length = other.run_length;
	return(*this);
}

yon_buffer_t& yon_gt_rcd::PrintVcf(yon_buffer_t& buffer, const uint8_t& n_ploidy) {
	if (this->allele[0] == 1) {
		buffer += '.';
		return(buffer);
	}
	if (this->allele[0] == 0) buffer += '.';
	else buffer.AddReadble(((this->allele[0] >> 1) - 2));

	for (uint32_t i = 1; i < n_ploidy; ++i) {
		if (this->allele[i] == 1) break;
		buffer += ((this->allele[i] & 1) ? '|' : '/');
		if (this->allele[i] == 0) buffer += '.';
		else buffer.AddReadble(((this->allele[i] >> 1) - 2));
	}
	return(buffer);
}

// Start GT

yon_gt::yon_gt() : eval_cont(0), add(0), global_phase(0), shift(0), p(0), m(0), method(0), n_s(0), n_i(0), n_o(0),
		   n_allele(0), ppa(nullptr), data(nullptr),
		   d_exp(nullptr), rcds(nullptr), n_i_occ(nullptr), d_occ(nullptr),
		   dirty(false)
{}

yon_gt::yon_gt(const yon_gt& other) :
	eval_cont(other.eval_cont), add(other.add), global_phase(other.global_phase), shift(other.shift),
	p(other.p), m(other.m), method(other.method), n_s(other.n_s), n_i(other.n_i), n_o(other.n_o),
	n_allele(other.n_allele), ppa(other.ppa), data(other.data),
	d_exp(nullptr), rcds(nullptr), n_i_occ(nullptr), d_occ(nullptr),
	dirty(other.dirty)
{
	if (other.d_exp != nullptr) {
		d_exp = new yon_gt_rcd[n_s];
		for (int i = 0; i < n_s; ++i) d_exp[i] = std::move(other.d_exp[i].Clone(m));
	}

	if (other.rcds != nullptr) {
		rcds = new yon_gt_rcd[n_i];
		for (int i = 0; i < n_i; ++i) rcds[i] = std::move(other.rcds[i].Clone(m));
	}

	if (n_o) {
		assert(other.n_i_occ != nullptr);
		n_i_occ = new uint32_t[n_o];
		for (int i = 0; i < n_o; ++i) n_i_occ[i] = other.n_i_occ[i];
		assert(other.d_occ != nullptr);
		d_occ = new yon_gt_rcd*[n_o];
		for (int i = 0; i < n_o; ++i) {
			d_occ[i] = new yon_gt_rcd[n_i_occ[i]];
			for (int j = 0; j < n_i_occ[i]; ++j)
				d_occ[i][j] = std::move(other.d_occ[i][j].Clone(m));
		}
	}
}

yon_gt& yon_gt::operator=(const yon_gt& other) {
	// Destroy local data.
	delete [] rcds;
	delete [] d_exp;
	if (d_occ != nullptr) {
		for (uint32_t i = 0; i < this->n_o; ++i)
			delete [] d_occ[i];
	}
	delete [] n_i_occ;
	delete [] d_occ;

	eval_cont = other.eval_cont; add = other.add; global_phase = other.global_phase;
	shift = other.shift; p = other.p; m = other.m; method = other.method;
	n_s = other.n_s; n_i = other.n_i; n_o = other.n_o;
	n_allele = other.n_allele; ppa = other.ppa; data = other.data;
	d_exp = nullptr; rcds = nullptr; n_i_occ = nullptr; d_occ = nullptr;
	dirty = other.dirty;

	if (other.d_exp != nullptr) {
		d_exp = new yon_gt_rcd[n_s];
		for (int i = 0; i < n_s; ++i) d_exp[i] = std::move(other.d_exp[i].Clone(m));
	}

	if (other.rcds != nullptr) {
		rcds = new yon_gt_rcd[n_i];
		for (int i = 0; i < n_i; ++i) rcds[i] = std::move(other.rcds[i].Clone(m));
	}

	if (n_o) {
		assert(other.n_i_occ != nullptr);
		n_i_occ = new uint32_t[n_o];
		for (int i = 0; i < n_o; ++i) n_i_occ[i] = other.n_i_occ[i];
		assert(other.d_occ != nullptr);
		d_occ = new yon_gt_rcd*[n_o];
		for (int i = 0; i < n_o; ++i) {
			d_occ[i] = new yon_gt_rcd[n_i_occ[i]];
			for (int j = 0; j < n_i_occ[i]; ++j)
				d_occ[i][j] = std::move(other.d_occ[i][j].Clone(m));
		}
	}

	return(*this);
}

yon_gt::yon_gt(yon_gt&& other) noexcept :
		eval_cont(other.eval_cont), add(other.add), global_phase(other.global_phase), shift(other.shift),
		p(other.p), m(other.m), method(other.method), n_s(other.n_s), n_i(other.n_i), n_o(other.n_o),
		n_allele(other.n_allele), ppa(other.ppa), data(other.data),
		d_exp(nullptr), rcds(nullptr), n_i_occ(nullptr), d_occ(nullptr),
		dirty(other.dirty)
{
	std::swap(d_exp, other.d_exp);
	std::swap(rcds, other.rcds);
	std::swap(n_i_occ, other.n_i_occ);
	std::swap(d_occ, other.d_occ);
}

yon_gt& yon_gt::operator=(yon_gt&& other) noexcept{
	if (this == &other) {
		// precautions against self-moves
		return *this;
	}

	// Destroy local data.
	delete [] rcds;
	delete [] d_exp;
	if (d_occ != nullptr) {
		for (uint32_t i = 0; i < this->n_o; ++i)
			delete [] d_occ[i];
	}
	delete [] n_i_occ;
	delete [] d_occ;

	eval_cont = other.eval_cont; add = other.add; global_phase = other.global_phase;
	shift = other.shift; p = other.p; m = other.m; method = other.method;
	n_s = other.n_s; n_i = other.n_i; n_o = other.n_o;
	n_allele = other.n_allele; ppa = other.ppa; data = other.data;
	d_exp = nullptr; rcds = nullptr; n_i_occ = nullptr; d_occ = nullptr;
	dirty = other.dirty;

	std::swap(d_exp, other.d_exp);
	std::swap(rcds, other.rcds);
	std::swap(n_i_occ, other.n_i_occ);
	std::swap(d_occ, other.d_occ);

	return(*this);
}


yon_gt::~yon_gt() {
	delete [] rcds;
	delete [] d_exp;
	if (d_occ != nullptr) {
		for (uint32_t i = 0; i < this->n_o; ++i)
			delete [] d_occ[i];
	}
	delete [] n_i_occ;
	delete [] d_occ;
}

bool yon_gt::Evaluate(void) {
	// Prevent double evaluation.
	if (this->eval_cont & YON_GT_UN_RCDS)
		return true;

	if (this->method == 1) return(this->EvaluateRecordsM1());
	else if (this->method == 2) return(this->EvaluateRecordsM2());
	else if (this->method == 4) return(this->EvaluateRecordsM4());
	else {
		std::cerr << "not implemented method " << (int)this->method << std::endl;
	}
	return false;
}

bool yon_gt::EvaluateRecordsM1() {
	// Prevent double evaluation.
	if (this->eval_cont & YON_GT_UN_RCDS)
		return true;

	switch(this->p) {
	case(1): return(this->EvaluateRecordsM1_<uint8_t>());
	case(2): return(this->EvaluateRecordsM1_<uint16_t>());
	case(4): return(this->EvaluateRecordsM1_<uint32_t>());
	case(8): return(this->EvaluateRecordsM1_<uint64_t>());
	default:
		std::cerr << "illegal primitive in EvaluateRecordsM1" << std::endl;
		exit(1);
	}
}

bool yon_gt::EvaluateRecordsM2() {
	// Prevent double evaluation.
	if (this->eval_cont & YON_GT_UN_RCDS)
		return true;

	switch(this->p) {
	case(1): return(this->EvaluateRecordsM2_<uint8_t>());
	case(2): return(this->EvaluateRecordsM2_<uint16_t>());
	case(4): return(this->EvaluateRecordsM2_<uint32_t>());
	case(8): return(this->EvaluateRecordsM2_<uint64_t>());
	default:
		std::cerr << "illegal primitive in EvaluateRecordsM2" << std::endl;
		exit(1);
	}
}

bool yon_gt::EvaluateRecordsM4() {
	// Prevent double evaluation.
	if (this->eval_cont & YON_GT_UN_RCDS)
		return true;

	switch(this->p) {
	case(1): return(this->EvaluateRecordsM4_<uint8_t>());
	case(2): return(this->EvaluateRecordsM4_<uint16_t>());
	case(4): return(this->EvaluateRecordsM4_<uint32_t>());
	case(8): return(this->EvaluateRecordsM4_<uint64_t>());
	default:
		std::cerr << "illegal primitive in EvaluateRecordsM1" << std::endl;
		exit(1);
	}
}

bool yon_gt::ExpandRecordsPpa(void) {
	if (this->eval_cont & YON_GT_UN_EXPAND)
		return true;

	if ((this->eval_cont & YON_GT_UN_RCDS) == false) {
		bool eval = this->Evaluate();
		if (eval == false) return false;
	}

	assert(this->rcds != nullptr);
	assert(this->ppa != nullptr);

	if (this->d_exp != nullptr) delete [] this->d_exp;
	this->d_exp = new yon_gt_rcd[this->n_s];

	uint64_t cum_sample = 0;
	uint64_t cum_offset = 0;
	for (uint32_t i = 0; i < this->n_i; ++i) {
		for (uint32_t j = 0; j < this->rcds[i].run_length; ++j, ++cum_sample) {
			const uint32_t& target_ppa = this->ppa->at(cum_sample);
			this->d_exp[target_ppa] = std::move(this->rcds[i].Clone(this->m));
			this->d_exp[target_ppa].run_length = 1;
			cum_offset += this->p * this->m;
		}
	}
	assert(cum_sample == this->n_s);
	assert(cum_offset == this->n_s * this->m * this->p);
	this->eval_cont |= YON_GT_UN_EXPAND;
	return true;
}

bool yon_gt::ExpandRecordsPpaExternal(yon_gt_rcd* d_expe) {
	if ((this->eval_cont & YON_GT_UN_RCDS) == false) {
		bool eval = this->Evaluate();
		if (eval == false) return false;
	}

	assert(this->rcds != nullptr);
	assert(this->ppa != nullptr);

	uint64_t cum_sample = 0;
	uint64_t cum_offset = 0;
	for (uint32_t i = 0; i < this->n_i; ++i) {
		for (uint32_t j = 0; j < this->rcds[i].run_length; ++j, ++cum_sample) {
			const uint32_t& target_ppa = this->ppa->at(cum_sample);
			d_expe[target_ppa] = std::move(this->rcds[i].Clone(this->m));
			d_expe[target_ppa].run_length = 1;
			cum_offset += this->p * this->m;
		}
	}
	assert(cum_sample == this->n_s);
	assert(cum_offset == this->n_s * this->m * this->p);
	return true;
}

bool yon_gt::Expand(void) {
	if (this->ppa != nullptr)
		return(this->ExpandRecordsPpa());
	else return(this->ExpandRecords());
}

bool yon_gt::ExpandExternal(yon_gt_rcd* d_expe) {
	if (this->ppa != nullptr)
		return(this->ExpandRecordsPpaExternal(d_expe));
	else return(this->ExpandRecordsExternal(d_expe));
}

bool yon_gt::ExpandRecords(void) {
	if ((this->eval_cont & YON_GT_UN_RCDS) == false) {
		bool eval = this->Evaluate();
		if (eval == false) return false;
	}

	if (this->eval_cont & YON_GT_UN_EXPAND)
		return true;

	assert(this->rcds != nullptr);

	if (this->d_exp != nullptr) delete [] this->d_exp;
	this->d_exp = new yon_gt_rcd[this->n_s];

	uint64_t cum_sample = 0;
	uint64_t cum_offset = 0;
	for (uint32_t i = 0; i < this->n_i; ++i) {
		for (uint32_t j = 0; j < this->rcds[i].run_length; ++j, ++cum_sample) {
			this->d_exp[cum_sample] = std::move(this->rcds[i].Clone(this->m));
			this->d_exp[cum_sample].run_length = 1;
			cum_offset += this->p * this->m;
		}
	}
	assert(cum_sample == this->n_s);
	assert(cum_offset == this->n_s * this->m * this->p);
	this->eval_cont |= YON_GT_UN_EXPAND;
	return true;
}

bool yon_gt::ExpandRecordsExternal(yon_gt_rcd* d_expe) {
	if ((this->eval_cont & YON_GT_UN_RCDS) == false) {
		bool eval = this->Evaluate();
		if (eval == false) return false;
	}

	assert(this->rcds != nullptr);

	uint64_t cum_sample = 0;
	uint64_t cum_offset = 0;
	for (uint32_t i = 0; i < this->n_i; ++i) {
		for (uint32_t j = 0; j < this->rcds[i].run_length; ++j, ++cum_sample) {
			d_expe[cum_sample] = std::move(this->rcds[i].Clone(this->m));
			d_expe[cum_sample].run_length = 1;
			cum_offset += this->p * this->m;
		}
	}
	assert(cum_sample == this->n_s);
	assert(cum_offset == this->n_s * this->m * this->p);
	return true;
}

bcf1_t* yon_gt::UpdateHtslibGenotypes(bcf1_t* rec, bcf_hdr_t* hdr) {
   // Prevent double evaluation.
   if ((this->eval_cont & YON_GT_UN_EXPAND)) {
	   bool check = this->Expand();
	   if (check == false) {
		   std::cerr << utility::timestamp("ERROR","GT") << "Failed to lazy-expand genotype records..." << std::endl;
		   return rec;
	   }
   }

   assert(this->d_exp != nullptr);

   int32_t* tmpi = new int32_t[this->n_s*this->m];
   uint32_t gt_offset = 0;
   for (uint32_t i = 0; i < this->n_s; ++i) {
	   for (uint32_t j = 0; j < this->m; ++j, ++gt_offset) {
		   if (this->d_exp[i].allele[j] == 0)      tmpi[gt_offset] = 0;
		   else if (this->d_exp[i].allele[j] == 1) tmpi[gt_offset] = 1;
		   else tmpi[gt_offset] = YON_GT_BCF1(this->d_exp[i].allele[j]);
	   }
   }
   assert(gt_offset == this->n_s*this->m);

   bcf_update_genotypes(hdr, rec, tmpi, this->n_s*this->m);

   delete [] tmpi;
   return(rec);
}

// summary record
yon_gt_summary_rcd::yon_gt_summary_rcd() :
	n_ploidy(0), n_ac_af(0), n_fs(0), ac(nullptr), af(nullptr),
	nm(0), npm(0), an(0), fs_a(nullptr), hwe_p(0),
	f_pic(0), heterozygosity(0)
{}

yon_gt_summary_rcd::yon_gt_summary_rcd(const yon_gt_summary_rcd& other) :
	n_ploidy(other.n_ploidy), n_ac_af(other.n_ac_af), n_fs(other.n_fs), ac(nullptr), af(nullptr),
	nm(other.nm), npm(other.npm), an(other.an), fs_a(nullptr), hwe_p(other.hwe_p),
	f_pic(other.f_pic), heterozygosity(other.heterozygosity)
{
	if (other.ac != nullptr) { ac = new uint32_t[n_ac_af]; memcpy(ac, other.ac, n_ac_af*sizeof(uint32_t)); }
	if (other.af != nullptr) { af = new float[n_ac_af]; memcpy(af, other.af, n_ac_af*sizeof(float)); }
	if (other.fs_a != nullptr) { fs_a = new float[n_fs]; memcpy(fs_a, other.fs_a, n_fs*sizeof(float)); }
}

yon_gt_summary_rcd& yon_gt_summary_rcd::operator=(const yon_gt_summary_rcd& other) {
	delete [] ac; ac = nullptr;
	delete [] af; af = nullptr;
	delete [] fs_a; fs_a = nullptr;
	n_ploidy = other.n_ploidy; n_ac_af = other.n_ac_af; n_fs = other.n_fs;
	nm = other.nm; npm = other.npm; an = other.an; hwe_p = other.hwe_p;
	f_pic = other.f_pic; heterozygosity = other.heterozygosity;
	if (other.ac != nullptr) { ac = new uint32_t[n_ac_af]; memcpy(ac, other.ac, n_ac_af*sizeof(uint32_t)); }
	if (other.af != nullptr) { af = new float[n_ac_af]; memcpy(af, other.af, n_ac_af*sizeof(float)); }
	if (other.fs_a != nullptr) { fs_a = new float[n_fs]; memcpy(fs_a, other.fs_a, n_fs*sizeof(float)); }

	return(*this);
}

yon_gt_summary_rcd::yon_gt_summary_rcd(yon_gt_summary_rcd&& other) noexcept :
	n_ploidy(other.n_ploidy), n_ac_af(other.n_ac_af), n_fs(other.n_fs), ac(nullptr), af(nullptr),
	nm(other.nm), npm(other.npm), an(other.an), fs_a(nullptr), hwe_p(other.hwe_p),
	f_pic(other.f_pic), heterozygosity(other.heterozygosity)
{
	std::swap(ac, other.ac);
	std::swap(af, other.af);
	std::swap(fs_a, other.fs_a);
}

yon_gt_summary_rcd& yon_gt_summary_rcd::operator=(yon_gt_summary_rcd&& other) noexcept {
	if (this == &other) {
		// precautions against self-moves
		return *this;
	}

	n_ploidy = other.n_ploidy; n_ac_af = other.n_ac_af; n_fs = other.n_fs;
	nm = other.nm; npm = other.npm; an = other.an; hwe_p = other.hwe_p;
	f_pic = other.f_pic; heterozygosity = other.heterozygosity;
	std::swap(ac, other.ac);
	std::swap(af, other.af);
	std::swap(fs_a, other.fs_a);

	return(*this);
}

yon_gt_summary_rcd::~yon_gt_summary_rcd() {
	delete [] ac; delete [] af;
	delete [] fs_a;
	// Do not delete ac_p as it is borrowed
}

// summary object
yon_gt_summary_obj::yon_gt_summary_obj() : n_cnt(0), children(nullptr) {}
yon_gt_summary_obj::yon_gt_summary_obj(yon_gt_summary_obj&& other) noexcept : n_cnt(other.n_cnt), children(nullptr) { std::swap(children, other.children); }
yon_gt_summary_obj& yon_gt_summary_obj::operator=(yon_gt_summary_obj&& other) noexcept{
	if (this == &other) return *this;
	delete [] children; children = nullptr;
	n_cnt = other.n_cnt;
	std::swap(children, other.children);
	return(*this);
}
yon_gt_summary_obj::~yon_gt_summary_obj() { delete [] this->children; }

// Clone helper ctor.
yon_gt_summary_obj::yon_gt_summary_obj(const uint8_t n_alleles, const uint64_t cnt, yon_gt_summary_obj* c) :
	n_cnt(cnt), children(nullptr)
{
	if (c != nullptr) {
		children = new yon_gt_summary_obj[n_alleles];
		for (int i = 0; i < n_alleles; ++i) {
			children[i] = std::move(c[i].Clone(n_alleles));
		}
	}
}

// genotype summary
yon_gt_summary::yon_gt_summary(void) :
	n_ploidy(0),
	n_alleles(0),
	alleles(nullptr),
	alleles_strand(nullptr),
	gt(nullptr),
	d(nullptr)
{

}

yon_gt_summary::yon_gt_summary(const uint8_t base_ploidy, const uint8_t n_alleles) :
	n_ploidy(base_ploidy),
	n_alleles(n_alleles + 2),
	alleles(new uint32_t[this->n_alleles]),
	alleles_strand(new uint32_t*[this->n_ploidy]),
	gt(new yon_gt_summary_obj[this->n_alleles]),
	d(nullptr)
{
	memset(this->alleles, 0, sizeof(uint32_t)*this->n_alleles);
	for (uint32_t i = 0; i < this->n_ploidy; ++i) {
		this->alleles_strand[i] = new uint32_t[this->n_alleles];
		memset(this->alleles_strand[i], 0, sizeof(uint32_t)*this->n_alleles);
	}

	// Add layers to the root node.
	for (uint32_t i = 0; i < this->n_alleles; ++i)
		this->AddGenotypeLayer(&this->gt[i], 1);
}

yon_gt_summary::yon_gt_summary(yon_gt_summary&& other) noexcept :
	n_ploidy(other.n_ploidy), n_alleles(other.n_alleles),
	alleles(nullptr), alleles_strand(nullptr), gt(nullptr), d(nullptr)
{
	std::swap(alleles, other.alleles);
	std::swap(alleles_strand, other.alleles_strand);
	std::swap(gt, other.gt);
	std::swap(d, other.d);
}

yon_gt_summary& yon_gt_summary::operator=(yon_gt_summary&& other) noexcept {
	if (this == &other) {
		// precautions against self-moves
		return *this;
	}

	// Clear local data without cosideration.
	delete [] this->alleles;
	delete [] this->gt;
	if (this->alleles_strand != nullptr) {
		for (uint32_t i = 0; i < this->n_ploidy; ++i)
			delete this->alleles_strand[i];
		delete [] this->alleles_strand;
	}
	delete this->d;

	// Move.
	n_ploidy = other.n_ploidy; n_alleles = other.n_alleles;
	alleles = nullptr; alleles_strand = nullptr; gt = nullptr; d = nullptr;
	std::swap(alleles, other.alleles);
	std::swap(alleles_strand, other.alleles_strand);
	std::swap(gt, other.gt);
	std::swap(d, other.d);

	return(*this);
}

yon_gt_summary::yon_gt_summary(const yon_gt_summary& other) :
	n_ploidy(other.n_ploidy), n_alleles(other.n_alleles),
	alleles(new uint32_t[n_alleles]),
	alleles_strand(new uint32_t*[n_ploidy]),
	gt(new yon_gt_summary_obj[n_alleles]),
	d(nullptr)
{
	for (int i = 0; i < n_alleles; ++i) alleles[i] = other.alleles[i];
	for (uint32_t i = 0; i < n_ploidy; ++i) {
		alleles_strand[i] = new uint32_t[n_alleles];
		for (int j = 0; j < n_alleles; ++j) alleles_strand[i][j] = other.alleles_strand[i][j];
	}
	if (other.d != nullptr) this->d = new yon_gt_summary_rcd(*other.d);
	for (int i = 0; i < n_alleles; ++i) gt[i] = std::move(other.gt[i].Clone(n_alleles));
}

yon_gt_summary& yon_gt_summary::operator=(const yon_gt_summary& other) {
	// Delete previous data without consideration.
	delete [] this->alleles;
	delete [] this->gt;
	if (this->alleles_strand != nullptr) {
		for (uint32_t i = 0; i < this->n_ploidy; ++i)
			delete this->alleles_strand[i];
		delete [] this->alleles_strand;
	}
	delete this->d; this->d = nullptr;

	n_ploidy = other.n_ploidy; n_alleles = other.n_alleles;
	alleles = new uint32_t[n_alleles];
	alleles_strand = new uint32_t*[n_ploidy];
	gt = new yon_gt_summary_obj[n_alleles];

	for (int i = 0; i < n_alleles; ++i) alleles[i] = other.alleles[i];
	for (uint32_t i = 0; i < this->n_ploidy; ++i) {
		alleles_strand[i] = new uint32_t[n_alleles];
		for (int j = 0; j < n_alleles; ++j) alleles_strand[i][j] = other.alleles_strand[i][j];
	}
	if (other.d != nullptr) this->d = new yon_gt_summary_rcd(*other.d);
	for (int i = 0; i < n_alleles; ++i) gt[i] = std::move(other.gt[i].Clone(n_alleles));

	return(*this);
}

void yon_gt_summary::Setup(const uint8_t base_ploidy, const uint8_t n_als) {
	n_ploidy  = base_ploidy;
	n_alleles = n_als + 2;
	// Destroy previous data without consideration.
	delete [] this->alleles;
	delete [] this->gt;
	if (this->alleles_strand != nullptr) {
		for (uint32_t i = 0; i < this->n_ploidy; ++i)
			delete this->alleles_strand[i];
		delete [] this->alleles_strand;
	}
	delete this->d;

	alleles = new uint32_t[n_alleles];
	alleles_strand = new uint32_t*[this->n_ploidy];
	gt = new yon_gt_summary_obj[this->n_alleles];
	d = nullptr;

	memset(this->alleles, 0, sizeof(uint32_t)*this->n_alleles);
	for (uint32_t i = 0; i < this->n_ploidy; ++i) {
		this->alleles_strand[i] = new uint32_t[this->n_alleles];
		memset(this->alleles_strand[i], 0, sizeof(uint32_t)*this->n_alleles);
	}

	// Add layers to the root node.
	for (uint32_t i = 0; i < this->n_alleles; ++i)
		this->AddGenotypeLayer(&this->gt[i], 1);
}

yon_gt_summary::~yon_gt_summary() {
	delete [] this->alleles;
	delete [] this->gt;
	if (this->alleles_strand != nullptr) {
		for (uint32_t i = 0; i < this->n_ploidy; ++i)
			delete this->alleles_strand[i];
		delete [] this->alleles_strand;
	}
	delete this->d;
}

bool yon_gt_summary::AddGenotypeLayer(yon_gt_summary_obj* target, uint8_t depth) {
	if (depth == this->n_ploidy) return false;
	if (target->children == nullptr) {
		target->children = new yon_gt_summary_obj[this->n_alleles];
	}

	for (uint32_t i = 0; i < this->n_alleles; ++i)
		this->AddGenotypeLayer(&target->children[i], depth+1);

	return false;
}

yon_gt_summary& yon_gt_summary::operator+=(const yon_gt& gt) {
	assert(gt.rcds != nullptr);

	// Iterate over available genotype records.
	for (uint32_t i = 0; i < gt.n_i; ++i) {
		// Target root node.
		yon_gt_summary_obj* target = &this->gt[gt.rcds[i].allele[0] >> 1];
		target->n_cnt += gt.rcds[i].run_length;

		// Iterate over alleles given the base ploidy.
		for (uint32_t j = 0; j < gt.m; ++j) {
			assert((gt.rcds[i].allele[j] >> 1) < this->n_alleles);
			// Add allelic counts.
			this->alleles[gt.rcds[i].allele[j] >> 1] += gt.rcds[i].run_length;
			// Add strand-specific (ploidy-aware) alleleic counts.
			this->alleles_strand[j][gt.rcds[i].allele[j] >> 1] += gt.rcds[i].run_length;
		}

		// Update remainder.
		for (uint32_t j = 1; j < gt.m; ++j) {
			target = &target->children[gt.rcds[i].allele[j] >> 1];
			target->n_cnt += gt.rcds[i].run_length;
		}
	}
	return(*this);
}

yon_gt_summary& yon_gt_summary::Add(const yon_gt& gt, const uint32_t n_items, const yon_gt_rcd* rcds) {
	assert(rcds != nullptr);

	// Iterate over available genotype records.
	uint32_t n_total_rle = 0;
	for (uint32_t i = 0; i < n_items; ++i) {
		// Target root node.
		yon_gt_summary_obj* target = &this->gt[rcds[i].allele[0] >> 1];
		target->n_cnt += rcds[i].run_length;

		// Iterate over alleles given the base ploidy.
		for (uint32_t j = 0; j < gt.m; ++j) {
			assert((rcds[i].allele[j] >> 1) < this->n_alleles);
			// Add allelic counts.
			this->alleles[rcds[i].allele[j] >> 1] += rcds[i].run_length;
			// Add strand-specific (ploidy-aware) alleleic counts.
			this->alleles_strand[j][rcds[i].allele[j] >> 1] += rcds[i].run_length;
			n_total_rle += rcds[i].run_length;
		}

		// Update remainder.
		for (uint32_t j = 1; j < gt.m; ++j) {
			target = &target->children[rcds[i].allele[j] >> 1];
			target->n_cnt += rcds[i].run_length;
		}
	}
	return(*this);
}

std::vector<uint64_t> yon_gt_summary::GetAlleleCounts(void) const {
	std::vector<uint64_t> c_allele(this->n_alleles, 0);
	for (uint32_t i = 0; i < this->n_alleles; ++i)
		c_allele[i] = this->alleles[i];

	return(c_allele);
}

std::vector< std::pair<uint64_t,double> > yon_gt_summary::GetAlleleCountFrequency(void) const {
	std::vector< std::pair<uint64_t,double> > c_allele(this->n_alleles);
	uint64_t n_total = 0;
	for (uint32_t i = 0; i < 2; ++i) {
		c_allele[i].first = this->alleles[i];
	}

	for (uint32_t i = 2; i < this->n_alleles; ++i) {
		c_allele[i].first = this->alleles[i];
		n_total += this->alleles[i];
	}

	if (n_total != 0) {
		for (uint32_t i = 0; i < this->n_alleles; ++i) {
			c_allele[i].second = (double)c_allele[i].first / n_total;
		}
	}

	return(c_allele);
}

std::vector< std::vector<uint64_t> > yon_gt_summary::GetAlleleStrandCounts(void) const {
	std::vector< std::vector<uint64_t> > c_allele(this->n_ploidy, std::vector<uint64_t>(this->n_alleles));
	for (uint32_t i = 0; i < this->n_ploidy; ++i) {
		for (uint32_t j = 0; j < this->n_alleles; ++j) {
			c_allele[i][j] = this->alleles_strand[i][j];
		}
	}

	return(c_allele);
}

bool yon_gt_summary::GetGenotype(std::vector<uint64_t>& data,
                                 yon_gt_summary_obj* target,
                                 uint8_t depth) const
{
	if (depth + 1 == this->n_ploidy) {
		assert(target->children != nullptr);
		for (uint32_t i = 0; i < this->n_alleles; ++i) {
			if (target->children[i].n_cnt != 0) {
				data.push_back(target->children[i].n_cnt);
			}
		}
		return false;
	}

	for (uint32_t i = 0; i < this->n_alleles; ++i) {
		this->GetGenotype(data, &target->children[i], depth + 1);
	}

	return(true);
}

// Todo: unfinished.
std::vector<yon_gt_rcd> yon_gt_summary::GetGenotypeCounts(bool drop_empty) const {
	std::vector<yon_gt_rcd> genotypes;

	// Traverse the trie and store records when hitting the leafs.
	std::vector<uint64_t> d;

	// If the target ploidy is greater than one (haploid) we collect
	// the genotypes by traversing the trie. Otherwise the genotypes
	// are the allele counts at the roots.
	if (this->n_ploidy > 1) {
		for (uint32_t i = 0; i < this->n_alleles; ++i) {
			//std::cerr << "outer " << i << "/" << (int)this->n_alleles << std::endl;
			this->GetGenotype(d, &this->gt[i], 1);
		}
	} else {
		for (uint32_t i = 0; i < this->n_alleles; ++i) {
			if (this->gt[i].n_cnt != 0) {
				d.push_back(this->gt[i].n_cnt);
			}
		}
	}
	uint64_t n_total = 0;
	for (uint32_t i = 0; i < d.size(); ++i)
		n_total += d[i];
	std::cerr << "collected: " << d.size() << " total = " << n_total << std::endl;
	return(genotypes);
}

std::vector<double> yon_gt_summary::GetStrandBiasAlleles(const bool phred_scale) const {
	if (this->n_ploidy != 2 || this->n_alleles + 2 < 2)
		return std::vector<double>();

	std::vector<double> strand_bias_p_values;
	double fisher_left_p, fisher_right_p, fisher_twosided_p;

	uint64_t n_cnt_fwd = 0;
	uint64_t n_cnt_rev = 0;
	for (uint32_t i = 2; i < this->n_alleles; ++i) {
		n_cnt_fwd += this->alleles_strand[0][i];
		n_cnt_rev += this->alleles_strand[1][i];
	}

	kt_fisher_exact(
	this->alleles_strand[0][2], // A: Allele on forward strand
	this->alleles_strand[1][2], // B: Allele on reverse strand
	n_cnt_fwd - this->alleles_strand[0][2], // C: Not allele on forward strand
	n_cnt_rev - this->alleles_strand[1][2], // D: Not allele on reverse strand
	&fisher_left_p, &fisher_right_p, &fisher_twosided_p);

	if (phred_scale) strand_bias_p_values.push_back(std::abs(-10 * log10(fisher_twosided_p)));
	else strand_bias_p_values.push_back(fisher_twosided_p);

	// If n_alleles = 2 then they are identical because of symmetry
	if (this->n_alleles - 2 > 2) {
		for (uint32_t p = 3; p < this->n_alleles; ++p) {
			kt_fisher_exact(
			this->alleles_strand[0][p], // A: Allele on forward strand
			this->alleles_strand[1][p], // B: Allele on reverse strand
			n_cnt_fwd - this->alleles_strand[0][p], // C: Not allele on forward strand
			n_cnt_rev - this->alleles_strand[1][p], // D: Not allele on reverse strand
			&fisher_left_p, &fisher_right_p, &fisher_twosided_p);

			if (phred_scale) strand_bias_p_values.push_back(std::abs(-10 * log10(fisher_twosided_p)));
			else strand_bias_p_values.push_back(fisher_twosided_p);
		}
	}
	return(strand_bias_p_values);
}

double yon_gt_summary::CalculateHardyWeinberg(void) const {
	if (this->n_ploidy != 2 || this->n_alleles - 2 != 2) return -1;

	uint64_t obs_hets = this->gt[2][3].n_cnt + this->gt[3][2].n_cnt; // alts
	uint64_t obs_hom1 = this->gt[2][2].n_cnt; // hom ref
	uint64_t obs_hom2 = this->gt[3][3].n_cnt; // hom alt

	uint64_t obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
	uint64_t obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

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
	for (curr_hets = mid; curr_hets > 1; curr_hets -= 2) {
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
	for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2) {
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
	for (i = 0; i <= rare_copies; i++) {
		if (het_probs[i] > het_probs[obs_hets])
			continue;

		p_hwe += het_probs[i];
	}

	p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

	delete [] het_probs;

	return(p_hwe);
}

bool yon_gt_summary::LazyEvaluate(void) {
	// Delete previous data without consideration.
	delete this->d;
	this->d = new yon_gt_summary_rcd;

	// Allele count and frequency
	this->d->n_ploidy = this->n_ploidy;
	this->d->n_ac_af  = this->n_alleles;
	this->d->ac       = new uint32_t[this->n_alleles];
	this->d->af       = new float[this->n_alleles];

	uint64_t n_total = 0;
	for (uint32_t i = 0; i < 2; ++i) this->d->ac[i] = this->alleles[i];
	for (uint32_t i = 2; i < this->n_alleles; ++i) {
		this->d->ac[i] = this->alleles[i];
		n_total += this->alleles[i];
	}

	if (n_total != 0) {
		for (uint32_t i = 0; i < this->n_alleles; ++i)
			this->d->af[i] = (double)this->d->ac[i] / n_total;

	} else {
		for (uint32_t i = 0; i < this->n_alleles; ++i)
			this->d->af[i] = 0;
	}

	this->d->npm = this->d->ac[0];
	this->d->nm  = n_total;
	this->d->an  = n_total + this->d->ac[0];

	if (n_total)
		this->d->hwe_p = this->CalculateHardyWeinberg();
	else
		this->d->hwe_p = 1;

	// Strand-specific bias and inbreeding coefficient (F-statistic)
	if (this->n_ploidy == 2) {
		// Total number of genotypes is the sum of the root
		// nodes excluding special missing and sentinel node
		// (0 and 1).
		const uint32_t n_total_gt = (this->gt[2][2].n_cnt + this->gt[2][3].n_cnt + this->gt[3][2].n_cnt + this->gt[3][3].n_cnt);

		if (n_total_gt)
			this->d->heterozygosity = ((double)this->gt[2][3].n_cnt + this->gt[3][2].n_cnt) / n_total_gt;
		else
			this->d->heterozygosity = 0;

		uint8_t n_fs_used = (this->n_alleles - 2 == 2 ? 1 : this->n_alleles - 2);
		this->d->n_fs = n_fs_used;
		this->d->fs_a = new float[n_fs_used];
		double fisher_left_p, fisher_right_p, fisher_twosided_p;
		uint64_t n_cnt_fwd = 0;
		uint64_t n_cnt_rev = 0;
		for (uint32_t i = 2; i < this->n_alleles; ++i) {
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
		if (this->n_alleles - 2 > 2) {
			for (uint32_t p = 3; p < this->n_alleles; ++p) {
				kt_fisher_exact(
				this->alleles_strand[0][p], // A: Allele on forward strand
				this->alleles_strand[1][p], // B: Allele on reverse strand
				n_cnt_fwd - this->alleles_strand[0][p], // C: Not allele on forward strand
				n_cnt_rev - this->alleles_strand[1][p], // D: Not allele on reverse strand
				&fisher_left_p, &fisher_right_p, &fisher_twosided_p);

				this->d->fs_a[pos++] = std::abs(-10 * log10(fisher_twosided_p));
			}
		}

		if (n_total_gt) {
			// Allele frequency of A
			const double p = ((double)2*this->gt[2][2].n_cnt + this->gt[2][3].n_cnt + this->gt[3][2].n_cnt) / (2*n_total_gt);
			// Genotype frequency of heterozyotes
			const double pg = ((double)this->gt[2][3].n_cnt + this->gt[3][2].n_cnt) / n_total_gt;
			// Expected heterozygosity
			const double exp = 2*p*(1-p);
			// Population inbreeding coefficient: F
			this->d->f_pic = (exp > 0 ? (exp-pg)/exp : 0);
		}
		else this->d->f_pic = 0;
	}

	return true;
}

bool yon_occ::ReadTable(const std::string file_name,
                        const yon_vnt_hdr_t& header,
                        const char delimiter)
{
	std::ifstream f;
	f.open(file_name);
	if (!f.good()) {
		std::cerr << utility::timestamp("ERROR") << "Stream is bad! Cannot open " << file_name << "..." << std::endl;
		return false;
	}

	uint32_t n_line = 0;
	std::string line;
	// Iterate over available lines in the input file.
	while (getline(f, line)) {
		//std::cerr << line << std::endl;

		// Tokenize string with delimiter.
		std::vector<std::string> params = utility::split(line, delimiter, false);

		// Assert that the first column is an existing sample name.
		const int32_t sample_id = header.GetSampleId(params[0]);
		if (sample_id < 0) {
			std::cerr << utility::timestamp("WARNING") << "Cannot find sample \"" << params[0] << "\" in groupings file..." << std::endl;
			continue;
		}

		// Iterate over tokens.
		for (uint32_t i = 1; i < params.size(); ++i) {
			map_type::const_iterator it = this->map.find(params[i]);
			if (it == this->map.end()) {
				// Not already set
				this->table.push_back(std::vector<uint32_t>(header.GetNumberSamples() + 1, 0));
				this->table.back()[sample_id + 1] = true;
				this->map[params[i]] = this->row_names.size();
				this->row_names.push_back(params[i]);
				//std::cerr << "Adding group: " << params[i] << " for " << this->table.back().size() << " samples" << std::endl;
			} else {
				// Already set
				this->table[it->second][sample_id + 1] = true;
			}
		}
	}

	return true;
}

bool yon_occ::BuildTable(void) {
	this->occ.clear();
	this->vocc.clear();
	if (this->table.size() == 0) {
		return false;
	}

	this->occ = std::vector< std::vector<uint32_t> >(this->table.size(), std::vector<uint32_t>( this->table[0].size(), 0));
	this->cum_sums = std::vector< uint32_t >( this->occ.size() );

	for (uint32_t i = 0; i < this->table.size(); ++i) {
		assert(this->table[i][0] == 0);
		for (uint32_t j = 1; j < this->occ[i].size(); ++j)
			this->occ[i][j] += this->occ[i][j-1] + this->table[i][j];

		this->cum_sums[i] = this->occ[i].back();
	}

	// Matrix transpose for faster random access lookups.
	this->vocc = std::vector< std::vector<uint32_t> >(this->table[0].size() , std::vector<uint32_t>(this->occ.size(), 0));
	for (int i = 0; i < this->table[0].size(); ++i) {
		for (int j = 0; j < this->occ.size(); ++j) {
			vocc[i][j] = occ[j][i];
		}
	}

	return true;
}

bool yon_occ::BuildTable(const yon_gt_ppa* ppa_p) {
	if (ppa_p == nullptr)
		return(this->BuildTable());

	this->occ.clear();
	this->vocc.clear();
	if (this->table.size() == 0) {
		return false;
	}

	// Convert pointer to reference.
	const yon_gt_ppa& ppa = *ppa_p;

	assert(ppa.n_s + 1 == this->table[0].size());

	this->occ = std::vector< std::vector<uint32_t> >(this->table.size(), std::vector<uint32_t>( this->table[0].size(), 0));
	this->cum_sums = std::vector< uint32_t >( this->occ.size() );

	for (uint32_t i = 0; i < this->table.size(); ++i) {
		assert(this->table[i][0] == 0);
		for (uint32_t j = 1; j < this->occ[i].size(); ++j)
			this->occ[i][j] += this->occ[i][j - 1] + this->table[i][ppa[j - 1] + 1];

		assert(this->occ[i][0] == 0);
		this->cum_sums[i] = this->occ[i].back();
	}

	// Matrix transpose for faster random access lookups.
	this->vocc = std::vector< std::vector<uint32_t> >(this->table[0].size() , std::vector<uint32_t>(this->occ.size(), 0));
	for (int i = 0; i < this->table[0].size(); ++i) {
		for (int j = 0; j < this->occ.size(); ++j) {
			vocc[i][j] = occ[j][i];
		}
	}

	return true;
}

// support
yon_vnt_cnt::yon_vnt_cnt(void) :
	gt_available(0),
	gt_has_missing(0),
	gt_phase_uniform(0),
	gt_has_mixed_phasing(0),
	gt_compression_type(0),
	gt_primtive_type(0),
	gt_mixed_ploidy(0),
	biallelic(0),
	simple_snv(0),
	diploid(0),
	alleles_packed(0),
	all_snv(0)
{}

void yon_vnt_cnt::operator=(const uint16_t& value) {
	const yon_vnt_cnt* const other = reinterpret_cast<const yon_vnt_cnt* const>(&value);
	this->gt_available      = other->gt_available;
	this->gt_has_missing    = other->gt_has_missing;
	this->gt_phase_uniform  = other->gt_phase_uniform;
	this->gt_has_mixed_phasing = other->gt_has_mixed_phasing;
	this->gt_compression_type  = other->gt_compression_type;
	this->gt_primtive_type = other->gt_primtive_type;
	this->gt_mixed_ploidy  = other->gt_mixed_ploidy;
	this->biallelic        = other->biallelic;
	this->simple_snv       = other->simple_snv;
	this->diploid          = other->diploid;
	this->alleles_packed   = other->alleles_packed;
	this->all_snv          = other->all_snv;
}

}
