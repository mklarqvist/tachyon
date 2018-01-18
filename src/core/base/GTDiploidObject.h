#ifndef CORE_GTDIPLOIDOBJECT_H_
#define CORE_GTDIPLOIDOBJECT_H_

namespace Tachyon{
namespace Core{

/**<
 * Primary higher-level abstraction object for interpreted
 * diploid GT primitives
 */
struct GTDiploidObject{
private:
	typedef GTDiploidObject self_type;
	typedef Core::MetaEntry meta_entry_type;

public:
	GTDiploidObject(const BYTE n_alleles) :
		phase(0),
		n_objects(0)
	{}
	~GTDiploidObject(){ }

	inline BYTE& operator[](const U32& p){ return(this->alleles[p]); }
	inline BYTE* getAlleles(void){ return(&this->alleles[0]); }
	inline const U64& size(void) const{ return(this->n_objects); }

	void operator()(const U64& gt_primitive, const meta_entry_type& meta_entry){
		const Core::TACHYON_GT_TYPE type = meta_entry.hot.getGenotypeType();
		if(type == Core::YON_GT_RLE_DIPLOID_BIALLELIC){
			const BYTE shift    = meta_entry.hot.controller.gt_anyMissing    ? 2 : 1;
			const BYTE add      = meta_entry.hot.controller.gt_mixed_phasing ? 1 : 0;

			if(add) this->phase = gt_primitive & 1;
			else    this->phase = meta_entry.hot.controller.gt_phase;

			this->alleles[0]    = (gt_primitive & ((1 << shift) - 1) << add) >> add;
			this->alleles[1]    = (gt_primitive & ((1 << shift) - 1) << (add+shift)) >> (add+shift);
			this->n_objects     = gt_primitive >> (2*shift + add);
		} else if(type == Core::YON_GT_RLE_DIPLOID_NALLELIC){
			const BYTE shift    = ceil(log2(meta_entry.cold.n_allele + meta_entry.hot.controller.gt_anyMissing)); // Bits occupied per allele, 1 value for missing
			const BYTE add      = meta_entry.hot.controller.gt_mixed_phasing ? 1 : 0;

			if(add) this->phase = gt_primitive & 1;
			else    this->phase = meta_entry.hot.controller.gt_phase;

			this->alleles[0]    = (gt_primitive & ((1 << shift) - 1) << add) >> add;
			this->alleles[1]    = (gt_primitive & ((1 << shift) - 1) << (add+shift)) >> (add+shift);
			this->n_objects     = gt_primitive >> (2*shift + add);
		} else {
			std::cerr << "not implemented" << std::endl;
			exit(1);
		}
	}

public:
	BYTE phase;
	BYTE alleles[2];
	U64  n_objects;
};

}
}

#endif /* GTDIPLOIDOBJECT_H_ */
