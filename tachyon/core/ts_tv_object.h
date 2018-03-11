#ifndef CORE_TS_TV_OBJECT_H_
#define CORE_TS_TV_OBJECT_H_

#include "../support/type_definitions.h"
#include "meta_entry.h"

namespace tachyon{
namespace core{

struct TsTvObject{
private:
	typedef TsTvObject self_type;
	typedef MetaEntry  meta_type;

public:
	TsTvObject() : n_transitions(0), n_transversions(0)
	{
		for(U32 i = 0; i < 5; ++i){
			this->base_conversions[i] = new U64[5];
			memset(&this->base_conversions[i][0], 0, sizeof(U64)*5);
		}
	}

	~TsTvObject(){
		for(U32 i = 0; i < 5; ++i)
			delete [] this->base_conversions[i];
	}

	self_type& operator+=(const self_type& other){
		this->n_transitions += other.n_transitions;
		this->n_transversions += other.n_transversions;
		for(U32 i = 0; i < 5; ++i){
			for(U32 j = 0; j < 5; ++j){
				this->base_conversions[i][j] += other.base_conversions[i][j];
			}
		}
		return(*this);
	}

	inline const double getTiTVRatio(void) const{
		if(this->n_transversions + this->n_transitions == 0)
			return 0;

		return((double)this->n_transitions / this->n_transversions);
	}

private:
	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		U64 n_total_variants = 0;
		out << entry.n_transversions << '\t' << entry.n_transitions << '\t' << entry.getTiTVRatio() << '\t';
		for(U32 i = 0; i < 4; ++i){
			for(U32 j = 0; j < 4; ++j){
				out << entry.base_conversions[i][j] << '\t';
				if(i != j) n_total_variants += entry.base_conversions[i][j];
			}
		}
		out << n_total_variants;
		return(out);
	}

public:
	// Base -> Base
	U64  n_transitions;
	U64  n_transversions;
	U64* base_conversions[5]; // {A,T,G,C,N} -> {A,T,G,C,N}. A 5x5 matrix
};

}
}



#endif /* CORE_TS_TV_OBJECT_H_ */
