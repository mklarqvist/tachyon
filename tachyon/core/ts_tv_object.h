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
	TsTvObject() :
		n_insertions(0),
		n_deletions(0),
		n_singletons(0),
		n_transitions(0),
		n_transversions(0)
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
		this->n_insertions    += other.n_insertions;
		this->n_deletions     += other.n_deletions;
		this->n_singletons    += other.n_singletons;
		this->n_transitions   += other.n_transitions;
		this->n_transversions += other.n_transversions;
		for(U32 i = 0; i < 5; ++i){
			for(U32 j = 0; j < 5; ++j){
				this->base_conversions[i][j] += other.base_conversions[i][j];
			}
		}
		return(*this);
	}

	inline const double getTiTVRatio(void) const{
		// Prevent division by 0
		if(this->n_transversions == 0) return 0;
		return((double)this->n_transitions / this->n_transversions);
	}

private:
	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		U64 n_total_variants = 0;
		out << entry.n_transversions << '\t' << entry.n_transitions << '\t' << entry.getTiTVRatio() << '\t';
		// Skip 5: encoding for N
		for(U32 i = 0; i < 4; ++i){
			for(U32 j = 0; j < 4; ++j){
				out << entry.base_conversions[i][j] << '\t';
				if(i != j) n_total_variants += entry.base_conversions[i][j];
			}
		}
		out << n_total_variants << '\t' << entry.n_insertions;
		return(out);
	}

public:
	// Base -> Base
	U32  n_insertions;
	U32  n_deletions;
	U32  n_singletons;
	U32  n_transitions;
	U32  n_transversions;
	U64* base_conversions[5]; // {A,T,G,C,N} -> {A,T,G,C,N}. A 5x5 matrix
};

}
}



#endif /* CORE_TS_TV_OBJECT_H_ */
