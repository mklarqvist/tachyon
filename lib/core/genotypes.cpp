#include "genotypes.h"

namespace tachyon{

yon_gt::~yon_gt(){
	delete [] d_bcf;
	delete [] d_bcf_ppa,
		delete [] rcds;
	delete [] d_exp;
	delete itree;
}

}
