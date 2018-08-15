#include "vcf_utils.h"

namespace tachyon {
namespace io {

const bcf_hrec_t* GetPopulatedHrec(const bcf_idpair_t& idPair) {
	for (int i = 0; i < 3; i++) {
		const bcf_hrec_t* hrec = idPair.val->hrec[i];
		if (hrec != nullptr) {
			return hrec;
		}
	}
	std::cerr << "No populated hrec in idPair. Error in htslib." << std::endl;
	return nullptr;
}

}
}
