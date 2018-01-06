#ifndef CORE_ITERATOR_GENOTYPEITERATOR_H_
#define CORE_ITERATOR_GENOTYPEITERATOR_H_

namespace Tachyon{
namespace Core{
namespace Iterator{

/**
 * We have several different GT representations
 *
 * Case diploid, bi-allelic, no NA, no missing, and run-length encoded
 * Case diploid, bi-allelic, no NA, has missing, and run-length encoded
 * Case diploid, n-allelic or has NA
 * Case diploid, n-allelic or has NA run-length encoded
 * Case m-ploid, n-allelic
 */
class GenotypeIteratorIterface{
private:
	typedef GenotypeIteratorIterface self_type;

public:
	GenotypeIteratorIterface();
	virtual ~GenotypeIteratorIterface();

	// std::vector<some_abstract_genotype_object> toVector(void) const;

private:
	virtual void toVCFString(std::ostream& stream) const =0;

};

class GenotypeIteratorDiploidBiallelicRLE : public GenotypeIteratorIterface{}
class GenotypeIteratorDiploidRLE : public GenotypeIteratorIterface{}
class GenotypeIteratorDiploid : public GenotypeIteratorIterface{}

}
}
}


#endif /* CORE_ITERATOR_GENOTYPEITERATOR_H_ */
