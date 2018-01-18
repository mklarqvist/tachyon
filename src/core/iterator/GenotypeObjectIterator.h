#ifndef ITERATOR_GENOTYPEOBJECTITERATOR_H_
#define ITERATOR_GENOTYPEOBJECTITERATOR_H_

namespace Tachyon{
namespace Iterator{

/**<
 * Primary iterator class for genotype objects. At the
 * moment the primary usecase involves passing U64
 * primitives to the iterator
 */
class GenotypeObjectIteratorInterface{
private:
	typedef GenotypeObjectIteratorInterface self_type;
	typedef Core::StreamContainer           container_type;
	typedef Core::BlockEntry                block_type;
	typedef MetaIterator                    meta_iterator_type;
	typedef ContainerIterator               container_iterator_type;
	typedef Core::MetaEntry                 meta_entry_type;
	typedef IntegerIterator                 integer_iterator_type;

public:
	GenotypeObjectIteratorInterface();
	virtual ~GenotypeObjectIteratorInterface();

	virtual BYTE getAllele(const U32& p) const =0;
	virtual void countAlleles() =0;
	virtual void countAllelesGroup() =0;
	virtual void getSummary() =0;
	virtual void getSummaryGroup() =0;

	//virtual void operator[](const U32& p) const =0;
	//virtual void at(const U32& p) const =0;
	//virtual std::vector<void> getGenotypes(void) =0;

	inline const U32 size(void) const{ return(this->n_width); }
	inline const U32 capacity(void) const{ return(this->n_capacity); }
	inline const U64 getNumberSamples(void) const{ return(this->n_samples); }

	// Todo: want to have a MetaIterator here with only the meta objects for the current GT type!!!!

private:
	template <class T>
	bool __updatePrimitive(const T& gt_primitive);

private:
	U64 n_samples;
	U32 current_position;
	U32 n_capacity;
	U32 n_width;
	// groupings
	U32 n_groups;

	// if RLE we have to iterate over each object
	U32 current_position_object;
	U32 n_current_object_length;
};

class GenotypeDiploidIterator : public GenotypeObjectIteratorInterface{
private:
	typedef GenotypeDiploidIterator self_type;
	typedef Core::GTDiploidObject   gt_type;

public:
	GenotypeDiploidIterator();
	~GenotypeDiploidIterator();

	void countAlleles();
	void countAllelesGroup();
	void getSummary();

private:
	gt_type* __gt_entries;
};

}
}

#endif /* GENOTYPEOBJECTITERATOR_H_ */
