#ifndef INDEX_INDEXBLOCKENTRYBITVECTOR_H_
#define INDEX_INDEXBLOCKENTRYBITVECTOR_H_

namespace Tachyon{
namespace Index{

// Size of entries in these records are
// inferred from the number of INFO/FORMAT/FILTER
// entries in all the records in a block
struct IndexBlockEntryBitvector{
	typedef IndexBlockEntryBitvector self_type;

public:
	IndexBlockEntryBitvector() : fields_set(0), bit_bytes(nullptr){}
	~IndexBlockEntryBitvector(){ delete [] this->bit_bytes; }

	inline void update(const BYTE& value, const U32& pos){ this->bit_bytes[pos] = value; }
	inline void allocate(const U32& w){ this->bit_bytes = new BYTE[w]; memset(this->bit_bytes, 0, w); }

	template <class T>
	inline const bool operator[](const T& p) const{ return((this->bit_bytes[p / 8] & (1 << (p % 8))) >> (p % 8)); }

public:
	U16 fields_set; // when reading only
	BYTE* bit_bytes;
};

}
}

#endif /* INDEX_INDEXBLOCKENTRYBITVECTOR_H_ */
