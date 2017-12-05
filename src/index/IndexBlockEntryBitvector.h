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
	IndexBlockEntryBitvector() : n_keys(0), fields_set(0), keys(nullptr), bit_bytes(nullptr){}
	~IndexBlockEntryBitvector(){
		delete [] this->keys;
		delete [] this->bit_bytes;
	}

	inline void update(const BYTE& value, const U32& pos){ this->bit_bytes[pos] = value; }
	inline void allocate(const U32& n_keys, const U32& w){
		delete [] this->keys;
		delete [] this->bit_bytes;
		this->keys = new U32[n_keys];
		this->bit_bytes = new BYTE[w];
		memset(this->bit_bytes, 0, w);
		this->n_keys = n_keys;
	}

	inline void allocate(const U32& w){
		delete [] this->bit_bytes;
		this->bit_bytes = new BYTE[w];
		memset(this->bit_bytes, 0, w);
	}

	template <class T>
	inline const bool operator[](const T& p) const{ return((this->bit_bytes[p / 8] & (1 << (p % 8))) >> (p % 8)); }

	inline const U32* const firstKey(void) const{ return(&this->keys[0]); }
	inline const U32* const lastKey(void) const {
		if(this->n_keys == 0) return(this->firstKey());
		return(&this->keys[this->n_keys - 1]);
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
			stream.write(reinterpret_cast<const char*>(&entry.n_keys), sizeof(U32));
			for(U32 i = 0; i < entry.n_keys; ++i)
				stream.write(reinterpret_cast<const char*>(&entry.keys[i]), sizeof(U32));

			return(stream);
		}

		friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
			stream.read(reinterpret_cast<char*>(&entry.n_keys), sizeof(U32));
			entry.keys = new U32[entry.n_keys];
			for(U32 i = 0; i < entry.n_keys; ++i)
				stream.read(reinterpret_cast<char*>(&entry.keys[i]), sizeof(U32));

			return(stream);
		}

		inline const U32 getBaseSize(void) const{
			return(sizeof(U32) + sizeof(U32)*this->n_keys);
		}

public:
	U32 n_keys;
	U16 fields_set; // when reading only
	U32* keys;
	BYTE* bit_bytes;
};

}
}

#endif /* INDEX_INDEXBLOCKENTRYBITVECTOR_H_ */
