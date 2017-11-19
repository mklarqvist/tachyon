#ifndef CORE_HASHCONTAINER_H_
#define CORE_HASHCONTAINER_H_

#include "../algorithm/OpenHashTable.h"

namespace Tachyon{
namespace Core{
namespace Support{

class HashContainer{
	typedef HashContainer self_type;
	typedef Hash::HashTable<U32, U32> hash_table;
	typedef std::vector<U32> id_vector;

public:
	HashContainer() : htable(16384){}
	~HashContainer(){}

	inline const bool get(const U32& value, U32& ret){
		U32* ret2 = nullptr;
		if(this->htable.GetItem(&value, ret2, sizeof(U32))){
			ret = this->data[*ret2];
			return true;
		}
		return false;
	}

	inline bool getRaw(const U32& value, U32& ret){
		U32* ret2 = nullptr;
		if(this->htable.GetItem(&value, ret2, sizeof(U32))){
			ret = *ret2;
			return true;
		}
		return false;
	}

	inline bool set(const U32& value){
		U32 tot = this->data.size();
		this->htable.SetItem(&value, tot, sizeof(U32));
		this->data.push_back(value);
		return true;
	}

	inline U32 setGet(const U32& value){
		U32* ret = nullptr;
		if(this->htable.GetItem(&value, ret, sizeof(U32))){
			return(*ret);
		} else {
			U32 tot = this->data.size();
			this->htable.SetItem(&value, tot, sizeof(U32));
			this->data.push_back(value);
			return(tot);
		}
	}

	inline const size_t size(void) const{ return(this->data.size()); }
	inline const U32& operator[](const U32& p) const{ return(this->data[p]); }
	inline void clear(void){
		this->htable.clear();
		this->data.clear();
	}

private:
	id_vector data;
	hash_table htable;
};

class HashVectorContainer{
	typedef HashVectorContainer self_type;
	typedef Hash::HashTable<U64, U32> hash_table;
	typedef std::vector< std::vector<U32> > id_vector;

public:
	HashVectorContainer() : htable(16384){}
	~HashVectorContainer(){}

	inline bool get(const U64& value, std::vector<U32>& ret){
		U32* ret2 = nullptr;
		if(this->htable.GetItem(&value, ret2, sizeof(U64))){
			ret = this->data[*ret2];
			return true;
		}
		return false;
	}

	inline bool getRaw(const U64& value, U32& ret){
		U32* ret2 = nullptr;
		if(this->htable.GetItem(&value, ret2, sizeof(U64))){
			ret = *ret2;
			return true;
		}
		return false;
	}

	inline bool set(const std::vector<U32>& value, const U64& hashValue){
		U32 tot = this->data.size();
		this->htable.SetItem(&hashValue, tot, sizeof(U64));
		this->data.push_back(value);
		return true;
	}

	inline U32 setGet(const std::vector<U32>& value, const U64& hashValue){
		U32* ret = nullptr;
		if(this->htable.GetItem(&hashValue, ret, sizeof(U64))){
			return(*ret);
		} else {
			U32 tot = this->data.size();
			this->htable.SetItem(&hashValue, tot, sizeof(U64));
			this->data.push_back(value);
			return(tot);
		}
	}

	inline const size_t size(void) const{ return(this->data.size()); }
	inline const std::vector<U32>& operator[](const U32& p) const{ return(this->data[p]); }
	inline void clear(void){
		this->htable.clear();
		this->data.clear();
	}

private:
	id_vector data;
	hash_table htable;
};

}
}
}



#endif /* CORE_HASHCONTAINER_H_ */
