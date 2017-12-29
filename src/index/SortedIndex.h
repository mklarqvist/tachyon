#ifndef INDEX_SORTEDINDEX_H_
#define INDEX_SORTEDINDEX_H_

#include "IndexEntry.h"
#include "IndexIndexEntry.h"

namespace Tachyon{
namespace Index{

class SortedIndex{
private:
	typedef SortedIndex self_type;
	typedef IndexEntry index_entry_type;
	typedef IndexIndexEntry index_index_type;

public:
	SortedIndex() : n_superindex(0), n_index(0){}
	~SortedIndex(){}

	bool buildSuperIndex(void){
		index_index_type indexindex;
		indexindex(this->index[0]);
		for(U32 i = 1; i < this->index.size(); ++i){
			if(indexindex == this->index[i]) indexindex += this->index[i];
			else {
				//std::cerr << indexindex.contigID << ':' << indexindex.minPosition << "-" << indexindex.maxPosition << std::endl;
				this->indexindex.push_back(indexindex);
				indexindex(this->index[i]);
			}
		}

		if(this->indexindex.size() == 0){
			this->indexindex.push_back(indexindex);
		}
		else if(indexindex != this->indexindex.back()){
			this->indexindex.push_back(indexindex);
		}

		return true;
	}

	inline const U32& size(void) const{ return(this->n_index); }
	inline void operator+=(const index_entry_type& entry){ this->index.push_back(entry); }

	const index_entry_type* const at(const U32& p)   const{ return(&this->index[p]); }
	const index_entry_type& operator[](const U32& p) const{ return( this->index[p]); }

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.n_superindex), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.n_index),      sizeof(U32));
		for(U32 i = 0; i < entry.n_superindex; ++i)
			stream << entry.indexindex[i];

		for(U32 i = 0; i < entry.n_index; ++i)
			stream << entry.index[i];

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.n_superindex), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.n_index),      sizeof(U32));

		entry.indexindex.resize(entry.n_superindex);
		entry.index.resize(entry.n_index);

		for(U32 i = 0; i < entry.n_superindex; ++i)
			stream >> entry.indexindex[i];

		for(U32 i = 0; i < entry.n_index; ++i)
			stream >> entry.index[i];

		return(stream);
	}

private:
	U32 n_superindex;
	U32 n_index;

	// Basic
	index_entry_type current_index_entry;
	std::vector<index_entry_type> index;
	std::vector<index_index_type> indexindex;
};

}
}

#endif /* INDEX_SORTEDINDEX_H_ */
