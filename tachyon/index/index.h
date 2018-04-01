#ifndef INDEX_INDEX_H_
#define INDEX_INDEX_H_
#include <bitset>

#include "index_entry_container.h"
#include "index_meta_container.h"

namespace tachyon{
namespace index{

class Index{
private:
	typedef Index                 self_type;
    typedef std::size_t           size_type;
	typedef IndexEntryContainer   container_type;
	typedef IndexMetaContainer    container_meta_type;
	typedef IndexEntry            entry_type;
	typedef IndexIndexEntry       entry_meta_type;

public:
	Index(){}
	~Index(){}

	/**<
	 * Builds the meta-index of index entries when
	 * the input data is sorted
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool buildSuperIndex(void){
		if(this->index_.size() == 0)
			return false;

		entry_meta_type indexindex;
		indexindex(this->index_[0]);
		for(U32 i = 1; i < this->index_.size(); ++i){
			if(indexindex == this->index_[i])
				indexindex += this->index_[i];
			else {
				this->index_meta_ += indexindex;
				indexindex(this->index_[i]);
			}
		}

		if(this->index_meta_.size() == 0){
			this->index_meta_ += indexindex;
		}
		else if(indexindex != this->index_meta_.back()){
			this->index_meta_.back() = indexindex;
		}

		return true;
	}

	// Capacity
	inline const bool empty(void) const{ return(this->index_.empty()); }
	const size_t size(void) const{ return(this->index_.size()); }
	const size_t sizeMeta(void) const{ return(this->index_meta_.size()); }

	inline void operator+=(const entry_type& entry){ this->index_ += entry; }
	inline void operator+=(const entry_meta_type& entry){ this->index_meta_ += entry; }

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.index_;
		stream << entry.index_meta_;
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.index_;
		stream >> entry.index_meta_;
		return(stream);
	}

private:
	container_type      index_;
	container_meta_type index_meta_;
};

}
}

#endif /* INDEX_INDEX_H_ */
