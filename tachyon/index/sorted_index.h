#ifndef INDEX_SORTED_INDEX_H_
#define INDEX_SORTED_INDEX_H_
#include <bitset>

#include "index_entry.h"
#include "index_index_entry.h"

namespace tachyon{
namespace index{

class SortedIndex{
private:
	typedef SortedIndex       self_type;
    typedef std::size_t       size_type;
    typedef IndexEntry        value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
	typedef IndexIndexEntry   meta_type;
	typedef meta_type*        pointer_meta;

public:
	SortedIndex() :
		n_superindex(0),
		n_index(0),
		n_capacity(1000),
		index(new value_type[1000]),
		indexindex(new meta_type[1000])
	{}
	~SortedIndex(){}

	/**<
	 * Builds the meta-index of index entries when
	 * the input data is sorted
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool buildSuperIndex(void){
		if(this->size() == 0)
			return false;

		meta_type indexindex;
		indexindex(this->index[0]);
		for(U32 i = 1; i < this->size(); ++i){
			if(indexindex == this->index[i])
				indexindex += this->index[i];
			else {
				this->indexindex[this->n_superindex++] = indexindex;
				indexindex(this->index[i]);
			}
		}

		if(this->size_meta() == 0){
			this->indexindex[this->n_superindex - 1] = indexindex;
		}
		else if(indexindex != this->indexindex[this->n_superindex]){
			this->indexindex[this->n_superindex - 1] = indexindex;
		}

		return true;
	}

	inline self_type& operator+=(const const_reference index_entry){
		if(this->n_index + 1 == this->n_capacity)
			this->resize();


		this->index[this->n_index++] = index_entry;
		return(*this);
	}

	void resize(void){
		pointer temp = this->index;
		this->index = new value_type[this->n_capacity*2];
		for(U32 i = 0; i < this->n_index; ++i)
			this->index[i] = temp[i];

		delete [] temp;
		this->n_capacity *= 2;
	}

	class iterator{
	private:
		typedef iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		reference operator*() const{ return *ptr_; }
		pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	class const_iterator{
	private:
		typedef const_iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		const_iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		const_reference operator*() const{ return *ptr_; }
		const_pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	// Element access
	inline reference at(const size_type& position){ return(this->index[position]); }
	inline const_reference at(const size_type& position) const{ return(this->index[position]); }
	inline reference operator[](const size_type& position){ return(this->index[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->index[position]); }
	inline pointer data(void){ return(this->index); }
	inline const_pointer data(void) const{ return(this->index); }
	inline reference front(void){ return(this->index[0]); }
	inline const_reference front(void) const{ return(this->index[0]); }
	inline reference back(void){ return(this->index[this->n_index - 1]); }
	inline const_reference back(void) const{ return(this->index[this->n_index - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_index == 0); }
	inline const size_type& size(void) const{ return(this->n_index); }
	inline const size_type& size_meta(void) const{ return(this->n_superindex); }

	// Iterator
	inline iterator begin(){ return iterator(&this->index[0]); }
	inline iterator end(){ return iterator(&this->index[this->n_index - 1]); }
	inline const_iterator begin() const{ return const_iterator(&this->index[0]); }
	inline const_iterator end() const{ return const_iterator(&this->index[this->n_index - 1]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->index[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->index[this->n_index - 1]); }

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		//const U32 n_superindex = entry.indexindex.size();
		//const U32 n_index      = entry.index.size();

		stream.write(reinterpret_cast<const char*>(&entry.n_superindex), sizeof(size_type));
		stream.write(reinterpret_cast<const char*>(&entry.n_index),      sizeof(size_type));

		for(U32 i = 0; i < entry.n_superindex; ++i)
			stream << entry.indexindex[i];

		for(U32 i = 0; i < entry.n_index; ++i)
			stream << entry.index[i];

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.n_superindex), sizeof(size_type));
		stream.read(reinterpret_cast<char*>(&entry.n_index),      sizeof(size_type));

		entry.index = new value_type[entry.n_index];
		entry.indexindex = new meta_type[entry.n_superindex];
		//entry.indexindex.resize(entry.n_superindex);
		//entry.index.resize(entry.n_index);

		for(U32 i = 0; i < entry.n_superindex; ++i)
			stream >> entry.indexindex[i];

		for(U32 i = 0; i < entry.n_index; ++i)
			stream >> entry.index[i];

		return(stream);
	}

private:
	size_type n_superindex;
	size_type n_index;
	size_type n_capacity;

	// Basic
	pointer index;
	pointer_meta indexindex;
};

}
}

#endif /* INDEX_SORTED_INDEX_H_ */
