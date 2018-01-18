#ifndef CORE_ITERATOR_METAHOTITERATOR_H_
#define CORE_ITERATOR_METAHOTITERATOR_H_

namespace Tachyon{
namespace Iterator{

class MetaHotIterator{
	typedef MetaHotIterator self_type;
	typedef Core::StreamContainer container_type;
	typedef Core::MetaHot entry_type;

public:
	MetaHotIterator() :
		__flag_data_literal(false),
		n_entries(0),
		current_position(0),
		entries(nullptr)
	{}

	MetaHotIterator(const container_type& container) :
		__flag_data_literal(false),
		n_entries(0),
		current_position(0),
		entries(nullptr)
	{
		this->setup(container);
	}

	MetaHotIterator(const container_type& container, const Core::TACHYON_GT_TYPE gt_retain_type) :
		__flag_data_literal(false),
		n_entries(0),
		current_position(0),
		entries(nullptr)
	{

	}

	~MetaHotIterator(){
		if(this->__flag_data_literal)
			delete [] this->entries;
	}

	bool setupFilterGenotypePrimitive(const container_type& container, const Core::TACHYON_GT_TYPE gt_retain_type){
		if(container.buffer_data_uncompressed.pointer == 0)
			return false;

		this->clear();
		this->__flag_data_literal = true;

		const entry_type* const entries_all = reinterpret_cast<const entry_type* const>(container.buffer_data_uncompressed.data);
		assert((container.buffer_data_uncompressed.pointer % sizeof(entry_type)) == 0);
		const U64 n_entries_all = container.buffer_data_uncompressed.pointer / sizeof(entry_type);

		for(U32 i = 0; i < n_entries_all; ++i){
			if(entries_all[i].getGenotypeType() == gt_retain_type)
				++this->n_entries;
		}
		this->entries = new entry_type[this->n_entries];

		if(gt_retain_type == Core::YON_GT_RLE_DIPLOID_BIALLELIC){
			for(U32 i = 0; i < n_entries_all; ++i){
				if(entries_all[i].getGenotypeType() == Core::YON_GT_RLE_DIPLOID_BIALLELIC)
					this->entries[i] = entries_all[i];
			}
		} else if(gt_retain_type == Core::YON_GT_RLE_DIPLOID_NALLELIC){
			for(U32 i = 0; i < n_entries_all; ++i){
				if(entries_all[i].getGenotypeType() == Core::YON_GT_RLE_DIPLOID_NALLELIC)
					this->entries[i] = entries_all[i];
			}
		} else {
			std::cerr << "not implemented" << std::endl;
			exit(1);
		}
		std::cerr << "Total entries: " << n_entries_all << '\t' << this->n_entries << std::endl;

		return(true);
	}

	bool setup(const container_type& container){
		if(container.buffer_data_uncompressed.pointer == 0)
			return false;

		this->clear();
		this->entries = reinterpret_cast<entry_type* const>(container.buffer_data_uncompressed.data);
		assert((container.buffer_data_uncompressed.pointer % sizeof(entry_type)) == 0);
		this->n_entries = container.buffer_data_uncompressed.pointer / sizeof(entry_type);
		return true;
	}

	inline bool operator()(const container_type& container){
		return(this->setup(container));
	}

	void clear(void){
		this->n_entries = 0;
		this->current_position = 0;
		if(this->__flag_data_literal)
			delete [] this->entries;
		this->__flag_data_literal = false;
		entries = nullptr;
	}

	bool nextEntry(const entry_type*& entry){
		if(this->current_position == this->n_entries)
			return false;

		entry = &this->entries[this->current_position];
		++this->current_position;
		return true;
	}

	bool previousEntry(const entry_type*& entry){
		if(this->current_position == 0)
			return false;

		entry = &this->entries[this->current_position];
		--this->current_position;
		return true;
	}

	inline void operator++(void){ if(this->current_position != this->n_entries) ++this->current_position; }
	inline void operator--(void){ if(this->current_position != 0) --this->current_position; }
	inline void operator+=(const S32& p){ if(this->current_position + p <= this->n_entries) this->current_position += p; }
	inline void operator-=(const S32& p){ if(this->current_position - p >= 0) this->current_position -= p; }

	inline const entry_type& operator[](const U32& p) const{ return(this->entries[p]); }
	inline const size_t size(void) const{ return(this->n_entries); }
	inline const entry_type& first(void) const{ return(this->entries[0]); }
	inline const entry_type& last(void) const{
		if(this->n_entries == 0) return(this->entries[0]);
		return(this->entries[this->n_entries - 1]);
	}

private:
	bool        __flag_data_literal;
	S32         n_entries;
	S32         current_position;
	entry_type* entries;
};

}
}

#endif /* CORE_ITERATOR_METAHOTITERATOR_H_ */
