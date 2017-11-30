#ifndef CORE_DECORATOR_METAHOTDECORATOR_H_
#define CORE_DECORATOR_METAHOTDECORATOR_H_

namespace Tachyon{
namespace Core{
namespace Decorator{

class MetaHotDecorator{
	typedef MetaHotDecorator self_type;
	typedef Core::StreamContainer container_type;
	typedef Core::MetaHot entry_type;

public:
	MetaHotDecorator() : n_entries(0), pos(0), curpos(0), position_offset(0), entries(nullptr){}
	MetaHotDecorator(container_type& container, const U64& offset) : n_entries(0), pos(0), curpos(0), position_offset(offset), entries(nullptr){
		this->set(container, offset);
	}
	~MetaHotDecorator(){}
	bool set(container_type& container, const U64& offest){
		if(container.buffer_data_uncompressed.pointer == 0)
			return false;

		this->clear();
		std::cerr << container.buffer_data_uncompressed.pointer << std::endl;
		this->entries = reinterpret_cast<const entry_type* const>(container.buffer_data_uncompressed.data);
		std::cout << this->entries[0] << std::endl;
		std::cout << &this->entries[1] - this->entries << '\t' << sizeof(entry_type) << std::endl;
		std::cout << this->entries[1] << std::endl;
		std::cout << (container.buffer_data_uncompressed.pointer % sizeof(entry_type)) << std::endl;
		this->position_offset = offest;
		assert((container.buffer_data_uncompressed.pointer % sizeof(entry_type)) == 0);
		this->n_entries = container.buffer_data_uncompressed.pointer / sizeof(entry_type);
		return true;
	}
	inline bool operator()(container_type& container, const U64& offset){
		return(this->set(container,offset));
	}
	void clear(void){
		this->n_entries = 0;
		this->pos = 0;
		this->curpos = 0;
		this->position_offset = 0;
		entries = nullptr;
	}

	bool nextEntry(const entry_type*& entry){
		if(this->pos == this->n_entries)
			return false;

		entry = &this->entries[this->pos];
		++this->pos;
		return true;
	}
	bool previousEntry(const entry_type*& entry){
		if(this->pos == 0)
			return false;

		entry = &this->entries[this->pos];
		--this->pos;
		return true;
	}

	inline void operator++(void){ if(this->pos != this->n_entries) ++this->pos; }
	inline void operator--(void){ if(this->pos != 0) --this->pos; }
	inline void operator+=(const S32& p){ if(this->pos + p <= this->n_entries) this->pos += p; }
	inline void operator-=(const S32& p){ if(this->pos - p >= 0) this->pos -= p; }

	inline const entry_type& operator[](const U32& p) const{ return(this->entries[p]); }
	inline const size_t size(void) const{ return(this->n_entries); }
	inline const entry_type& first(void) const{ return(this->entries[0]); }
	inline const entry_type& last(void) const{
		if(this->n_entries == 0) return(this->entries[0]);
		return(this->entries[this->n_entries - 1]);
	}

private:
	size_t n_entries;
	size_t pos;
	size_t curpos;
	U64 position_offset;
	const entry_type* entries;
};

}
}
}

#endif /* CORE_DECORATOR_METAHOTDECORATOR_H_ */
