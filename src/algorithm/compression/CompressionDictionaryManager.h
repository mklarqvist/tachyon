#ifndef ALGORITHM_COMPRESSION_COMPRESSIONDICTIONARYMANAGER_H_
#define ALGORITHM_COMPRESSION_COMPRESSIONDICTIONARYMANAGER_H_

#include "zstd.h"
#include "dictBuilder/zdict.h"

namespace Tachyon{
namespace Algorithm{
namespace Compression{

struct CompressionDictionaryManagerPair{
private:
	typedef CompressionDictionaryManagerPair self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Core::StreamContainer stream_container;

public:
	CompressionDictionaryManagerPair() : available(false){}
	CompressionDictionaryManagerPair(const U32 size) : available(false), buffer_data(size){}
	~CompressionDictionaryManagerPair(){
		this->buffer_data.deleteAll();
		this->buffer_dictionary.deleteAll();
	}

	void resize(const U32 size){ this->buffer_data.resize(size); }

	inline const bool isAvailable(void) const{ return(this->available); }

	void operator+=(const stream_container& container){
		if(container.header.controller.uniform == false){
			if(this->buffer_data.pointer < 100e6){
				this->buffer_data += container.buffer_data;
				this->lengths.push_back(container.buffer_data.pointer);
			}
		}

	}

	bool buildDictionary(void){
		if(this->buffer_data.pointer < 100e3){
			std::cerr << "no data" << std::endl;
			return false;
		}
		this->buffer_dictionary.resize(this->buffer_data.pointer * 1.2);


		const size_t ret = ZDICT_trainFromBuffer(this->buffer_dictionary.data, this->buffer_dictionary.capacity(), this->buffer_data.data, &this->lengths[0], this->lengths.size());
		if(ZSTD_isError(ret)){
			std::cerr << "error zstd: " << ZSTD_getErrorString(ZSTD_getErrorCode(ret)) << std::endl;
			return false;
		}
		this->buffer_dictionary.pointer = ret;
		if((float)this->buffer_dictionary.pointer/this->buffer_data.pointer < 0.1)
			this->available = true;

		std::cerr << "Dictionary: " << this->buffer_data.pointer << '/' << this->buffer_dictionary.pointer << '\t' << (float)this->buffer_dictionary.pointer/this->buffer_data.pointer << "\t" << (this->available ? "PASS" : "SKIP") << std::endl;

		return true;
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		const U32 size = entry.field_name.size();
		stream.write(reinterpret_cast<const char* const>(&size), sizeof(U32));
		stream.write(&entry.field_name[0], entry.field_name.size());
		stream.write(reinterpret_cast<const char* const>(&entry.buffer_dictionary.pointer), sizeof(U64));
		stream.write(entry.buffer_dictionary.data, entry.buffer_dictionary.pointer);
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		entry.available = true;
		U32 size;
		stream.read(reinterpret_cast<char*>(&size), sizeof(U32));
		entry.field_name.resize(size);
		stream.read(&entry.field_name[0], size);
		stream.read(reinterpret_cast<char*>(&entry.buffer_dictionary.pointer), sizeof(U64));
		stream.read(entry.buffer_dictionary.data, entry.buffer_dictionary.pointer);
		return(stream);
	}

public:
	bool available;
	std::string field_name;
	std::vector<size_t> lengths;
	buffer_type buffer_data;
	buffer_type buffer_dictionary;
};

class CompressionDictionaryManager{
private:
	typedef CompressionDictionaryManager self_type;
	typedef IO::BasicBuffer buffer_type;
	typedef CompressionDictionaryManagerPair pair_type;

public:
	CompressionDictionaryManager() : load_byte_limit(100e6), min_byte_limit(100e3), n_map_fields(0), map_fields(nullptr){}
	~CompressionDictionaryManager(){
		delete [] this->map_fields;
	}

	void resizeAll(const U32 size){
		this->meta_hot.resize(size);
		this->meta_cold.resize(size);
		this->gt_rle.resize(size);
		this->gt_simple.resize(size);
		if(this->map_fields != nullptr){
			for(U32 i = 0; i < this->n_map_fields; ++i)
				this->map_fields[i].resize(size);
		}
	}

	void allocateFields(const U32 n_fields, const U32 size){
		delete [] this->map_fields;
		this->n_map_fields = n_fields;
		this->map_fields = new pair_type[n_fields];
		for(U32 i = 0; i < n_fields; ++i)
			this->map_fields[i].resize(size);
	}

public:
	U32 load_byte_limit;
	U32 min_byte_limit;
	U32 n_map_fields;
	pair_type meta_hot;
	pair_type meta_cold;
	pair_type gt_rle;
	pair_type gt_simple;
	pair_type* map_fields;
};

}
}
}

#endif /* ALGORITHM_COMPRESSION_COMPRESSIONDICTIONARYMANAGER_H_ */
