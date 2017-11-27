#ifndef CORE_TACHYONREADER_H_
#define CORE_TACHYONREADER_H_

#include "BlockEntry.h"

namespace Tachyon{
namespace Core{

class TachyonReader{
public:

	TachyonReader() : filesize(0){}
	TachyonReader(const std::string& filename) : input_file(filename), filesize(0){}
	~TachyonReader(){}

	bool open(void);
	bool open(const std::string& filename){
		if(filename.size() == 0){
			std::cerr << "no filename" << std::endl;
			return false;
		}
		this->stream.open(filename, std::ios::binary | std::ios::in | std::ios::ate);
		this->filesize = (U64)this->stream.tellg();
		if(!this->stream.good()){
			std::cerr << "failed to read file" << std::endl;
			return false;
		}
		this->stream.seekg(0);
		if(!this->stream.good()){
			std::cerr << "failed to rewrind" << std::endl;
			return false;
		}
		return true;
	}

	bool nextBlock(){
		if(this->stream.tellg() == this->filesize){
			std::cerr << "eof all done" << std::endl;
			return false;
		}

		this->stream >> this->block;
		return true;
	}
	bool seekBlock(const U32& b);

public:
	std::string input_file;
	std::ifstream stream;
	U64 filesize;
	Core::BlockEntry block;
};

}
}

/*
#include "base/EntryHotMeta.h"
#include "base/EntryColdMeta.h"
#include "base/GTRecords.h"

namespace Tachyon{

template <class T>
class TachyonBlockIterator{
	typedef TachyonBlockIterator self_type;
	typedef Support::EntryHotMeta<T> meta_type;
	typedef Support::EntryColdMeta meta_complex_type;
	typedef Totempole::IndexEntry totempole_entry_type;

public:
	TachyonBlockIterator(char* data, const U64& size, const totempole_entry_type& totempole) :
		position(0),
		p_rle(0),
		p_simple(0),
		pointer(0),
		upper_limit(0),
		width(size),
		totempole(totempole),
		data(data),
		meta(reinterpret_cast<const meta_type*>(data)),
		encoding_RLE(&data[0]),
		encoding_simple(&data[0]),
		meta_complex(&data[0])
	{
		this->upper_limit = this->meta[0].n_runs;
	}

	~TachyonBlockIterator(){}

	bool operator++(void){
		if(this->position + 1 == this->totempole.n_variants) return false;
		++this->position;
		this->pointer = 0;

		this->upper_limit = this->getMeta().n_runs;
		if(this->getMeta().controller.biallelic){
			this->encoding_RLE = &this->data[0];
			++this->p_rle;
		}
		else {
			this->encoding_simple = &this->data[0];
			++this->p_simple;
		}

		return true;
	}

	inline const bool isRLE(void) const{ return(this->meta[this->position].isRLE()); }
	inline const U32& size(void) const{ return(this->upper_limit); }

	template <class S>
	const bool nextRun(const S*& run){
		if(this->pointer == this->upper_limit)
			return false;

		//run = this->encoding_RLE;
		run = reinterpret_cast<const S*>(this->encoding_RLE);
		++this->pointer;
		this->encoding_RLE += sizeof(S);
		return true;
	}

	template <class S, BYTE missing = 1>
	const bool nextRun(const Support::TachyonRun<S, missing>*& run){
		if(this->pointer == this->upper_limit)
			return false;

		//run = this->encoding_RLE;
		run = reinterpret_cast<const Support::TachyonRun<S, missing>*>(this->encoding_RLE);
		++this->pointer;
		this->encoding_RLE += sizeof(S);
		return true;
	}

	template <class S, BYTE missing = 1>
	const bool nextRun(const Support::TachyonRunNoPhase<S, missing>*& run){
		if(this->pointer == this->upper_limit)
			return false;

		//run = this->encoding_RLE;
		run = reinterpret_cast<const Support::TachyonRunNoPhase<S, missing>*>(this->encoding_RLE);
		++this->pointer;
		this->encoding_RLE += sizeof(S);
		return true;
	}

	template <class S>
	const bool nextRunSimple(const Support::TachyonRunSimple<S>*& field){
		if(this->pointer == this->upper_limit)
			return false;

		field = reinterpret_cast<const Support::TachyonRunSimple<S>*>(this->encoding_simple);
		++this->pointer;
		this->encoding_simple += sizeof(S);
		return true;
	}

	template <class S>
	const bool nextRunSimple(const S*& field){
		if(this->pointer == this->upper_limit)
			return false;

		field = reinterpret_cast<const S*>(this->encoding_simple);
		++this->pointer;
		this->encoding_simple += sizeof(S);
		return true;
	}

	inline const meta_type& getMeta(void) const{ return(this->meta[this->position]); }
	inline meta_complex_type& getMetaComplex(void){
		return(*reinterpret_cast<meta_complex_type*>(&this->meta_complex[this->meta[this->position].virtual_offset_cold_meta]));
	}

	bool countGenotypes(void);
	bool countGenotypesGroup(void);

private:
	bool __countGenotypesRLE(void);
	bool __countGenotypesRLEGroup(void);
	bool __countGenotypesSimple(void);
	bool __countGenotypesSimpleGroup(void);

private:
	U32 position;    // current meta position
	U32 p_rle;       // position RLE
	U32 p_simple;    // position simple
	U32 pointer;     // iterator pointer
	U32 upper_limit; // relative upper bounds in iterator
	const U64 width;
	const totempole_entry_type& totempole;
	const char* const data;
	const meta_type* meta;
	const char* encoding_RLE;
	const char* encoding_simple;
	char* meta_complex;
};

}
*/

#endif /* CORE_TACHYONREADER_H_ */
