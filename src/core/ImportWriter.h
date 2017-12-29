#ifndef CORE_TOMAHAWKIMPORTWRITER_H_
#define CORE_TOMAHAWKIMPORTWRITER_H_

#include <cassert>
#include <fstream>

#include "../support/TypeDefinitions.h"
#include "../index/IndexEntry.h"

namespace Tachyon {

inline bool bytePreprocessor(const U32* const data, const size_t& size, char* destination){
	if(size == 0) return false;

	BYTE* targetData = reinterpret_cast<BYTE*>(destination);
	BYTE* s1 = &targetData[size*3];
	BYTE* s2 = &targetData[size*2];
	BYTE* s3 = &targetData[size*1];
	BYTE* s4 = &targetData[size*0];

	for(U32 i = 0; i < size; ++i){
		const BYTE* const p = reinterpret_cast<const BYTE* const>(&data[i]);
		s1[i] = (p[0] & 255);
		s2[i] = (p[1] & 255);
		s3[i] = (p[2] & 255);
		s4[i] = (p[3] & 255);
		const U32 x = s4[i] << 24 | s3[i] << 16 | s2[i] << 8 | s1[i];
		assert(x == data[i]);
	}

	return true;
}

inline bool bytePreprocessorRevert(const char* const data, const size_t& size, char* destination){
	if(size == 0) return false;

	const BYTE* d = reinterpret_cast<const BYTE*>(data);
	U32* dest = reinterpret_cast<U32*>(destination);
	const BYTE* const s1 = &d[size*3];
	const BYTE* const s2 = &d[size*2];
	const BYTE* const s3 = &d[size*1];
	const BYTE* const s4 = &d[size*0];

	for(U32 i = 0; i < size; ++i){
		const U32 x = s4[i] << 24 | s3[i] << 16 | s2[i] << 8 | s1[i];
		dest[i] = x;
	}

	return true;
}

class ImportWriter {
private:
	typedef Index::IndexEntry index_entry_type;
	typedef Index::IndexIndexEntry index_index_type;

public:
	ImportWriter();
	~ImportWriter();

	bool Open(const std::string output);
	bool WriteHeader(void);
	void WriteFinal(void);
	void CheckOutputNames(const std::string& input);

public:
	// Stream information
	std::ofstream stream;   // stream for data
	std::string filename;
	std::string basePath;
	std::string baseName;

	// Basic
	index_entry_type current_index_entry;
	std::vector<index_entry_type> index;
	std::vector<index_index_type> indexindex;
	U64 n_blocksWritten;            // number of blocks written
	U64 n_variants_written;         // number of variants written
};

} /* namespace Tomahawk */

#endif /* CORE_TOMAHAWKIMPORTWRITER_H_ */
