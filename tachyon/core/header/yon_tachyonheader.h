#ifndef CORE_BASE_HEADER_YON_TACHYONHEADER_H_
#define CORE_BASE_HEADER_YON_TACHYONHEADER_H_

#include "../../support/type_definitions.h"
#include "../../support/MagicConstants.h"
#include "../../algorithm/OpenHashTable.h"
#include "header_contig.h"
#include "header_map_entry.h"
#include "header_sample.h"

namespace tachyon{
namespace core{

class TachyonHeader{
private:
	typedef TachyonHeader                     self_type;
	typedef core::HeaderContig                contig_type;
	typedef core::HeaderMapEntry              map_entry_type;
	typedef core::HeaderSample                sample_type;
	typedef hash::HashTable<std::string, S32> hash_table_type;

public:
	explicit TachyonHeader(void);
	~TachyonHeader();

	inline const contig_type& getContig(const U32& p) const{ return(this->contigs[p]); }
	inline const sample_type& getSample(const U32& p) const{ return(this->samples[p]); }
	inline const map_entry_type& getEntry(const U32& p) const{ return(this->entries[p]); }

	const bool getContig(const std::string& p, contig_type*& target) const;
	const bool getSample(const std::string& p, sample_type*& target) const;
	const bool getEntry(const std::string& p, map_entry_type*& target) const;

private:
	bool buildMapTable(void);
	bool buildHashTables(void);

	friend std::ifstream& operator<<(std::ifstream& stream, self_type& entry){
		entry.file_header_string.resize(constants::FILE_HEADER.size());
		stream.read(&entry.file_header_string[0], constants::FILE_HEADER.size());
		stream.read(reinterpret_cast<char*>(&entry.version_major),sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.version_minor),sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.version_patch),sizeof(U16));

		stream.read(reinterpret_cast<char*>(&entry.n_contigs),sizeof(U32));
		entry.contigs = new contig_type[entry.n_contigs];
		for(U32 i = 0; i < entry.n_contigs; ++i)
			stream >> entry.contigs[i];

		stream.read(reinterpret_cast<char*>(&entry.n_samples),sizeof(U32));
		entry.samples = new sample_type[entry.n_samples];
		for(U32 i = 0; i < entry.n_samples; ++i)
			stream >> entry.samples[i];

		stream.read(reinterpret_cast<char*>(&entry.n_entries),sizeof(U32));
		entry.entries = new map_entry_type[entry.n_entries];
		for(U32 i = 0; i < entry.n_entries; ++i)
			stream >> entry.entries[i];

		entry.buildMapTable();
		entry.buildHashTables();

		return(stream);
	}

public:
	std::string      file_header_string;
	U16              version_major;
	U16              version_minor;
	U16              version_patch;
	U32              n_contigs;
	U32              n_samples;
	U32              n_entries;
	contig_type*     contigs;
	sample_type*     samples;
	map_entry_type*  entries;
	U32*             mapTable;

	// Constructed during run-time
	hash_table_type* htable_contigs; // hash table for contig names
	hash_table_type* htable_samples; // hash table for sample names
	hash_table_type* htable_entries; // hash map from name to identifier
};

}
}



#endif /* CORE_BASE_HEADER_YON_TACHYONHEADER_H_ */
