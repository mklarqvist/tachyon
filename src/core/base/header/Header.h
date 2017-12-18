#ifndef CORE_BASE_HEADER_HEADER_H_
#define CORE_BASE_HEADER_HEADER_H_

namespace Tachyon{
namespace Core{

class Header{
private:
	typedef Header self_type;
	typedef Core::HeaderContig contig_type;
	typedef Core::HeaderMapEntry map_entry_type;
	typedef Core::HeaderSample sample_type;
	typedef Hash::HashTable<std::string, S32> hash_table_type;

public:
	explicit Header(void) :
		version_major(0),
		version_minor(0),
		version_patch(0),
		n_contigs(0),
		n_samples(0),
		n_entries(0),
		contigs(nullptr),
		samples(nullptr),
		entries(nullptr),
		mapTable(nullptr),
		htable_contigs(nullptr),
		htable_samples(nullptr),
		htable_entries(nullptr)
	{}

	~Header(){
		delete [] this->contigs;
		delete [] this->samples;
		delete [] this->entries;
		delete [] this->mapTable;
		delete this->htable_contigs;
		delete this->htable_samples;
		delete this->htable_entries;
	}

	inline const contig_type& getContig(const U32& p) const{ return(this->contigs[p]); }
	inline const sample_type& getSample(const U32& p) const{ return(this->samples[p]); }
	inline const map_entry_type& getEntry(const U32& p) const{ return(this->entries[p]); }

	bool buildMapTable(void){
		if(this->n_entries == 0)
			return false;

		delete [] this->mapTable;

		S32 largest_idx = -1;
		for(U32 i = 0; i < this->n_entries; ++i){
			if(this->entries[i].IDX > largest_idx)
				largest_idx = this->entries[i].IDX;
		}

		this->mapTable = new U32[largest_idx + 1];
		memset(this->mapTable, 0, sizeof(U32)*(largest_idx+1));
		S32 localID = 0;
		for(U32 i = 0; i < this->n_entries; ++i){
			this->mapTable[this->entries[i].IDX] = localID++;
			std::cerr << i << "->" << this->mapTable[i] << std::endl;
		}

		return true;
	}

private:
	friend std::ifstream& operator<<(std::ifstream& stream, self_type& entry){
		entry.file_header_string.resize(Constants::FILE_HEADER.size());
		stream.read(&entry.file_header_string[0], Constants::FILE_HEADER.size());
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

		return(stream);
	}

public:
	std::string file_header_string;
	U16 version_major;
	U16 version_minor;
	U16 version_patch;
	U32 n_contigs;
	U32 n_samples;
	U32 n_entries;
	contig_type* contigs;
	sample_type* samples;
	map_entry_type* entries;
	U32* mapTable;

	// Constructed during run-time
	hash_table_type* htable_contigs; // hash table for contig names
	hash_table_type* htable_samples; // hash table for sample names
	hash_table_type* htable_entries; // hash map from name to identifier
};

}
}



#endif /* CORE_BASE_HEADER_HEADER_H_ */
