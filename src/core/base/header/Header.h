#ifndef CORE_BASE_HEADER_HEADER_H_
#define CORE_BASE_HEADER_HEADER_H_

#include "../../../support/TypeDefinitions.h"
#include "../../../support/MagicConstants.h"
#include "HeaderContig.h"
#include "HeaderMapEntry.h"
#include "HeaderSample.h"
#include "../../../algorithm/OpenHashTable.h"

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

	inline const bool getContig(const std::string& p, contig_type*& target) const{
		if(this->htable_contigs == nullptr) return false;
		S32* ret = nullptr;
		if(this->htable_contigs->GetItem(&p[0], &p, ret, p.size())){
			target = &this->contigs[*ret];
			return true;
		}
		return false;
	}

	inline const bool getSample(const std::string& p, sample_type*& target) const{
		if(this->htable_samples == nullptr) return false;
		S32* ret = nullptr;
		if(this->htable_samples->GetItem(&p[0], &p, ret, p.size())){
			target = &this->samples[*ret];
			return true;
		}
		return false;
	}

	inline const bool getEntry(const std::string& p, map_entry_type*& target) const{
		if(this->htable_entries == nullptr) return false;
		S32* ret = nullptr;
		if(this->htable_entries->GetItem(&p[0], &p, ret, p.size())){
			target = &this->entries[*ret];
			return true;
		}
		return false;
	}

private:
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
			//std::cerr << i << "->" << this->mapTable[i] << std::endl;
		}

		return true;
	}

	bool buildHashTables(void){
		if(this->n_contigs){
			if(this->n_contigs*2 < 5012){
				this->htable_contigs = new hash_table_type(5012);
			} else
				this->htable_contigs = new hash_table_type(this->n_contigs*2);

			for(S32 i = 0; i < this->n_contigs; ++i){
				this->htable_contigs->SetItem(&this->contigs[i].name[0], &this->contigs[i].name, i, this->contigs[i].name.size());
			}
		}

		if(this->n_samples){
			if(this->n_samples*2 < 5012){
				this->htable_samples = new hash_table_type(5012);
			} else
				this->htable_samples = new hash_table_type(this->n_samples*2);

			for(S32 i = 0; i < this->n_samples; ++i){
				this->htable_samples->SetItem(&this->samples[i].name[0], &this->samples[i].name, i, this->samples[i].name.size());
			}
		}

		if(this->n_entries){
			if(this->n_entries*2 < 5012){
				this->htable_entries = new hash_table_type(5012);
			} else
				this->htable_entries = new hash_table_type(this->n_entries*2);

			for(S32 i = 0; i < this->n_entries; ++i){
				this->htable_entries->SetItem(&this->entries[i].ID[0], &this->entries[i].ID, i, this->entries[i].ID.size());
			}
		}

		return true;
	}

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
		entry.buildHashTables();

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
