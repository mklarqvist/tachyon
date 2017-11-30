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
	explicit Header(void) : n_contigs(0), n_samples(0), n_entries(0), contigs(nullptr), samples(nullptr), entries(nullptr), htable_contigs(nullptr), htable_samples(nullptr), htable_entries(nullptr){}
	~Header(){
		delete [] this->contigs;
		delete [] this->samples;
		delete [] this->entries;
		delete this->htable_contigs;
		delete this->htable_samples;
		delete this->htable_entries;
	}

	inline const contig_type& getContig(const U32& p) const{ return(this->contigs[p]); }
	inline const sample_type& getSample(const U32& p) const{ return(this->samples[p]); }
	inline const map_entry_type& getEntry(const U32& p) const{ return(this->entries[p]); }

private:
	friend std::ifstream& operator<<(std::ifstream& stream, self_type& entry){
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

		return(stream);
	}

public:
	U32 n_contigs;
	U32 n_samples;
	U32 n_entries;
	contig_type* contigs;
	sample_type* samples;
	map_entry_type* entries;

	// Constructed during run-time
	hash_table_type* htable_contigs; // hash table for contig names
	hash_table_type* htable_samples; // hash table for sample names
	hash_table_type* htable_entries; // hash map from name to identifier
};

}
}



#endif /* CORE_BASE_HEADER_HEADER_H_ */
