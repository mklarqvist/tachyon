#ifndef CORE_BASE_HEADER_YON_TACHYONHEADER_H_
#define CORE_BASE_HEADER_YON_TACHYONHEADER_H_

#include "../../support/type_definitions.h"
#include "../../support/MagicConstants.h"
#include "../../algorithm/OpenHashTable.h"
#include "header_contig.h"
#include "header_map_entry.h"
#include "header_sample.h"
#include "header_magic.h"
#include "../../io/vcf/VCFHeader.h"

namespace tachyon{
namespace core{

/**<
 * This class describes data mandatory data in the
 * Tachyon header
 */
class TachyonHeader{
private:
	typedef TachyonHeader                     self_type;
	typedef HeaderMagic                       magic_type;
	typedef core::HeaderContig                contig_type;
	typedef core::HeaderMapEntry              map_entry_type;
	typedef core::HeaderSample                sample_type;
	typedef hash::HashTable<std::string, S32> hash_table_type;
	typedef vcf::VCFHeader                    vcf_header_type;

public:
	explicit TachyonHeader(void);
	TachyonHeader(const vcf_header_type& vcf_header);
	~TachyonHeader();

	inline const contig_type& getContig(const U32& position) const{ return(this->contigs[position]); }
	inline const sample_type& getSample(const U32& position) const{ return(this->samples[position]); }
	inline const map_entry_type& getEntry(const U32& position) const{ return(this->entries[position]); }

	const bool getContig(const std::string& p, contig_type*& target) const;
	const bool getSample(const std::string& p, sample_type*& target) const;
	const bool getEntry(const std::string& p, map_entry_type*& target) const;
	const map_entry_type* getEntry(const std::string& p) const;

	inline const U64& getSampleNumber(void) const{ return(this->header_magic.n_samples); }
	inline U64& getSampleNumber(void){ return(this->header_magic.n_samples); }
	inline const U32& getContigNumber(void) const{ return(this->header_magic.n_contigs); }
	inline U32& getContigNumber(void){ return(this->header_magic.n_contigs); }

	/**<
	 * Interconvert a VCF header (provided during import) to
	 * a tachyon header
	 * @param vcf_header Target input VCF header
	 */
	void operator=(const vcf_header_type& vcf_header){
		this->header_magic.n_contigs      = vcf_header.contigs.size();
		this->header_magic.n_samples      = vcf_header.sampleNames.size();
		this->header_magic.n_declarations = vcf_header.map.size();

		for(U32 i = 0; i < vcf_header.literal_lines.size(); ++i)
			this->literals += vcf_header.literal_lines[i];

		this->header_magic.l_literals     = this->literals.size();

		// Cleanup previous
		delete [] this->contigs;
		delete [] this->samples;
		delete [] this->entries;

		this->contigs = new contig_type[this->header_magic.getNumberContigs()];
		for(U32 i = 0; i < this->header_magic.getNumberContigs(); ++i)
			this->contigs[i] = vcf_header.contigs[i];

		this->samples = new sample_type[this->header_magic.getNumberSamples()];
		for(U32 i = 0; i < this->header_magic.getNumberSamples(); ++i)
			this->samples[i] = vcf_header.sampleNames[i];

		this->entries = new map_entry_type[this->header_magic.n_declarations];
		for(U32 i = 0; i < this->header_magic.n_declarations; ++i)
			this->entries[i] = vcf_header.map[i];

		this->buildMapTable();
		this->buildHashTables();
	}

	// write
	std::ofstream& write(std::ofstream& stream){
		stream << this->header_magic;
		for(U32 i = 0; i < this->header_magic.n_contigs; ++i)
			stream << this->contigs[i];

		for(U32 i = 0; i < this->header_magic.n_samples; ++i)
			stream << this->samples[i];

		for(U32 i = 0; i < this->header_magic.n_declarations; ++i)
			stream << this->entries[i];

		stream.write(&this->literals[0], this->literals.size());
		return(stream);
	}

	const bool has_format_field(const std::string& field_name) const{
		map_entry_type* match = nullptr;
		if(this->getEntry(field_name, match))
			return true;

		return false;
	}

	const bool has_info_field(const std::string& field_name) const{
		map_entry_type* match = nullptr;
		if(this->getEntry(field_name, match))
			return true;

		return false;
	}

	const bool has_filter_field(const std::string& field_name) const{
		map_entry_type* match = nullptr;
		if(this->getEntry(field_name, match))
			return true;

		return(false);
	}

private:
	bool buildMapTable(void);
	bool buildHashTables(void);

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.header_magic;

		entry.contigs = new contig_type[entry.header_magic.n_contigs];
		for(U32 i = 0; i < entry.header_magic.n_contigs; ++i)
			stream >> entry.contigs[i];

		entry.samples = new sample_type[entry.header_magic.n_samples];
		for(U32 i = 0; i < entry.header_magic.n_samples; ++i)
			stream >> entry.samples[i];

		entry.entries = new map_entry_type[entry.header_magic.n_declarations];
		for(U32 i = 0; i < entry.header_magic.n_declarations; ++i)
			stream >> entry.entries[i];

		entry.literals.resize(entry.header_magic.l_literals);
		stream.read(&entry.literals[0], entry.header_magic.l_literals);

		entry.buildMapTable();
		entry.buildHashTables();

		return(stream);
	}

public:
	magic_type       header_magic;
	std::string      literals;
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
