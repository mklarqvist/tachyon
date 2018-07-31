#ifndef CORE_BASE_HEADER_YON_TACHYONHEADER_H_
#define CORE_BASE_HEADER_YON_TACHYONHEADER_H_

#include "header_contig.h"
#include "header_magic.h"
#include "header_map_entry.h"
#include "header_sample.h"
#include "support/type_definitions.h"
#include "support/magic_constants.h"
#include "algorithm/OpenHashTable.h"
#include "io/basic_buffer.h"

namespace tachyon{
namespace core{

/**<
 * This class describes data that is mandatory in the Tachyon
 * file-format
 */
class VariantHeader{
private:
	typedef VariantHeader                     self_type;
	typedef HeaderMagic                       magic_type;
	typedef core::HeaderContig                contig_type;
	typedef core::HeaderMapEntry              map_entry_type;
	typedef core::HeaderSample                sample_type;
	typedef hash::HashTable<std::string, S32> hash_table_type;

public:
	explicit VariantHeader(void);
	VariantHeader(const self_type& other);
	~VariantHeader();

	inline const contig_type& getContig(const U32& position) const{ return(this->contigs[position]); }
	inline const sample_type& getSample(const U32& position) const{ return(this->samples[position]); }

	bool getContig(const std::string& p, contig_type*& target) const;
	bool getSample(const std::string& p, sample_type*& target) const;
	bool getInfoField(const std::string& p, map_entry_type*& target) const;
	bool getFormatField(const std::string& p, map_entry_type*& target) const;
	bool getFilterField(const std::string& p, map_entry_type*& target) const;
	const map_entry_type* getInfoField(const std::string& p) const;
	const map_entry_type* getFormatField(const std::string& p) const;
	const map_entry_type* getFilterField(const std::string& p) const;

	inline const U64& getSampleNumber(void) const{ return(this->header_magic.n_samples); }
	inline U64& getSampleNumber(void){ return(this->header_magic.n_samples); }
	inline const U32& getContigNumber(void) const{ return(this->header_magic.n_contigs); }
	inline U32& getContigNumber(void){ return(this->header_magic.n_contigs); }

	// write
	std::ostream& write(std::ostream& stream);

	bool has_format_field(const std::string& field_name) const{
		map_entry_type* match = nullptr;
		if(this->getFormatField(field_name, match))
			return true;

		return false;
	}

	bool has_info_field(const std::string& field_name) const{
		map_entry_type* match = nullptr;
		if(this->getInfoField(field_name, match))
			return true;

		return false;
	}

	bool has_filter_field(const std::string& field_name) const{
		map_entry_type* match = nullptr;
		if(this->getFilterField(field_name, match))
			return true;

		return(false);
	}

	std::ostream& writeHeaderVCF(std::ostream& stream, const bool showFormat = true) const;

private:
	bool buildHashTables(void);

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const self_type& header){
		buffer << header.header_magic;
		for(U32 i = 0; i < header.header_magic.n_contigs; ++i) buffer << header.contigs[i];
		for(U32 i = 0; i < header.header_magic.n_samples; ++i) buffer << header.samples[i];
		for(U32 i = 0; i < header.header_magic.n_info_values; ++i)   buffer << header.info_fields[i];
		for(U32 i = 0; i < header.header_magic.n_format_values; ++i) buffer << header.format_fields[i];
		for(U32 i = 0; i < header.header_magic.n_filter_values; ++i) buffer << header.filter_fields[i];

		buffer.Add(&header.literals[0], header.literals.size());
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& header){
		buffer >> header.header_magic;
		delete [] header.contigs;
		delete [] header.samples;
		delete [] header.info_fields;
		delete [] header.format_fields;
		delete [] header.filter_fields;
		header.contigs = new contig_type[header.header_magic.n_contigs];
		header.samples = new sample_type[header.header_magic.n_samples];
		header.info_fields   = new map_entry_type[header.header_magic.n_info_values];
		header.format_fields = new map_entry_type[header.header_magic.n_format_values];
		header.filter_fields = new map_entry_type[header.header_magic.n_filter_values];
		for(U32 i = 0; i < header.header_magic.n_contigs; ++i) buffer >> header.contigs[i];
		for(U32 i = 0; i < header.header_magic.n_samples; ++i) buffer >> header.samples[i];
		for(U32 i = 0; i < header.header_magic.n_info_values; ++i)   buffer >> header.info_fields[i];
		for(U32 i = 0; i < header.header_magic.n_format_values; ++i) buffer >> header.format_fields[i];
		for(U32 i = 0; i < header.header_magic.n_filter_values; ++i) buffer >> header.filter_fields[i];

		header.literals.resize(header.header_magic.l_literals);
		buffer.read(&header.literals[0], header.header_magic.l_literals);

		header.buildHashTables();

		return(buffer);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.header_magic;

		delete [] entry.contigs;
		entry.contigs = new contig_type[entry.header_magic.n_contigs];
		for(U32 i = 0; i < entry.header_magic.n_contigs; ++i) stream >> entry.contigs[i];

		delete [] entry.samples;
		entry.samples = new sample_type[entry.header_magic.n_samples];
		for(U32 i = 0; i < entry.header_magic.n_samples; ++i) stream >> entry.samples[i];

		delete [] entry.info_fields;
		entry.info_fields = new map_entry_type[entry.header_magic.n_info_values];
		for(U32 i = 0; i < entry.header_magic.n_info_values; ++i) stream >> entry.info_fields[i];

		delete [] entry.format_fields;
		entry.format_fields = new map_entry_type[entry.header_magic.n_format_values];
		for(U32 i = 0; i < entry.header_magic.n_format_values; ++i) stream >> entry.format_fields[i];

		delete entry.filter_fields;
		entry.filter_fields = new map_entry_type[entry.header_magic.n_filter_values];
		for(U32 i = 0; i < entry.header_magic.n_filter_values; ++i) stream >> entry.filter_fields[i];

		entry.literals.resize(entry.header_magic.l_literals);
		stream.read(&entry.literals[0], entry.header_magic.l_literals);

		entry.buildHashTables();

		return(stream);
	}

public:
	magic_type       header_magic;
	std::string      literals;
	contig_type*     contigs;
	sample_type*     samples;
	map_entry_type*  info_fields;
	map_entry_type*  format_fields;
	map_entry_type*  filter_fields;

	// Constructed during run-time
	hash_table_type* htable_contigs; // hash table for contig names
	hash_table_type* htable_samples; // hash table for sample names
	hash_table_type* htable_info_fields; // hash map from name to identifier
	hash_table_type* htable_format_fields; // hash map from name to identifier
	hash_table_type* htable_filter_fields; // hash map from name to identifier
};

}
}



#endif /* CORE_BASE_HEADER_YON_TACHYONHEADER_H_ */
