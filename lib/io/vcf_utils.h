#ifndef IO_HTSLIB_INTEGRATION_H_
#define IO_HTSLIB_INTEGRATION_H_

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

#include "htslib/kstring.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"
#include "htslib/hts.h"

#include "support/helpers.h"
#include "io/basic_buffer.h"

namespace tachyon {
namespace io {

const std::vector<std::string> BCF_TYPE_LOOKUP = {"NULL","INT8","INT16","INT32",
                                                  "ERROR","FLOAT","ERROR","CHAR"};

struct VcfContig {
public:
	VcfContig() : idx(0), n_bases(0){}
	~VcfContig() = default;

public:
	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The name of the contig. Canonically this is the first
	// non-whitespace-containing string after the > marker in a FASTA file.
	// For example, the line:
	//      >chr1 more info here
	// has a name of "chr1" and a description of "more info here"
	std::string name;

	// Ideally this record is filled in as described above, but not all FASTA
	// readers capture the description information after the name. Since a
	// description is not required by the FASTA spec, we cannot distinguish cases
	// where a description was not present and where a parser ignored it.
	std::string description;

	// The length of this contig in basepairs.
	int64_t n_bases;

	// Additional information used when reading and writing VCF headers. An
	// example map of key-value extra fields would transform an input line
	// containing 'assembly=B36,taxonomy=x,species="Homo sapiens"' to a map with
	// "assembly" -> "B36", "taxonomy" -> "x", "species" -> "Homo sapiens". We
	// never use this information internally, other than reading it in so we can
	// write the contig out again.
	std::vector< std::pair<std::string, std::string> > extra;
};

// Temp declare
struct VcfInfo{
public:
	VcfInfo() : idx(0){}
	~VcfInfo() = default;
public:
	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The unique ID of the INFO field. Examples include "MQ0" or "END".
	std::string id;

	// Required. The number of values included with the info field. This should be
	// the string representation of the number, e.g. "1" for a single entry, "2"
	// for a pair of entries, etc. Special cases arise when the number of entries
	// depend on attributes of the Variant or are unknown in advance, and include:
	// "A": The field has one value per alternate allele.
	// "R": The field has one value per allele (including the reference).
	// "G": The field has one value for each possible genotype.
	// ".": The number of values varies, is unknown, or is unbounded.
	std::string number;

	// Required. The type of the INFO field. Valid values are "Integer", "Float",
	// "Flag", "Character", and "String".
	std::string type;

	// Required by VCF. The description of the field.
	std::string description;

	// Optional. The annotation source used to generate the field.
	std::string source;

	// Optional. The version of the annotation source used to generate the field.
	std::string version;
};

struct VcfFormat{
public:
	VcfFormat() : idx(0){}
	~VcfFormat() = default;

public:
	// Required. The unique ID of the FORMAT field. Examples include "GT", "PL".
	std::string id;

	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The number of entries expected. See description above in the
	// VcfInfo message.
	std::string number;

	// Required. The type of the field. Valid values are "Integer", "Float",
	// "Character", and "String" (same as INFO except "Flag" is not supported).
	std::string type;

	// Required by VCF. The description of the field.
	std::string description;
};

struct VcfFilter{
public:
	VcfFilter() : idx(0){}
	~VcfFilter() = default;

	friend std::ostream& operator<<(std::ostream& stream, const VcfFilter& flt){
		stream.write((const char*)&flt.idx, sizeof(uint32_t));

		utility::SerializeString(flt.id, stream);
		utility::SerializeString(flt.description, stream);

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, VcfFilter& flt){
		stream.read((char*)&flt.idx, sizeof(uint32_t));

		utility::DeserializeString(flt.id, stream);
		utility::DeserializeString(flt.description, stream);

		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VcfFilter& flt){
		io::SerializePrimitive(flt.idx, buffer);
		io::SerializeString(flt.id, buffer);
		io::SerializeString(flt.description, buffer);
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VcfFilter& flt){
		io::DeserializePrimitive(flt.idx, buffer);
		io::DeserializeString(flt.id, buffer);
		io::DeserializeString(flt.description, buffer);
		return(buffer);
	}

public:
	// Required. The internal identifier for this field
	uint32_t idx;

	// Required. The unique ID of the filter. Examples include "PASS", "RefCall".
	std::string id;

	// Required by VCF. The description of the filter.
	std::string description;
};

// This record type is a catch-all for other types of headers. For example,
// ##pedigreeDB=http://url_of_pedigrees
// The VcfExtra message would represent this with key="pedigreeDB",
// value="http://url_of_pedigrees".
struct VcfExtra{
public:
	VcfExtra() = default;
	VcfExtra(const std::string& key, const std::string& value) :
		key(key),
		value(value)
	{}

	~VcfExtra() = default;

	friend std::ostream& operator<<(std::ostream& stream, const VcfExtra& extra){
		utility::SerializeString(extra.key, stream);
		utility::SerializeString(extra.value, stream);
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, VcfExtra& extra){
		utility::DeserializeString(extra.key, stream);
		utility::DeserializeString(extra.value, stream);
		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VcfExtra& extra){
		io::SerializeString(extra.key, buffer);
		io::SerializeString(extra.value, buffer);
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VcfExtra& extra){
		io::DeserializeString(extra.key, buffer);
		io::DeserializeString(extra.value, buffer);
		return(buffer);
	}

public:
  // Required by VCF. The key of the extra header field. Note that this key does
  // not have to be unique within a VcfHeader.
  std::string key;

  // Required by VCF. The value of the extra header field.
  std::string value;
};

// This record type is a catch-all for other headers containing multiple
// key-value pairs. For example, headers may have META lines that provide
// metadata about the VCF as a whole, e.g.
// ##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
// The VcfStructuredExtra message would represent this with key="META",
// and fields mapping "ID" -> "Assay", "Type" -> "String", etc.
struct VcfStructuredExtra{
public:
	VcfStructuredExtra() = default;
	~VcfStructuredExtra() = default;

	friend std::ostream& operator<<(std::ostream& stream, const VcfStructuredExtra& extra){
		utility::SerializeString(extra.key, stream);
		size_t l_extra = extra.fields.size();
		stream.write((const char*)&l_extra, sizeof(size_t));
		for(U32 i = 0; i < extra.fields.size(); ++i)
			stream << extra.fields[i];

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, VcfStructuredExtra& extra){
		utility::DeserializeString(extra.key, stream);
		size_t l_extra;
		stream.read((char*)&l_extra, sizeof(size_t));
		extra.fields.resize(l_extra);
		for(U32 i = 0; i < extra.fields.size(); ++i)
			stream >> extra.fields[i];

		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VcfStructuredExtra& extra){
		io::SerializeString(extra.key, buffer);
		size_t l_extra = extra.fields.size();
		io::SerializePrimitive(l_extra, buffer);
		for(U32 i = 0; i < extra.fields.size(); ++i)
			buffer << extra.fields[i];

		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VcfStructuredExtra& extra){
		io::DeserializeString(extra.key, buffer);
		size_t l_extra;
		io::DeserializePrimitive(l_extra, buffer);
		extra.fields.resize(l_extra);
		for(U32 i = 0; i < extra.fields.size(); ++i)
			buffer >> extra.fields[i];

		return(buffer);
	}

public:
	// Required by VCF. The key of the extra header field. Note that this key does
	// not have to be unique within a VcfHeader.
	std::string key;

	// Required by VCF. The key=value pairs contained in the structure.
	std::vector<VcfExtra> fields;
};

//
// -----------------------------------------------------------------------------
// VCF type encoding utilities
template<class T>
struct VcfType {
  // Predicates for checking missing and sentinel entries.  Use these, not ==.
  // Is argument the "missing" value?
  static bool IsMissing(T);
  // Is argument the vector end sentinel value?
  static bool IsVectorEnd(T);
};

// See interface description comment above.
template<>
struct VcfType<int8_t> {
  static bool IsMissing(int8_t v)    { return (v == bcf_int8_missing); }
  static bool IsVectorEnd(int8_t v)  { return (v == bcf_int8_vector_end); }
};

// See interface description comment above.
template<>
struct VcfType<int16_t> {
  static bool IsMissing(int16_t v)    { return (v == bcf_int16_missing); }
  static bool IsVectorEnd(int16_t v)  { return (v == bcf_int16_vector_end); }
};

// See interface description comment above.
template<>
struct VcfType<int> {
  static bool IsMissing(int v)    { return (v == bcf_int32_missing); }
  static bool IsVectorEnd(int v)  { return (v == bcf_int32_vector_end); }
};

// See interface description comment above.
template<>
struct VcfType<float> {
  static bool IsMissing(float v)    { return bcf_float_is_missing(v); }
  static bool IsVectorEnd(float v)  { return bcf_float_is_vector_end(v); }
};


template <class T>
struct VcfGenotype {
	// Predicates for checking missing and sentinel entries.
	static bool IsMissing(const T& value){ return(value == bcf_gt_missing); }
};

// Genotype helper
struct VcfGenotypeSummary{
public:
	VcfGenotypeSummary(void) :
		base_ploidy(0),
		phase_if_uniform(0),
		mixed_phasing(false),
		invariant(false),
		n_missing(0),
		n_vector_end(0)
	{}

	~VcfGenotypeSummary() = default;

	/**<
	 * Gathers summary statistics for a vector of genotypes
	 * at a given site. Collects information regarding the
	 * number of missing genotypes and count of sentinel
	 * nodes, checks if the phasing is uniform and whether
	 * all the genotypes are identical.
	 * Todo: Use hashes to check for uniformity of genotypes.
	 * @param n_samples Total number of samples in the input vector.
	 *                  This is equivalent to the samples in the file.
	 * @param fmt       The target htslib format structure.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	template<class T> bool evaluate(const size_t& n_samples, const bcf_fmt_t& fmt){
		if(fmt.p_len == 0) return true;
		assert(fmt.size/fmt.n == sizeof(T));

		// Set the base ploidy. This corresponds to the LARGEST
		// ploidy observed of ANY individual sample at the given
		// locus. If a genotype has a ploidy < base ploidy then
		// it is trailed with the sentinel node symbol to signal
		// that the remainder of the vector is NULL.
		this->base_ploidy = fmt.n;

		// Find first phase
		this->phase_if_uniform = 0;
		// Iterate over genotypes to find the first valid phase
		// continuing the search if the current value is the
		// sentinel node symbol.
		int j = fmt.n - 1;
		for(uint32_t i = 0; i < n_samples; ++i){
			if(VcfGenotype<T>::IsMissing(fmt.p[j]) == true
			   || VcfType<T>::IsVectorEnd(fmt.p[j]) == true)
				j += fmt.n;
			else {
				this->phase_if_uniform = fmt.p[j] & 1;
				break;
			}
		}

		// Iterate over genotypes to compute summary statistics
		// regarding missingness, number of special sentinel
		// symbols and assess uniformity of phasing.
		j = fmt.n - 1;
		for(uint32_t i = 0; i < n_samples; ++i){
			if(VcfGenotype<T>::IsMissing(fmt.p[j]) == false
			   && VcfType<int8_t>::IsVectorEnd(fmt.p[j]) == false
			   && (fmt.p[j] & 1) != this->phase_if_uniform)
			{
				this->mixed_phasing = true;
			}

			// Iterate over the number of chromosomes / individual
			for(int k = 0; k < fmt.n; ++k, ++j){
				this->n_missing    += VcfGenotype<T>::IsMissing(fmt.p[j]);
				this->n_vector_end += VcfType<T>::IsVectorEnd(fmt.p[j]);
			}
		}

		return true;
	}

	inline bool isBaseDiploid(void) const{ return(this->base_ploidy == 2); }

public:
	uint8_t  base_ploidy;
	bool     phase_if_uniform;
	bool     mixed_phasing;
	bool     invariant;
	uint64_t n_missing;
	uint64_t n_vector_end;
};

// Returns the hrec that contains information or nullptr if none does.
const bcf_hrec_t* GetPopulatedHrec(const bcf_idpair_t& idPair);

class VcfHeader{
public:
	typedef VcfHeader self_type;
	typedef VcfContig contig_type;
	typedef bcf_hdr_t hts_vcf_header;
	typedef VcfFormat format_type;
	typedef VcfInfo   info_type;
	typedef VcfFilter filter_type;
	typedef VcfStructuredExtra structured_extra_type;
	typedef VcfExtra  extra_type;
	typedef std::unordered_map<std::string, uint32_t> map_type;
	typedef std::unordered_map<uint32_t, uint32_t>    map_reverse_type;

public:
	VcfHeader() = default;
	VcfHeader(const VcfHeader& other) :
		fileformat_string_(other.fileformat_string_),
		literals_(other.literals_),
		samples_(other.samples_),
		contigs_(other.contigs_),
		info_fields_(other.info_fields_),
		format_fields_(other.format_fields_),
		filter_fields_(other.filter_fields_),
		structured_extra_fields_(other.structured_extra_fields_),
		extra_fields_(other.extra_fields_)
	{
		this->BuildMaps();
		this->BuildReverseMaps();
	}

	~VcfHeader() = default;

	inline size_t GetNumberSamples(void) const{ return(this->samples_.size()); }
	inline size_t GetNumberContigs(void) const{ return(this->contigs_.size()); }

	// Adds Contig information from the idPair to the ContigInfo object.
	void AddContigInfo(const bcf_idpair_t& idPair) {
	  // ID and length are special-cased in the idPair.
		//std::cerr << "Contig: " << pos_in_fasta << "\t" << std::string(idPair.key) << ": " << idPair.val->info[0] << std::endl;
		VcfContig c;
		c.name = idPair.key;
		c.n_bases = idPair.val->info[0];

	  const bcf_hrec_t* hrec0 = GetPopulatedHrec(idPair);
	  if (hrec0 != nullptr) {
	    for (int j = 0; j < hrec0->nkeys; j++) {
	    	const std::string current_key(hrec0->keys[j]);
			// Add any non-ID and non-length info to the structured map of additional
			// information.
	    	if (current_key == "ID" ||
	    	    current_key == "length")
	    	{
	    		//continue;
	    	} else if(current_key == "IDX"){
	    		c.idx = atoi(hrec0->vals[j]);
	    	} else {
	    		c.extra.push_back(std::pair<std::string,std::string>(current_key, std::string(hrec0->vals[j])));
	    	}
	    }
	  } else {
		  std::cerr << "hrec error" << std::endl;
		  return;
	  }

	  // Add current contig to map
	  if(this->contigs_.size() == 0){
		  this->contigs_.push_back(c);
		  this->contigs_map_[c.name] = 0;
		  return;
	  }

	  if(this->contigs_map_.find(c.name) == this->contigs_map_.end()){
		  this->contigs_map_[c.name] = this->contigs_.size();
		  this->contigs_.push_back(c);
	  } else {
		  std::cerr << "illegal: duplicated contig name" << std::endl;
	  }
	}

	// Adds FILTER information from the bcf_hrec_t to the VcfFilterInfo object.
	void AddFilterInfo(const bcf_hrec_t* hrec) {
	  if (hrec->nkeys >= 2 && std::string(hrec->keys[0]) == "ID" &&
	      std::string(hrec->keys[1]) == "Description")
	  {
		VcfFilter f;
	    f.id = std::string(hrec->vals[0]);
	    f.description = std::string(hrec->vals[1]);
	    for(int i = 2; i < hrec->nkeys; ++i){
	    	if(std::string(hrec->keys[i]) == "IDX"){
	    		f.idx = atoi(hrec->vals[i]);
	    	}
	    }

	    // Add current filter field to map.
	    if(this->filter_fields_.size() == 0){
			this->filter_fields_.push_back(f);
			this->filter_fields_map_[f.id] = 0;
			return;
		}

		if(this->filter_fields_map_.find(f.id) == this->filter_fields_map_.end()){
			this->filter_fields_map_[f.id] = this->filter_fields_.size();
			this->filter_fields_.push_back(f);
		} else {
			std::cerr << "illegal: duplicated filter name: " << f.id << std::endl;
		}

	  } else {
	    std::cerr << "Malformed FILTER field detected in header, leaving this "
	                 "filter empty" << std::endl;
	  }
	}

	// Adds INFO information from the bcf_hrec_t to the VcfInfo object.
	void AddInfo(const bcf_hrec_t* hrec) {
	  if (hrec->nkeys >= 4 && std::string(hrec->keys[0]) == "ID" &&
	      std::string(hrec->keys[1]) == "Number" && std::string(hrec->keys[2]) == "Type" &&
	      std::string(hrec->keys[3]) == "Description")
	  {
		VcfInfo f;
	    f.id = std::string(hrec->vals[0]);
	    f.number = std::string(hrec->vals[1]);
	    f.type = std::string(hrec->vals[2]);
	    f.description = std::string(hrec->vals[3]);
	    for (int i = 4; i < hrec->nkeys; i++) {
	      const std::string current_key = std::string(hrec->keys[i]);
	      if (current_key == "Source") {
	        f.source = std::string(hrec->vals[i]);
	      } else if (current_key == "Version") {
	        f.version = std::string(hrec->vals[i]);
	      } else if (current_key == "IDX") {
	    	  f.idx = atoi(hrec->vals[i]);
	      }
	    }

	    // Add current info field to map.
	    if(this->info_fields_.size() == 0){
			this->info_fields_.push_back(f);
			this->info_fields_map_[f.id] = 0;
			return;
		}

		if(this->info_fields_map_.find(f.id) == this->info_fields_map_.end()){
			this->info_fields_map_[f.id] = this->info_fields_.size();
			this->info_fields_.push_back(f);
		} else {
			std::cerr << "illegal: duplicated info name: " << f.id << std::endl;
		}

	  } else {
	    std::cerr << "Malformed INFO field detected in header, leaving this "
	                 "info empty" << std::endl;
	  }
	}

	// Adds FORMAT information from the bcf_hrec_t to the VcfFormatInfo object.
	void AddFormatInfo(const bcf_hrec_t* hrec) {
	  if (hrec->nkeys >= 4 && std::string(hrec->keys[0]) == "ID" &&
	      std::string(hrec->keys[1]) == "Number" && std::string(hrec->keys[2]) == "Type" &&
	      std::string(hrec->keys[3]) == "Description")
	  {
		VcfFormat f;
	    f.id          = std::string(hrec->vals[0]);
	    f.number      = std::string(hrec->vals[1]);
	    f.type        = std::string(hrec->vals[2]);
	    f.description = std::string(hrec->vals[3]);
	    for (int i = 4; i < hrec->nkeys; i++) {
	    	const std::string current_key = std::string(hrec->keys[i]);
	    	if (current_key == "IDX") {
				  f.idx = atoi(hrec->vals[i]);
			  }
	    }

	    // Add current format field to map.
	    if(this->format_fields_.size() == 0){
			this->format_fields_.push_back(f);
			this->format_fields_map_[f.id] = 0;
			return;
		}

		if(this->format_fields_map_.find(f.id) == this->format_fields_map_.end()){
			this->format_fields_map_[f.id] = this->format_fields_.size();
			this->format_fields_.push_back(f);
		} else {
			std::cerr << "illegal: duplicated format name: " << f.id << std::endl;
		}

	  } else {
	    std::cerr << "Malformed FORMAT field detected in header, leaving this "
	                    "format empty" << std::endl;
	  }
	}

	// Adds structured information from the bcf_hrec_t to the VcfStructuredExtra.
	void AddStructuredExtra(const bcf_hrec_t* hrec) {
	  VcfStructuredExtra f;
	  f.key = std::string(hrec->key);
	  for (int i = 0; i < hrec->nkeys; i++)
		  f.fields.push_back(VcfExtra(std::string(hrec->keys[i]), std::string(hrec->vals[i])));

	  this->structured_extra_fields_.push_back(f);
	}

	// Adds unstructured information from the bcf_hrec_t to the VcfExtra object.
	void AddExtra(const bcf_hrec_t* hrec) {
	  VcfExtra f;
	  f.key   = std::string(hrec->key);
	  f.value = std::string(hrec->value);
	  this->extra_fields_.push_back(f);
	}

	void AddSample(const std::string& sample_name) {
		if(this->samples_.size() == 0){
			this->samples_.push_back(sample_name);
			this->samples_map_[sample_name] = 0;
			return;
		}

		if(this->samples_map_.find(sample_name) == this->samples_map_.end()){
			this->samples_map_[sample_name] = this->samples_.size();
			this->samples_.push_back(sample_name);
		} else {
			std::cerr << "illegal: duplicated sample name: " << sample_name << std::endl;
		}
	}

	VcfContig* GetContig(const std::string& name) {
		map_type::const_iterator it = this->contigs_map_.find(name);
		if(it != this->contigs_map_.end()) return(&this->contigs_[it->second]);
		return(nullptr);
	}

	VcfContig* GetContig(const int& idx) {
		map_reverse_type::const_iterator it = this->contigs_reverse_map_.find(idx);
		if(it != this->contigs_reverse_map_.end()) return(&this->contigs_[it->second]);
		return(nullptr);
	}

	VcfInfo* GetInfo(const std::string& name) {
		map_type::const_iterator it = this->info_fields_map_.find(name);
		if(it != this->info_fields_map_.end()) return(&this->info_fields_[it->second]);
		return(nullptr);
	}

	VcfInfo* GetInfo(const int& idx) {
		map_reverse_type::const_iterator it = this->info_fields_reverse_map_.find(idx);
		if(it != this->info_fields_reverse_map_.end()) return(&this->info_fields_[it->second]);
		return(nullptr);
	}

	VcfFormat* GetFormat(const std::string& name) {
		map_type::const_iterator it = this->format_fields_map_.find(name);
		if(it != this->format_fields_map_.end()) return(&this->format_fields_[it->second]);
		return(nullptr);
	}

	VcfFormat* GetFormat(const int& idx) {
		map_reverse_type::const_iterator it = this->format_fields_reverse_map_.find(idx);
		if(it != this->format_fields_reverse_map_.end()) return(&this->format_fields_[it->second]);
		return(nullptr);
	}

	VcfFilter* GetFilter(const std::string& name) {
		map_type::const_iterator it = this->filter_fields_map_.find(name);
		if(it != this->filter_fields_map_.end()) return(&this->filter_fields_[it->second]);
		return(nullptr);
	}

	VcfFilter* GetFilter(const int& idx) {
		map_reverse_type::const_iterator it = this->filter_fields_reverse_map_.find(idx);
		if(it != this->filter_fields_reverse_map_.end()) return(&this->filter_fields_[it->second]);
		return(nullptr);
	}

	std::string* GetSample(const std::string& name) {
		map_type::const_iterator it = this->samples_map_.find(name);
		if(it != this->samples_map_.end()) return(&this->samples_[it->second]);
		return(nullptr);
	}

	bool BuildReverseMaps(void){
		this->contigs_reverse_map_.clear();
		this->info_fields_reverse_map_.clear();
		this->format_fields_reverse_map_.clear();
		this->filter_fields_reverse_map_.clear();

		for(uint32_t i = 0; i < this->contigs_.size(); ++i)       this->contigs_reverse_map_[this->contigs_[i].idx] = i;
		for(uint32_t i = 0; i < this->info_fields_.size(); ++i)   this->info_fields_reverse_map_[this->info_fields_[i].idx] = i;
		for(uint32_t i = 0; i < this->format_fields_.size(); ++i) this->format_fields_reverse_map_[this->format_fields_[i].idx] = i;
		for(uint32_t i = 0; i < this->filter_fields_.size(); ++i) this->filter_fields_reverse_map_[this->filter_fields_[i].idx] = i;

		return true;
	}

	bool BuildMaps(void){
		this->info_fields_map_.clear();
		this->format_fields_map_.clear();
		this->filter_fields_map_.clear();
		this->contigs_map_.clear();

		for(uint32_t i = 0; i < this->contigs_.size(); ++i)       this->contigs_map_[this->contigs_[i].name] = i;
		for(uint32_t i = 0; i < this->info_fields_.size(); ++i)   this->info_fields_map_[this->info_fields_[i].id] = i;
		for(uint32_t i = 0; i < this->format_fields_.size(); ++i) this->format_fields_map_[this->format_fields_[i].id] = i;
		for(uint32_t i = 0; i < this->filter_fields_.size(); ++i) this->filter_fields_map_[this->filter_fields_[i].id] = i;

		return true;
	}

	/**<
	 * Recodes the internal IDX field for contig info, INFO, FORMAT, and FILTER
	 * from any range to the range [0, 1, ..., n-1] as desired in Tachyon.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool RecodeIndices(void){
		for(uint32_t i = 0; i < this->contigs_.size(); ++i)       this->contigs_[i].idx = i;
		for(uint32_t i = 0; i < this->info_fields_.size(); ++i)   this->info_fields_[i].idx = i;
		for(uint32_t i = 0; i < this->format_fields_.size(); ++i) this->format_fields_[i].idx = i;
		for(uint32_t i = 0; i < this->filter_fields_.size(); ++i) this->filter_fields_[i].idx = i;

		if(this->BuildMaps() == false) return false;
		if(this->BuildReverseMaps() == false) return false;
		return true;
	}

	/**<
	* Converts this header object into a hts_vcf_header object from the
	* internally stored literal string. This object is required for
	* writing out VCF/BCF files.
	* @return
	*/
	hts_vcf_header* ConvertVcfHeader(void){
		std::string internal = this->literals_;
		internal += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
		if(this->samples_.size()){
			internal += "\tFORMAT\t";
			internal += this->samples_[0];
			for(size_t i = 1; i < this->samples_.size(); ++i)
				internal += "\t" + this->samples_[i];
		}
		internal += "\n";

		hts_vcf_header* hdr = bcf_hdr_init("r");
		int ret = bcf_hdr_parse(hdr, (char*)internal.c_str());
		if(ret != 0){
			std::cerr << "failed to get bcf header from literals" << std::endl;
			bcf_hdr_destroy(hdr);
			return(nullptr);
		}

		return(hdr);
	}

	// Append a string to the literal string
	inline void AppendLiteralString(const std::string& literal_addition){ this->literals_ += literal_addition; }

public:
	// VCF file version string.
	std::string fileformat_string_;
	// Literal string for VcfHeader data. Contains all of the Vcf header data up
	// to the start of the main header line ("#CHROM"...). As such, sample names
	// are not available in this string and needs to be appended before converting
	// back into a htslib vcf header.
	std::string literals_;

	// Vcf header lines parse into:
	// Samples:   Individual sample names.
	// VcfContig: Information relating to the interpretation of a contig. Data
	//            include its name, length in bases, its internal index identifier
	//            and optional additional information.
	// VcfInfo:   Data specifying a given INFO field
	// VcfFormat: Data specifying a given FORMAT field
	// VcfFilter: Data specifying a given FILTER field
	// VcfStructuredExtra:
	std::vector<std::string>        samples_;
	std::vector<VcfContig>          contigs_;
	std::vector<VcfInfo>            info_fields_;
	std::vector<VcfFormat>          format_fields_;
	std::vector<VcfFilter>          filter_fields_;
	std::vector<VcfStructuredExtra> structured_extra_fields_;
	std::vector<VcfExtra>           extra_fields_;

	// Utility members
	//
	// Hash tables allowing the mapping from the unique identifier string
	// (such as contig name) to the relative index offset of that object.
	// This approach requires another layer of indirection when mapping
	// from the index to the actual target. For example:
	//
	// contigs[contigs_map_["chr20"].second] <- maps to the actual target
	//
	// The reverse maps allows the mapping from a unique IDX identifier
	// to the relative index offset of that object. As above, this requires
	// an addition indirect lookup to access the desired object. For example
	// mapping the first occuring FORMAT field to its name:
	//
	// reader->vcf_header_.format_fields_[reader->vcf_header_.format_fields_reverse_map_[container.at(0)->d.fmt[0].id]].id
	//
	// map_type hash tables permits mapping string name -> index offset
	// map_reverse_type hash tables permits mapping integer IDX -> index offset
	map_type samples_map_;
	map_type contigs_map_;
	map_type info_fields_map_;
	map_type format_fields_map_;
	map_type filter_fields_map_;
	map_reverse_type contigs_reverse_map_;       // map IDX -> index offset
	map_reverse_type info_fields_reverse_map_;   // map IDX -> index offset
	map_reverse_type format_fields_reverse_map_; // map IDX -> index offset
	map_reverse_type filter_fields_reverse_map_; // map IDX -> index offset
};

class VcfReader{
private:
	typedef VcfReader self_type;

public:
	// Singleton design pattern for retrieving a guaranteed unique pointer
	// to a VcfReader. This choice is to prevent inadvertent writes to the
	// target file as another file-handle is accessing it.
	static std::unique_ptr<self_type> FromFile(const std::string& variants_path){
		htsFile* fp = hts_open(variants_path.c_str(), "r");
		if (fp == nullptr) {
			std::cerr << "Could not open " << variants_path << std::endl;
			return nullptr;
		}

		bcf_hdr_t* header = bcf_hdr_read(fp);
		if (header == nullptr){
			std::cerr << "Couldn't parse header for " << fp->fn << std::endl;
			return nullptr;
		}

		return std::unique_ptr<self_type>(new self_type(variants_path, fp, header));
	}

	bool next(const int unpack_level = BCF_UN_ALL){
		if (bcf_read(this->fp_, this->header_, this->bcf1_) < 0) {
			if (bcf1_->errcode) {
				std::cerr << "Failed to parse VCF record: " << bcf1_->errcode << std::endl;
				return false;
			} else {
				return false;
			}
		}

		bcf_unpack(this->bcf1_, unpack_level);
		return true;
	}

	bool next(bcf1_t* bcf_entry, const int unpack_level = BCF_UN_ALL){
		if (bcf_read(this->fp_, this->header_, bcf_entry) < 0) {
			if (bcf_entry->errcode) {
				std::cerr << "Failed to parse VCF record: " << bcf1_->errcode << std::endl;
				return false;
			} else {
				return false;
			}
		}

		bcf_unpack(bcf_entry, unpack_level);
		return true;
	}

	/**<
	 * Utility function that writes the VcfHeader literals string into
	 * a target output stream. The literals string does NOT contain
	 * sample information or the column header string ("#CHROM...").
	 * @param stream Target output stream
	 */
	inline void PrintLiterals(std::ostream& stream) const{ stream << this->vcf_header_.literals_ << std::endl; }

	/**<
	 * Utility function that writes a valid VCF header output string
	 * to the target stream.
	 * @param stream Target output stream
	 */
	void PrintVcfHeader(std::ostream& stream) const{
		this->PrintLiterals(stream);
		stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
		if(this->vcf_header_.samples_.size()){
			stream << "\tFORMAT\t";
			stream << this->vcf_header_.samples_[0];
			for(size_t i = 1; i < this->vcf_header_.samples_.size(); ++i)
				stream << "\t" + this->vcf_header_.samples_[i];
		}
		stream << "\n";
	}

private:
	VcfReader(const std::string& variants_path,
              htsFile* fp,
			  bcf_hdr_t* header) :
    fp_(fp),
    header_(header),
    bcf1_(bcf_init())
{
    if (this->header_->nhrec < 1) {
        std::cerr << "Empty header, not a valid VCF." << std::endl;
        return;
    }

    // Store the file-format header string
    if (std::string(this->header_->hrec[0]->key) != "fileformat") {
        std::cerr << "Not a valid VCF, fileformat needed: " << variants_path << std::endl;
    } else {
    	this->vcf_header_.fileformat_string_ = std::string(this->header_->hrec[0]->key);
    }

    // Fill in the contig info for each contig in the VCF header. Directly
    // accesses the low-level C struct because there are no indirection
    // macros/functions by htslib API.
    // BCF_DT_CTG: offset for contig (CTG) information in BCF dictionary (DT).
    const int n_contigs = this->header_->n[BCF_DT_CTG];
    for (int i = 0; i < n_contigs; ++i) {
        const bcf_idpair_t& idPair = this->header_->id[BCF_DT_CTG][i];
        this->vcf_header_.AddContigInfo(idPair);
    }

    // Iterate through all hrecs (except the first, which was 'fileformat') to
    // populate the rest of the headers.
    for (int i = 1; i < this->header_->nhrec; i++) {
    	const bcf_hrec_t* hrec0 = this->header_->hrec[i];
		switch (hrec0->type) {
		case BCF_HL_CTG:
			// Contigs are populated above, since they store length in the
			// bcf_idinfo_t* structure.
			break;
		case BCF_HL_FLT:
			this->vcf_header_.AddFilterInfo(hrec0);
			break;
		case BCF_HL_INFO:
			this->vcf_header_.AddInfo(hrec0);
			break;
		case BCF_HL_FMT:
			this->vcf_header_.AddFormatInfo(hrec0);
			break;
		case BCF_HL_STR:
			this->vcf_header_.AddStructuredExtra(hrec0);
			break;
		case BCF_HL_GEN:
			this->vcf_header_.AddExtra(hrec0);
			break;
		default:
			std::cerr << "Unknown hrec0->type: " << hrec0->type << std::endl;
			break;
		}
    }

    // Populate samples info.
    int n_samples = bcf_hdr_nsamples(this->header_);
    for (int i = 0; i < n_samples; i++)
    	this->vcf_header_.AddSample(std::string(this->header_->samples[i]));

    this->vcf_header_.BuildReverseMaps();

    // Build literal VCF header string for storage.
    kstring_t htxt = {0,0,0};
	bcf_hdr_format(this->header_, 0, &htxt);
	while (htxt.l && htxt.s[htxt.l-1] == '\0') --htxt.l; // kill trailing zeros
	std::string temp = std::string(htxt.s, htxt.l);
	size_t pos = temp.find("#CHROM"); // search for start of column header line
	temp = temp.substr(0, pos);
	this->vcf_header_.literals_ = temp;
	free(htxt.s);
}

public:
	~VcfReader() {
		bcf_destroy(this->bcf1_);
		bcf_hdr_destroy(this->header_);
		hts_close(this->fp_);
	}

public:
	// Contextual representation of vcf header
	VcfHeader vcf_header_;

	// A pointer to the htslib file used to access the VCF data.
	htsFile * fp_;

	// A htslib header data structure obtained by parsing the header of this VCF.
	bcf_hdr_t * header_;

	// htslib representation of a parsed vcf line.
	bcf1_t* bcf1_;
};

}

namespace containers {

class VcfContainer{
public:
	typedef VcfContainer       self_type;
	typedef bcf1_t             value_type;
	typedef value_type&        reference;
	typedef const value_type&  const_reference;
	typedef value_type*        pointer;
	typedef const value_type*  const_pointer;
	typedef std::ptrdiff_t     difference_type;
	typedef std::size_t        size_type;

public:
	VcfContainer(void) :
		n_carry_over_(0),
		n_entries_(0),
		n_capacity_(500),
		entries_(new pointer[500])
	{
		for(size_type i = 0; i < this->capacity(); ++i)
			this->entries_[i] = nullptr;
	}

	VcfContainer(const size_type& start_capacity) :
		n_carry_over_(0),
		n_entries_(0),
		n_capacity_(start_capacity),
		entries_(new pointer[start_capacity])
	{
		for(size_type i = 0; i < this->capacity(); ++i)
			this->entries_[i] = nullptr;
	}

	~VcfContainer(){
		if(this->entries_ != nullptr){
			for(std::size_t i = 0; i < this->n_entries_; ++i)
				bcf_destroy(this->entries_[i]);

			::operator delete[](static_cast<void*>(this->entries_));
		}
	}

	VcfContainer(const VcfContainer& other) = delete; // Disallow copy ctor

	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline size_type sizeWithoutCarryOver(void) const{ return(this->n_entries_ - this->n_carry_over_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }
	inline pointer front(void){ return(this->entries_[0]); }
	inline const_pointer front(void) const{ return(this->entries_[0]); }
	inline pointer back(void){ return(this->entries_[this->size() == 0 ? 0 : this->size() - 1 - this->n_carry_over_]); }
	inline const_pointer back(void) const{ return(this->entries_[this->size() == 0 ? 0 : this->size() - 1 - this->n_carry_over_]); }

	inline void operator+=(const pointer entry){ this->entries_[this->n_entries_++] = entry; }
	inline pointer operator[](const uint32_t position){ return(this->entries_[position]); }
	inline const_pointer operator[](const uint32_t position) const{ return(this->entries_[position]); }
	inline pointer at(const uint32_t position){ return(this->entries_[position]); }
	inline const_pointer at(const uint32_t position) const{ return(this->entries_[position]); }

	inline pointer end(void){ return(this->entries_[this->n_entries_]); }
	inline const_pointer end(void) const{ return(this->entries_[this->n_entries_]); }

	void resize(const size_t new_size){
		if(new_size < this->capacity()){
			for(size_t i = new_size; i < this->n_entries_; ++i)
				bcf_destroy(this->entries_[i]);

			if(this->n_entries_ >= new_size) this->n_entries_ = new_size;
			return;
		}

		pointer* temp = new pointer[new_size];
		for(size_t i = 0; i < this->size(); ++i)
			temp[i] = this->entries_[i];

		delete [] this->entries_;
		this->entries_ = temp;
		this->n_capacity_ = new_size;
	}

	bool getVariants(const int32_t n_variants, const int64_t n_bases, std::unique_ptr<io::VcfReader>& reader){
		if(this->size() + n_variants >= this->capacity())
			this->resize(this->size() + n_variants + 64);

		VcfContainer::pointer bcf1_ = this->end();
		if(bcf1_ == nullptr) bcf1_  = bcf_init();
		if(reader->next(bcf1_) == false)
			return false;

		*this += bcf1_;

		int64_t first_pos    = bcf1_->pos;
		int32_t first_contig = bcf1_->rid;
		if(this->size() != 1){
			first_pos    = this->entries_[0]->pos;
			first_contig = this->entries_[0]->rid;
		}

		if(bcf1_->pos - first_pos > n_bases || first_contig != bcf1_->rid){
			this->n_carry_over_ = 1;
			return(this->size() - 1);
		}

		for(int32_t i = 1; i < n_variants; ++i){
			bcf1_ = this->end();
			if(bcf1_ == nullptr) bcf1_  = bcf_init();

			if(reader->next(bcf1_) == false)
				return(this->size());

			*this += bcf1_;

			if(bcf1_->pos - first_pos > n_bases || first_contig != bcf1_->rid){
				this->n_carry_over_ = 1;
				return(this->size() - 1);
			}
		}

		return(this->size());
	}

	// Calculate genotype summary statistics from a lazy evaluated bcf1_t struct.
	// Warning: this function does NOT check if the FORMAT field GT exists either
	// in the header or in the structure itself. The assumption is that it does
	// exist and according to the Bcf specification has to be the first FORMAT
	// field set.
	io::VcfGenotypeSummary GetGenotypeSummary(const uint32_t position, const uint64_t& n_samples) const{
		io::VcfGenotypeSummary g;

		// If there are no FORMAT fields there cannot exist any
		// GT data.
		if(this->at(position)->n_fmt == 0)
			return(g);

		// Iterate through the allowed primitive types for genotypes to collect summary
		// statistics for genotypes at this loci. Information collected includes the
		// base ploidy, if there's any mixed phasing, the number of missing genotypes, and
		// the number of samples that has a special end-of-vector encoding.
		// Only the signed primitives int8_t, int16_t, and int32_t are valid for genotypes.
		switch(this->at(position)->d.fmt[0].type){
		case(BCF_BT_INT8):  g.evaluate<int8_t> (n_samples, this->at(position)->d.fmt[0]); break;
		case(BCF_BT_INT16): g.evaluate<int16_t>(n_samples, this->at(position)->d.fmt[0]); break;
		case(BCF_BT_INT32): g.evaluate<int32_t>(n_samples, this->at(position)->d.fmt[0]); break;
		case(BCF_BT_NULL):
		case(BCF_BT_FLOAT):
		case(BCF_BT_CHAR):
		default:
			std::cerr << "Illegal genotype primtive type: " << io::BCF_TYPE_LOOKUP[this->at(position)->d.fmt[0].type] << std::endl;
		}

		return(g);
	}

	void clear(void){
		uint32_t start_pos = 0;
		if(this->n_carry_over_){
			assert(this->size() != 0);
			std::swap(this->entries_[this->size() - 1], this->entries_[0]);
			start_pos = 1;
			this->n_carry_over_ = 0;
		}

		for(uint32_t i = start_pos; i < this->size(); ++i)
			bcf_clear(this->entries_[i]);

		this->n_entries_ = start_pos;
	}

public:
	uint32_t  n_carry_over_;
	size_type n_entries_;
	size_type n_capacity_;
	pointer*  entries_;
};

}
}



#endif /* IO_HTSLIB_INTEGRATION_H_ */
