#ifndef IO_HTSLIB_INTEGRATION_H_
#define IO_HTSLIB_INTEGRATION_H_

#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

#include "htslib/kstring.h"
#include "htslib/vcf.h"
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

	std::string ToVcfString(const bool is_bcf = false) const{
		// Template:
		// ##contig=<ID=GL000241.1,assembly=b37,length=42152>
		std::string ret = "##contig=<ID=" + this->name;
		if(extra.size()){
			ret += "," + this->extra[0].first + "=" + this->extra[0].second;
			for(uint32_t i = 1; i < this->extra.size(); ++i){
				ret += "," + this->extra[i].first + "=" + this->extra[i].second;
			}
		}
		if(this->description.size()) ret += ",Description=" + this->description;
		ret += ",length=" + std::to_string(this->n_bases);
		if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
		ret += ">";
		return(ret);
	}

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

	std::string ToVcfString(const bool is_bcf = false) const{
		// Template:
		// ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
		std::string ret = "##INFO=<ID=" + this->id;
		ret += ",Number=" + this->number;
		ret += ",Type=" + this->type;
		ret += ",Description=" + this->description;
		if(this->source.size()) ret += ",Source=" + this->source;
		if(this->source.size()) ret += ",Version=" + this->version;
		if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
		ret += ">";
		return(ret);
	}

	std::string ToVcfString(const uint32_t idx) const{
		// Template:
		// ##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
		std::string ret = "##INFO=<ID=" + this->id;
		ret += ",Number=" + this->number;
		ret += ",Type=" + this->type;
		ret += ",Description=" + this->description;
		if(this->source.size()) ret += ",Source=" + this->source;
		if(this->source.size()) ret += ",Version=" + this->version;
		ret += ",IDX=" + std::to_string(idx);
		ret += ">";
		return(ret);
	}

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

	std::string ToVcfString(const bool is_bcf = false) const{
		// Template:
		// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		std::string ret = "##FORMAT=<ID=" + this->id;
		ret += ",Number=" + this->number;
		ret += ",Type=" + this->type;
		ret += ",Description=" + this->description;
		if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
		ret += ">";
		return(ret);
	}

	std::string ToVcfString(const uint32_t idx) const{
		// Template:
		// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		std::string ret = "##FORMAT=<ID=" + this->id;
		ret += ",Number=" + this->number;
		ret += ",Type=" + this->type;
		ret += ",Description=" + this->description;
		ret += ",IDX=" + std::to_string(idx);
		ret += ">";
		return(ret);
	}

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

	std::string ToVcfString(const bool is_bcf = false) const{
		// Template:
		// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		std::string ret = "##FILTER=<ID=" + this->id;
		ret += ",Description=" + this->description;
		if(is_bcf) ret += ",IDX=" + std::to_string(this->idx);
		ret += ">";
		return(ret);
	}

	std::string ToVcfString(const uint32_t idx) const{
		// Template:
		// ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
		std::string ret = "##FILTER=<ID=" + this->id;
		ret += ",Description=" + this->description;
		ret += ",IDX=" + std::to_string(idx);
		ret += ">";
		return(ret);
	}

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

	std::string ToVcfString(void) const{
		// Template:
		// ##source=CombineGVCFs
		std::string ret = "##" + this->key + "=" + this->value;
		return(ret);
	}

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

	std::string ToVcfString(void) const{
		// Template:
		// ##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
		std::string ret = "##" + this->key + "=<";
		ret += this->fields[0].key + "=" + this->fields[0].value;
		for(uint32_t i = 1; i < this->fields.size(); ++i)
			ret += "," + this->fields[i].key + "=" + this->fields[i].value;
		ret += ">";
		return(ret);
	}

	friend std::ostream& operator<<(std::ostream& stream, const VcfStructuredExtra& extra){
		utility::SerializeString(extra.key, stream);
		size_t l_extra = extra.fields.size();
		stream.write((const char*)&l_extra, sizeof(size_t));
		for(uint32_t i = 0; i < extra.fields.size(); ++i)
			stream << extra.fields[i];

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, VcfStructuredExtra& extra){
		utility::DeserializeString(extra.key, stream);
		size_t l_extra;
		stream.read((char*)&l_extra, sizeof(size_t));
		extra.fields.resize(l_extra);
		for(uint32_t i = 0; i < extra.fields.size(); ++i)
			stream >> extra.fields[i];

		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VcfStructuredExtra& extra){
		io::SerializeString(extra.key, buffer);
		size_t l_extra = extra.fields.size();
		io::SerializePrimitive(l_extra, buffer);
		for(uint32_t i = 0; i < extra.fields.size(); ++i)
			buffer << extra.fields[i];

		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VcfStructuredExtra& extra){
		io::DeserializeString(extra.key, buffer);
		size_t l_extra;
		io::DeserializePrimitive(l_extra, buffer);
		extra.fields.resize(l_extra);
		for(uint32_t i = 0; i < extra.fields.size(); ++i)
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
	VcfHeader(const VcfHeader& other);
	~VcfHeader() = default;

	inline size_t GetNumberSamples(void) const{ return(this->samples_.size()); }
	inline size_t GetNumberContigs(void) const{ return(this->contigs_.size()); }

	// Adds Contig information from the idPair to the ContigInfo object.
	void AddContigInfo(const bcf_idpair_t& idPair);

	// Adds FILTER information from the bcf_hrec_t to the VcfFilterInfo object.
	void AddFilterInfo(const bcf_hrec_t* hrec);

	// Adds INFO information from the bcf_hrec_t to the VcfInfo object.
	void AddInfo(const bcf_hrec_t* hrec);

	// Adds FORMAT information from the bcf_hrec_t to the VcfFormatInfo object.
	void AddFormatInfo(const bcf_hrec_t* hrec);

	// Adds structured information from the bcf_hrec_t to the VcfStructuredExtra.
	void AddStructuredExtra(const bcf_hrec_t* hrec);

	// Adds unstructured information from the bcf_hrec_t to the VcfExtra object.
	void AddExtra(const bcf_hrec_t* hrec);

	void AddSample(const std::string& sample_name);

	VcfContig* GetContig(const std::string& name);
	VcfContig* GetContig(const int& idx);
	VcfInfo* GetInfo(const std::string& name);
	VcfInfo* GetInfo(const int& idx);
	VcfFormat* GetFormat(const std::string& name);
	VcfFormat* GetFormat(const int& idx);
	VcfFilter* GetFilter(const std::string& name);
	VcfFilter* GetFilter(const int& idx);
	std::string* GetSample(const std::string& name);

	bool BuildReverseMaps(void);
	bool BuildMaps(void);

	/**<
	 * Recodes the internal IDX field for contig info, INFO, FORMAT, and FILTER
	 * from any range to the range [0, 1, ..., n-1] as desired in Tachyon.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool RecodeIndices(void);

	/**<
	* Converts this header object into a hts_vcf_header object from the
	* internally stored literal string. This object is required for
	* writing out VCF/BCF files.
	* @return
	*/
	hts_vcf_header* ConvertVcfHeader(void);

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

class VcfReader {
public:
	typedef VcfReader self_type;

public:
	// Singleton design pattern for retrieving a guaranteed unique pointer
	// to a VcfReader. This choice is to prevent inadvertent writes to the
	// target file as another file-handle is accessing it.
	static std::unique_ptr<self_type> FromFile(const std::string& variants_path, uint32_t n_extra_threads = 0){
		htsFile* fp = hts_open(variants_path.c_str(), "r");
		if(n_extra_threads){
			int ret = hts_set_threads(fp, n_extra_threads);
			if(ret < 0){
				std::cerr << "failed to open multiple handles" << std::endl;
				return nullptr;
			}
		}

		if (fp == nullptr) {
			std::cerr << utility::timestamp("ERROR")  << "Could not open " << variants_path << std::endl;
			return nullptr;
		}

		bcf_hdr_t* header = bcf_hdr_read(fp);
		if (header == nullptr){
			std::cerr << utility::timestamp("ERROR") << "Couldn't parse header for " << fp->fn << std::endl;
			return nullptr;
		}

		return std::unique_ptr<self_type>(new self_type(variants_path, fp, header));
	}

	bool next(const int unpack_level = BCF_UN_ALL){
		if (bcf_read(this->fp_, this->header_, this->bcf1_) < 0) {
			if (bcf1_->errcode) {
				std::cerr << utility::timestamp("ERROR") << "Failed to parse VCF record: " << bcf1_->errcode << std::endl;
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
				std::cerr << utility::timestamp("ERROR") << "Failed to parse VCF record: " << bcf1_->errcode << std::endl;
				return false;
			} else {
				//std::cerr << utility::timestamp("ERROR") << "Failed to retrieve a htslib bcf1_t record!" << std::endl;
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
        std::cerr << utility::timestamp("ERROR") << "Empty header, not a valid VCF." << std::endl;
        return;
    }

    // Store the file-format header string
    if (std::string(this->header_->hrec[0]->key) != "fileformat") {
        std::cerr << utility::timestamp("ERROR") << "Not a valid VCF, fileformat needed: " << variants_path << std::endl;
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
			std::cerr << utility::timestamp("ERROR") << "Unknown hrec0->type: " << hrec0->type << std::endl;
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

}



#endif /* IO_HTSLIB_INTEGRATION_H_ */
