#ifndef CORE_BASE_HEADER_YON_TACHYONHEADER_H_
#define CORE_BASE_HEADER_YON_TACHYONHEADER_H_

#include "support/enums.h"
#include "support/helpers.h"
#include "algorithm/OpenHashTable.h"
#include "io/basic_buffer.h"

#include "io/vcf_utils.h"
#include <unordered_map>

namespace tachyon{

struct YonContig : public io::VcfContig {
public:
	YonContig() : n_blocks(0){}
	YonContig(const io::VcfContig& vcf_contig) : io::VcfContig(vcf_contig), n_blocks(0){}
	~YonContig() = default;

	std::string ToVcfString(const bool is_bcf = false) const{ return(io::VcfContig::ToVcfString(is_bcf)); }
	std::string ToVcfString(const uint32_t idx) const{ return(io::VcfContig::ToVcfString(idx)); }

	inline void operator++(void){ ++this->n_blocks; }
	inline void operator--(void){ --this->n_blocks; }
	template <class T> inline void operator+=(const T value){ this->n_blocks += value; }
	template <class T> inline void operator-=(const T value){ this->n_blocks -= value; }

	friend std::ostream& operator<<(std::ostream& stream, const YonContig& contig){
		utility::SerializePrimitive(contig.idx, stream);
		utility::SerializePrimitive(contig.n_bases, stream);
		utility::SerializePrimitive(contig.n_blocks, stream);
		utility::SerializeString(contig.name, stream);
		utility::SerializeString(contig.description, stream);

		size_t size_helper = contig.extra.size();
		utility::SerializePrimitive(size_helper, stream);
		for(uint32_t i = 0; i < contig.extra.size(); ++i){
			utility::SerializeString(contig.extra[i].first, stream);
			utility::SerializeString(contig.extra[i].second, stream);
		}
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, YonContig& contig){
		utility::DeserializePrimitive(contig.idx, stream);
		utility::DeserializePrimitive(contig.n_bases, stream);
		utility::DeserializePrimitive(contig.n_blocks, stream);
		utility::DeserializeString(contig.name, stream);
		utility::DeserializeString(contig.description, stream);

		size_t l_extra;
		utility::DeserializePrimitive(l_extra, stream);
		contig.extra.resize(l_extra);
		for(uint32_t i = 0; i < contig.extra.size(); ++i){
			utility::DeserializeString(contig.extra[i].first, stream);
			utility::DeserializeString(contig.extra[i].second, stream);
		}
		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const YonContig& contig){
		io::SerializePrimitive(contig.idx, buffer);
		io::SerializePrimitive(contig.n_bases, buffer);
		io::SerializePrimitive(contig.n_blocks, buffer);
		io::SerializeString(contig.name, buffer);
		io::SerializeString(contig.description, buffer);

		size_t size_helper = contig.extra.size();
		io::SerializePrimitive(size_helper, buffer);
		for(uint32_t i = 0; i < contig.extra.size(); ++i){
			io::SerializeString(contig.extra[i].first, buffer);
			io::SerializeString(contig.extra[i].second, buffer);
		}
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, YonContig& contig){
		io::DeserializePrimitive(contig.idx, buffer);
		io::DeserializePrimitive(contig.n_bases, buffer);
		io::DeserializePrimitive(contig.n_blocks, buffer);
		io::DeserializeString(contig.name, buffer);
		io::DeserializeString(contig.description, buffer);

		size_t l_extra;
		io::DeserializePrimitive(l_extra, buffer);
		contig.extra.resize(l_extra);
		for(uint32_t i = 0; i < contig.extra.size(); ++i){
			io::DeserializeString(contig.extra[i].first, buffer);
			io::DeserializeString(contig.extra[i].second, buffer);
		}
		return(buffer);
	}

public:
	// Number of Tachyon blocks associated with this contig
	uint32_t n_blocks;
};

class YonInfo : public io::VcfInfo {
public:
	YonInfo() : yon_type(YON_VCF_HEADER_FLAG){}
	YonInfo(const io::VcfInfo& vcf_info) : io::VcfInfo(vcf_info), yon_type(YON_VCF_HEADER_FLAG){
		this->EvaluateType();
	}
	~YonInfo() = default;

	std::string ToVcfString(const bool is_bcf = false) const{ return(io::VcfInfo::ToVcfString(is_bcf)); }
	std::string ToVcfString(const uint32_t idx) const{ return(io::VcfInfo::ToVcfString(idx)); }

	bool EvaluateType(void){
		if(this->type == "Integer") this->yon_type = YON_VCF_HEADER_INTEGER;
		else if(this->type == "Float") this->yon_type = YON_VCF_HEADER_FLOAT;
		else if(this->type == "Flag") this->yon_type = YON_VCF_HEADER_FLAG;
		else if(this->type == "Character") this->yon_type = YON_VCF_HEADER_CHARACTER;
		else if(this->type == "String") this->yon_type = YON_VCF_HEADER_STRING;
		else {
			std::cerr << "Illegal header type: " << this->type << std::endl;
			return false;
 		}
 		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const YonInfo& info){
		utility::SerializePrimitive(info.idx, stream);
		utility::SerializeString(info.id, stream);
		utility::SerializeString(info.number, stream);
		utility::SerializeString(info.type, stream);
		utility::SerializeString(info.description, stream);
		utility::SerializeString(info.source, stream);
		utility::SerializeString(info.version, stream);

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, YonInfo& info){
		utility::DeserializePrimitive(info.idx, stream);
		utility::DeserializeString(info.id, stream);
		utility::DeserializeString(info.number, stream);
		utility::DeserializeString(info.type, stream);
		utility::DeserializeString(info.description, stream);
		utility::DeserializeString(info.source, stream);
		utility::DeserializeString(info.version, stream);
		info.EvaluateType();

		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const YonInfo& info){
		io::SerializePrimitive(info.idx, buffer);
		io::SerializeString(info.id, buffer);
		io::SerializeString(info.number, buffer);
		io::SerializeString(info.type, buffer);
		io::SerializeString(info.description, buffer);
		io::SerializeString(info.source, buffer);
		io::SerializeString(info.version, buffer);

		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, YonInfo& info){
		io::DeserializePrimitive(info.idx, buffer);
		io::DeserializeString(info.id, buffer);
		io::DeserializeString(info.number, buffer);
		io::DeserializeString(info.type, buffer);
		io::DeserializeString(info.description, buffer);
		io::DeserializeString(info.source, buffer);
		io::DeserializeString(info.version, buffer);
		info.EvaluateType();

		return(buffer);
	}

public:
	TACHYON_VARIANT_HEADER_FIELD_TYPE yon_type;
};

class YonFormat : public io::VcfFormat {
public:
	YonFormat() : yon_type(YON_VCF_HEADER_FLAG){}
	YonFormat(const io::VcfFormat& vcf_format) : io::VcfFormat(vcf_format), yon_type(YON_VCF_HEADER_FLAG){
		this->EvaluateType();
	}
	~YonFormat() = default;

	std::string ToVcfString(const bool is_bcf = false) const{ return(io::VcfFormat::ToVcfString(is_bcf)); }
	std::string ToVcfString(const uint32_t idx) const{ return(io::VcfFormat::ToVcfString(idx)); }

	bool EvaluateType(void){
		if(this->type == "Integer") this->yon_type = YON_VCF_HEADER_INTEGER;
		else if(this->type == "Float") this->yon_type = YON_VCF_HEADER_FLOAT;
		else if(this->type == "Character") this->yon_type = YON_VCF_HEADER_CHARACTER;
		else if(this->type == "String") this->yon_type = YON_VCF_HEADER_STRING;
		else {
			std::cerr << "Illegal header type: " << this->type << std::endl;
			return false;
		}
		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const YonFormat& fmt){
		utility::SerializePrimitive(fmt.idx, stream);
		utility::SerializeString(fmt.id, stream);
		utility::SerializeString(fmt.number, stream);
		utility::SerializeString(fmt.type, stream);
		utility::SerializeString(fmt.description, stream);

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, YonFormat& fmt){
		utility::DeserializePrimitive(fmt.idx, stream);
		utility::DeserializeString(fmt.id, stream);
		utility::DeserializeString(fmt.number, stream);
		utility::DeserializeString(fmt.type, stream);
		utility::DeserializeString(fmt.description, stream);
		fmt.EvaluateType();

		return(stream);
	}

	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const YonFormat& fmt){
		io::SerializePrimitive(fmt.idx, buffer);
		io::SerializeString(fmt.id, buffer);
		io::SerializeString(fmt.number, buffer);
		io::SerializeString(fmt.type, buffer);
		io::SerializeString(fmt.description, buffer);

		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, YonFormat& fmt){
		io::DeserializePrimitive(fmt.idx, buffer);
		io::DeserializeString(fmt.id, buffer);
		io::DeserializeString(fmt.number, buffer);
		io::DeserializeString(fmt.type, buffer);
		io::DeserializeString(fmt.description, buffer);
		fmt.EvaluateType();

		return(buffer);
	}

public:
	TACHYON_VARIANT_HEADER_FIELD_TYPE yon_type;
};

class VariantHeader{
public:
	typedef VariantHeader self_type;
	typedef bcf_hdr_t     hts_vcf_header;
	typedef std::unordered_map<std::string, uint32_t> map_type;
	typedef std::unordered_map<uint32_t, uint32_t>    map_reverse_type;

public:
	VariantHeader() = default;
	VariantHeader(const VariantHeader& other) :
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

	VariantHeader(const io::VcfHeader& other) :
		fileformat_string_(other.fileformat_string_),
		literals_(other.literals_),
		samples_(other.samples_),
		filter_fields_(other.filter_fields_),
		structured_extra_fields_(other.structured_extra_fields_),
		extra_fields_(other.extra_fields_)
	{
		this->BuildMaps();
		this->BuildReverseMaps();

		this->contigs_.resize(other.contigs_.size());
		for(uint32_t i = 0; i < other.contigs_.size(); ++i)
			this->contigs_[i] = other.contigs_[i];

		this->info_fields_.resize(other.info_fields_.size());
		for(uint32_t i = 0; i < other.info_fields_.size(); ++i)
			this->info_fields_[i] = other.info_fields_[i];

		this->format_fields_.resize(other.format_fields_.size());
		for(uint32_t i = 0; i < other.format_fields_.size(); ++i)
			this->format_fields_[i] = other.format_fields_[i];

		this->RecodeIndices();
	}

	~VariantHeader() = default;

	inline size_t GetNumberSamples(void) const{ return(this->samples_.size()); }
	inline size_t GetNumberContigs(void) const{ return(this->contigs_.size()); }

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

	const YonContig* GetContig(const std::string& name) const {
		map_type::const_iterator it = this->contigs_map_.find(name);
		if(it != this->contigs_map_.end()) return(&this->contigs_[it->second]);
		return(nullptr);
	}

	const YonContig* GetContig(const int& idx) const {
		map_reverse_type::const_iterator it = this->contigs_reverse_map_.find(idx);
		if(it != this->contigs_reverse_map_.end()) return(&this->contigs_[it->second]);
		return(nullptr);
	}

	const YonInfo* GetInfo(const std::string& name) const {
		map_type::const_iterator it = this->info_fields_map_.find(name);
		if(it != this->info_fields_map_.end()) return(&this->info_fields_[it->second]);
		return(nullptr);
	}

	const YonInfo* GetInfo(const int& idx) const {
		map_reverse_type::const_iterator it = this->info_fields_reverse_map_.find(idx);
		if(it != this->info_fields_reverse_map_.end()) return(&this->info_fields_[it->second]);
		return(nullptr);
	}

	const YonFormat* GetFormat(const std::string& name) const {
		map_type::const_iterator it = this->format_fields_map_.find(name);
		if(it != this->format_fields_map_.end()) return(&this->format_fields_[it->second]);
		return(nullptr);
	}

	const YonFormat* GetFormat(const int& idx) const {
		map_reverse_type::const_iterator it = this->format_fields_reverse_map_.find(idx);
		if(it != this->format_fields_reverse_map_.end()) return(&this->format_fields_[it->second]);
		return(nullptr);
	}

	const io::VcfFilter* GetFilter(const std::string& name) const {
		map_type::const_iterator it = this->filter_fields_map_.find(name);
		if(it != this->filter_fields_map_.end()) return(&this->filter_fields_[it->second]);
		return(nullptr);
	}

	const io::VcfFilter* GetFilter(const int& idx) const {
		map_reverse_type::const_iterator it = this->filter_fields_reverse_map_.find(idx);
		if(it != this->filter_fields_reverse_map_.end()) return(&this->filter_fields_[it->second]);
		return(nullptr);
	}

	const std::string* GetSample(const std::string& name) const {
		map_type::const_iterator it = this->samples_map_.find(name);
		if(it != this->samples_map_.end()) return(&this->samples_[it->second]);
		return(nullptr);
	}

	int32_t GetSampleId(const std::string& name) const {
		map_type::const_iterator it = this->samples_map_.find(name);
		if(it != this->samples_map_.end()) return(it->second);
		return(-1);
	}

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
	hts_vcf_header* ConvertVcfHeaderLiterals(const bool add_format = true);
	hts_vcf_header* ConvertVcfHeader(const bool add_format = true);

	void AddGenotypeAnnotationFields(void);
	void AddGenotypeAnnotationFields(const std::vector<std::string>& group_names);

	// Append a string to the literal string
	inline void AppendLiteralString(const std::string& literal_addition){ this->literals_ += literal_addition; }

	// Print the literals and the column header.
	std::ostream& PrintVcfHeader(std::ostream& stream) const;

	std::string ToString(const bool is_bcf = false) const;

	friend std::ostream& operator<<(std::ostream& stream, const VariantHeader& header);
	friend std::istream& operator>>(std::istream& stream, VariantHeader& header);
	friend io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VariantHeader& header);
	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VariantHeader& header);

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
	std::vector<std::string>            samples_;
	std::vector<YonContig>              contigs_;
	std::vector<YonInfo>                info_fields_;
	std::vector<YonFormat>              format_fields_;
	std::vector<io::VcfFilter>          filter_fields_;
	std::vector<io::VcfStructuredExtra> structured_extra_fields_;
	std::vector<io::VcfExtra>           extra_fields_;

	// Utility members
	//
	// Hash tables allowing the mapping from the unique identifier string
	// (such as contig name) to the relative index offset of that object.
	// This approach requires another layer of indirection when mapping
	// from the index to the actual target. For example:
	//
	// contigs_[contigs_map_["chr20"].second] <- maps to the actual target
	//
	// The reverse maps allows the mapping from a unique IDX identifier
	// to the relative index offset of that object. As above, this requires
	// an addition indirect lookup to access the desired object. For example
	// mapping the first occuring FORMAT field to its name:
	//
	// format_fields_[format_fields_reverse_map_[container.at(0)->d.fmt[0].id].second].id
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

}



#endif /* CORE_BASE_HEADER_YON_TACHYONHEADER_H_ */
