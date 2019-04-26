/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TACHYON_VARIANT_HEADER_H_
#define TACHYON_VARIANT_HEADER_H_

#include <unordered_map>

#include "tachyon.h"
#include "utility.h"
#include "buffer.h"
#include "support_vcf.h"

#define YON_FOOTER_LENGTH ((TACHYON_FILE_EOF_LENGTH) + sizeof(uint64_t)*3 + sizeof(uint16_t))

namespace tachyon {

class yon_vnt_hdr_t {
public:
	typedef yon_vnt_hdr_t self_type;
	typedef bcf_hdr_t     hts_vcf_header;
	typedef std::unordered_map<std::string, uint32_t> map_type;
	typedef std::unordered_map<uint32_t, uint32_t>    map_reverse_type;

public:
	yon_vnt_hdr_t() = default;
	yon_vnt_hdr_t(const yon_vnt_hdr_t& other);
	~yon_vnt_hdr_t() = default;

	inline size_t GetNumberSamples(void) const { return(this->samples_.size()); }
	inline size_t GetNumberContigs(void) const { return(this->contigs_.size()); }

	void AddSample(const std::string& sample_name);

	const YonContig* GetContig(const std::string& name) const;
	const YonContig* GetContig(const int& idx) const;
	const YonInfo* GetInfo(const std::string& name) const;
	const YonInfo* GetInfo(const int& idx) const;
	const YonFormat* GetFormat(const std::string& name) const;
	const YonFormat* GetFormat(const int& idx) const;
	const VcfFilter* GetFilter(const std::string& name) const;
	const VcfFilter* GetFilter(const int& idx) const;
	const std::string* GetSample(const std::string& name) const;
	int32_t GetSampleId(const std::string& name) const;

	bool BuildReverseMaps(void);
	bool BuildMaps(void);

	/**<
	 * Recodes the internal IDX field for contig info, INFO, FORMAT, and FILTER
	 * from any range to the range [0, 1, ..., n-1] as desired in Tachyon.
	 * @return Returns TRUE upon success or FALSE otherwise.
	 */
	bool RecodeIndices(void);

	void AddGenotypeAnnotationFields(void);
	void AddGenotypeAnnotationFields(const std::vector<std::string>& group_names);

	// Append a string to the literal string
	inline void AppendLiteralString(const std::string& literal_addition) { this->literals_ += literal_addition; }

	// Print the literals and the column header.
	std::ostream& PrintVcfHeader(std::ostream& stream) const;

	std::string ToString(const bool is_bcf = false) const;

	friend std::ostream& operator<<(std::ostream& stream, const yon_vnt_hdr_t& header);
	friend std::istream& operator>>(std::istream& stream, yon_vnt_hdr_t& header);
	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const yon_vnt_hdr_t& header);
	friend yon_buffer_t& operator>>(yon_buffer_t& buffer, yon_vnt_hdr_t& header);

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
	std::vector<VcfFilter>              filter_fields_;
	std::vector<VcfStructuredExtra>     structured_extra_fields_;
	std::vector<VcfExtra>               extra_fields_;

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

struct yon_ftr_t {
public:
	typedef yon_ftr_t self_type;

public:
	yon_ftr_t();
	yon_ftr_t(const char* const data);
	yon_ftr_t(const self_type& other);
	~yon_ftr_t() = default;

	inline const uint64_t& GetEODOffset(void) const { return(this->offset_end_of_data); }
	inline uint64_t& GetEODOffset(void) { return(this->offset_end_of_data); }
	inline const uint64_t& GetNumberBlocks(void) const { return(this->n_blocks); }
	inline uint64_t& GetNumberBlocks(void) { return(this->n_blocks); }
	inline const uint64_t& GetNumberVariants(void) const { return(this->n_variants); }
	inline uint64_t& GetNumberVariants(void) { return(this->n_variants); }
	inline const uint16_t& GetController(void) const { return(this->controller); }
	inline uint16_t& GetController(void) { return(this->controller); }

	bool Validate(void) const;

	friend std::ostream& operator<<(std::ostream& stream, const self_type& yon_ftr_t);
	friend std::istream& operator>>(std::istream& stream, self_type& yon_ftr_t);

public:
	uint64_t  offset_end_of_data;
	uint64_t  n_blocks;
	uint64_t  n_variants;
	uint16_t  controller;
    uint8_t   EOF_marker[TACHYON_FILE_EOF_LENGTH];
};

}

#endif

