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

#include "support/magic_constants.h"
#include "support/helpers.h"
#include "buffer.h"
#include "support_vcf.h"

namespace tachyon {

class VariantHeader {
public:
	typedef VariantHeader self_type;
	typedef bcf_hdr_t     hts_vcf_header;
	typedef std::unordered_map<std::string, uint32_t> map_type;
	typedef std::unordered_map<uint32_t, uint32_t>    map_reverse_type;

public:
	VariantHeader() = default;
	VariantHeader(const VariantHeader& other);
	//VariantHeader(const io::VcfHeader& other);
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

	const VcfFilter* GetFilter(const std::string& name) const {
		map_type::const_iterator it = this->filter_fields_map_.find(name);
		if(it != this->filter_fields_map_.end()) return(&this->filter_fields_[it->second]);
		return(nullptr);
	}

	const VcfFilter* GetFilter(const int& idx) const {
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
	friend yon_buffer_t& operator<<(yon_buffer_t& buffer, const VariantHeader& header);
	friend yon_buffer_t& operator>>(yon_buffer_t& buffer, VariantHeader& header);

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

}

#endif

