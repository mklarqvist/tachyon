#ifndef CORE_TACHYON_READER_H_
#define CORE_TACHYON_READER_H_

#include <cmath>

#include "zstd.h"
#include "zstd_errors.h"

#include "algorithm/compression/compression_manager.h"
#include "algorithm/encryption/EncryptionDecorator.h"
#include "algorithm/timer.h"
#include "containers/format_container.h"
#include "containers/format_container_string.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/info_container_string.h"
#include "containers/primitive_group_container.h"
#include "containers/meta_container.h"
#include "algorithm/digital_digest.h"
#include "containers/variantblock.h"
#include "core/genotype_object.h"
#include "core/footer/footer.h"
#include "core/header/variant_header.h"
#include "math/fisher_math.h"
#include "math/square_matrix.h"
#include "math/basic_vector_math.h"
#include "utility/support_vcf.h"
#include "index/index.h"

namespace tachyon{

struct VariantReaderSettings{
public:
	typedef VariantReaderSettings self_type;
	typedef index::IndexEntry     index_entry_type;

public:
	VariantReaderSettings() :
		drop_format(false),
		header_only(false),
		show_header(true),
		custom_delimiter(false),
		custom_delimiter_char(0),
		custom_output_format(false),
		filter_any(false),
		filter_all(false),
		annotate_genotypes(false),
		output_FORMAT_as_vector(false)
	{}
	~VariantReaderSettings() = default;

	/**<
	 * Validates the regex patterns of the given interval strings.
	 * Note that this function only checks if the strings have correct
	 * syntax and does NOT parse or check the existence of these regions
	 * @return Returns TRUE if all strings are valid interval strings or FALSE otherwise
	 */
	bool validateIntervalStrings(void){
		if(this->interval_strings.size() == 0)
			return true;

		for(U32 i = 0; i < this->interval_strings.size(); ++i){
			// scrub whitespace
			//this->interval_strings[i].erase(remove_if(this->interval_strings[i].begin(), this->interval_strings[i].end(), isspace), this->interval_strings[i].end());
			this->interval_strings[i] = utility::remove_whitespace(this->interval_strings[i]);

			if (std::regex_match (this->interval_strings[i], constants::YON_REGEX_CONTIG_ONLY )){
				//std::cerr << "chromosome onlu" << std::endl;
			} else if (std::regex_match (this->interval_strings[i], constants::YON_REGEX_CONTIG_POSITION )){
				//std::cerr << "chromosome pos" << std::endl;
			} else if (std::regex_match (this->interval_strings[i], constants::YON_REGEX_CONTIG_RANGE )){
				//std::cerr << "chromosome pos - pos" << std::endl;
			} else {
				std::cerr << utility::timestamp("ERROR") << "Uninterpretable interval string: " << this->interval_strings[i] << std::endl;
				return false;
			}
		}

		return true;
	}

	/**<
	 * Construct a string with the internal interpreted parameters
	 * @return Returns a string
	 */
	inline std::string get_settings_string(void) const{
		return(std::string(
		"##tachyon_viewCommandSettings={"
		"\"input\":" + (this->input.length() ? "\"" + this->input + "\"" : "null") +
		",\"output\":" + (this->output.length() ? "\"" + this->output + "\"" : "null") +
		",\"keychain_file\":" + (this->keychain_file.length() ? "\"" + this->keychain_file + "\"" : "null") +
		",\"annotate_genotypes\":" + (this->annotate_genotypes ? "true" : "false") +
		",\"drop_format\":" + (this->drop_format ? "true" : "false") +
		"}; timestamp=" + tachyon::utility::datetime()
 		));
	}

	inline const bool hasBlockList(void) const{ return(this->yon_blocks_list.size()); }
	inline std::vector<index_entry_type>& getBlockList(void){ return(this->yon_blocks_list); }
	inline const std::vector<index_entry_type>& getBlockList(void) const{ return(this->yon_blocks_list); }

public:
	bool drop_format; // drop FORMAT fields
	bool header_only; // show only the VCF header
	bool show_header; // show the VCF header
	bool custom_delimiter; // output uses a custom delimiter
	char custom_delimiter_char; // what is the custom delimiter
	bool custom_output_format;  // output has a custom format
	bool filter_any;  // filter output
	bool filter_all;  // filter
	bool annotate_genotypes;
	bool output_FORMAT_as_vector;
	std::string input;
	std::string output;
	std::string keychain_file;
	std::string output_type;
	std::vector<std::string> interval_strings;
	std::string sample_names_file;
	std::vector<std::string> sample_names;
	std::vector<index_entry_type> yon_blocks_list;
};

struct VariantReaderObjects{
public:
	typedef VariantReaderObjects                 self_type;
	typedef containers::MetaContainer            meta_container_type;
	typedef containers::GenotypeContainer        gt_container_type;
	typedef containers::InfoContainerInterface   info_interface_type;
	typedef containers::FormatContainerInterface format_interface_type;
	typedef containers::GenotypeSummary          genotype_summary_type;

public:
	VariantReaderObjects() :
		loaded_genotypes(false),
		loaded_meta(false),
		n_loaded_info(0),
		n_loaded_format(0),
		meta(nullptr),
		genotypes(nullptr),
		genotype_summary(nullptr),
		info_fields(nullptr),
		format_fields(nullptr)
	{}

	~VariantReaderObjects(){
		delete this->meta;
		delete this->genotypes;
		delete this->genotype_summary;

		for(U32 i = 0; i < this->n_loaded_info; ++i) delete this->info_fields[i];
		delete [] this->info_fields;

		for(U32 i = 0; i < this->n_loaded_format; ++i) delete this->format_fields[i];
		delete [] this->format_fields;
	}

public:
	bool loaded_genotypes;
	bool loaded_meta;
	size_t n_loaded_info;
	size_t n_loaded_format;
	std::vector<U32> info_keep;
	std::vector<U32> format_keep;
	std::vector< U16 > additional_info_execute_flag_set;
	std::vector< std::vector<U32> > local_match_keychain_info;
	std::vector< std::vector<U32> > local_match_keychain_format;

	std::vector<std::string> info_field_names;
	std::vector<std::string> format_field_names;

	meta_container_type*     meta;
	gt_container_type*       genotypes;
	genotype_summary_type*   genotype_summary;
	info_interface_type**    info_fields;
	format_interface_type**  format_fields;
};

template <class ValueClass>
struct VariantReaderFiltersTuple{
public:
	typedef VariantReaderFiltersTuple<ValueClass> self_type;
	typedef bool (self_type::*filter_function)(const ValueClass& target, const ValueClass& limit) const;

public:
	VariantReaderFiltersTuple() :
		filter(false),
		r_value(0),
		comparator(&self_type::__filterGreaterEqual)
	{}

	VariantReaderFiltersTuple(const ValueClass& r_value) :
		filter(true),
		r_value(r_value),
		comparator(&self_type::__filterGreaterEqual)
	{}

	VariantReaderFiltersTuple(const ValueClass& r_value, const TACHYON_COMPARATOR_TYPE& comparator) :
		filter(true),
		r_value(r_value),
		comparator(nullptr)
	{
		switch(comparator){
		case(YON_CMP_GREATER):       this->comparator = &self_type::__filterGreater;      break;
		case(YON_CMP_GREATER_EQUAL): this->comparator = &self_type::__filterGreaterEqual; break;
		case(YON_CMP_LESS):          this->comparator = &self_type::__filterLesser;       break;
		case(YON_CMP_LESS_EQUAL):    this->comparator = &self_type::__filterLesserEqual;  break;
		case(YON_CMP_EQUAL):         this->comparator = &self_type::__filterEqual;        break;
		case(YON_CMP_NOT_EQUAL):     this->comparator = &self_type::__filterNotEqual;     break;
		}
	}

	void operator()(const ValueClass& r_value){
		this->filter = true;
		this->r_value = r_value;
	}

	void operator()(const ValueClass& r_value, const TACHYON_COMPARATOR_TYPE& comparator){
		this->filter  = true;
		this->r_value = r_value;

		switch(comparator){
		case(YON_CMP_GREATER):       this->comparator = &self_type::__filterGreater;      break;
		case(YON_CMP_GREATER_EQUAL): this->comparator = &self_type::__filterGreaterEqual; break;
		case(YON_CMP_LESS):          this->comparator = &self_type::__filterLesser;       break;
		case(YON_CMP_LESS_EQUAL):    this->comparator = &self_type::__filterLesserEqual;  break;
		case(YON_CMP_EQUAL):         this->comparator = &self_type::__filterEqual;        break;
		case(YON_CMP_NOT_EQUAL):     this->comparator = &self_type::__filterNotEqual;     break;
		}
	}

	~VariantReaderFiltersTuple() = default;

	inline bool applyFilter(const ValueClass& l_value) const{ return((this->*comparator)(l_value, r_value)); }

	// Comparator functions
	inline bool __filterLesser(const ValueClass& target, const ValueClass& limit) const{return(target < limit);}
	inline bool __filterLesserEqual(const ValueClass& target, const ValueClass& limit) const{return(target <= limit);}
	inline bool __filterGreater(const ValueClass& target, const ValueClass& limit) const{return(target > limit);}
	inline bool __filterGreaterEqual(const ValueClass& target, const ValueClass& limit) const{return(target >= limit);}
	inline bool __filterEqual(const ValueClass& target, const ValueClass& limit) const{return(target == limit);}
	inline bool __filterNotEqual(const ValueClass& target, const ValueClass& limit) const{return(target != limit);}

public:
	bool            filter;
	ValueClass      r_value;
	filter_function comparator;
};

struct VariantReaderFilters{
public:
	typedef VariantReaderFilters self_type;
	typedef VariantReaderObjects objects_type;
	typedef bool (self_type::*filter_function)(const objects_type& objects, const U32& position) const;
	typedef bool (self_type::*family_filter_function)(void) const;

	typedef bool (self_type::*filter_float)(const float& target, const float& limit) const;
	typedef bool (self_type::*filter_integer)(const S32& target, const S32& limit) const;

	// example:
	// TACHYON_COMPARATOR_TYPE::YON_CMP_EQUAL

public:
	VariantReaderFilters() :
		target_intervals(false)
	{

	}

	~VariantReaderFilters() = default;

	// Has mixed phasing
	inline bool filterMixedPhasing(const objects_type& objects, const U32& position) const{
		assert(objects.meta != nullptr);
		return(objects.meta->at(position).isGTMixedPhasing() == true);
	}

	// GT data matches this
	inline bool filterUniformMatchPhase(const objects_type& objects,
												 const U32& position) const
	{
		assert(objects.meta != nullptr);
		return(objects.meta->at(position).isGTMixedPhasing() == false &&
			   objects.meta->at(position).controller.gt_phase == this->filter_uniform_phase.r_value);
	}

	bool filterUncalled(const objects_type& objects, const U32& position) const;
	bool filterPloidy(const objects_type& objects, const U32& position) const;
	bool filterSampleList(const objects_type& objects, const U32& position) const;

	// Use custom AVL tree for position
	bool filterKnownNovel(const objects_type& objects, const U32& position) const;

	// BCFtools calculate this as the SUM of all ALT counts
	// We filter based on ANY ALT frequency OPERATOR the target frequency
	bool filterAlleleFrequency(const objects_type& objects, const U32& position) const{
		const std::vector<double> af = objects.genotype_summary->calculateAlleleFrequency(objects.meta->at(position));
		for(U32 i = 1; i < af.size(); ++i){
			if(this->filter_af.applyFilter(af[i]))
				return true;
			//if((this->*(this->allele_frequency_comparator))(af[i], this->allele_frequency))
			//	return(true);
		}
		return(false);
	}

	bool filterVariantClassification(const objects_type& objects, const U32& position) const;
	bool filterUnseenAlternativeAlleles(const objects_type& objects, const U32& position) const;
	bool filterRegions(const objects_type& objects) const; // Filter by target intervals
	bool filterFILTER(const objects_type& objects) const;  // Filter by desired FILTER values
	bool filterINFO(const objects_type& objects) const;    // custom filter. e.g. AC<1024

	inline bool filterAlternativeAlleles(const objects_type& object, const U32& position) const{
		// Remove one to total count as REF is counted here
		// Recast as signed integer to avoid possible underflowing issues
		return(this->filter_n_alts.applyFilter(object.meta->at(position).getNumberAlleles() - 1));

		//return((this->*(this->n_alt_alleles_comparator))((S32)object.meta->at(position).getNumberAlleles() - 1, this->n_alt_alleles));
	}

	inline bool filterHasMissingGenotypes(const objects_type& object, const U32& position) const{
		return(this->filter_missing.applyFilter(object.meta->at(position).controller.gt_anyMissing));
	}

	/**<
	 * Constructs the filter pointer vector given the fields that have been set
	 * @return Returns TRUE if passing construction or FALSE otherwise
	 */
	bool build(void){
		this->filters.clear();
		if(this->filter_n_alts.filter)    this->filters.push_back(&self_type::filterAlternativeAlleles);
		//if(this->mixed_phasing_only)      this->filters.push_back(&self_type::filterMixedPhasing);
		if(this->filter_missing.filter)   this->filters.push_back(&self_type::filterHasMissingGenotypes);
		if(this->filter_af.filter) this->filters.push_back(&self_type::filterAlleleFrequency);
		if(this->filter_uniform_phase.filter) this->filters.push_back(&self_type::filterUniformMatchPhase);
		return true;
	}

	/**<
	 * Iteratively apply filters in the filter pointer vector
	 * @param objects  Target objects container structure
	 * @param position Target position (relative loci) in the container
	 * @return         Returns TRUE if passes filtering or FALSE otherwise
	 */
	bool filter(const objects_type& objects, const U32 position) const{
		for(U32 i = 0 ; i < this->filters.size(); ++i){
			// Todo: invoke this only when necessary AND possible
			objects.genotypes->at(position).getSummary(*objects.genotype_summary);
			if((this->*(this->filters[i]))(objects, position) == false){
				return false;
			}
		}
		return true;
	}

public:
	// Cannot template this
	inline bool __filterLesser(const float& target, const float& limit) const{return(target < limit);}
	inline bool __filterLesserEqual(const float& target, const float& limit) const{return(target <= limit);}
	inline bool __filterGreater(const float& target, const float& limit) const{return(target > limit);}
	inline bool __filterGreaterEqual(const float& target, const float& limit) const{return(target >= limit);}
	inline bool __filterEqual(const float& target, const float& limit) const{return(target == limit);}
	inline bool __filterNotEqual(const float& target, const float& limit) const{return(target != limit);}

	inline bool __filterLesser(const S32& target, const S32& limit) const{return(target < limit);}
	inline bool __filterLesserEqual(const S32& target, const S32& limit) const{return(target <= limit);}
	inline bool __filterGreater(const S32& target, const S32& limit) const{return(target > limit);}
	inline bool __filterGreaterEqual(const S32& target, const S32& limit) const{return(target >= limit);}
	inline bool __filterEqual(const S32& target, const S32& limit) const{return(target == limit);}
	inline bool __filterNotEqual(const S32& target, const S32& limit) const{return(target != limit);}

public:
	bool target_intervals;
	// std::vector<intervals> intervals;

	std::vector<filter_function> filters;
	std::vector<family_filter_function> family_filters;

	VariantReaderFiltersTuple<bool>  filter_uniform_phase;
	VariantReaderFiltersTuple<SBYTE> filter_n_alts;
	VariantReaderFiltersTuple<bool>  filter_missing;
	VariantReaderFiltersTuple<float> filter_af;
};

class VariantReader{
	typedef VariantReader                          self_type;
	typedef io::BasicBuffer                        buffer_type;
	typedef core::VariantHeader                    header_type;
	typedef core::Footer                           footer_type;
	typedef algorithm::CompressionManager          codec_manager_type;
	typedef DataBlockSettings                      block_settings_type;
	typedef VariantReaderSettings                  settings_type;
	typedef index::Index                           index_type;
	typedef index::IndexEntry                      index_entry_type;
	typedef algorithm::VariantDigitalDigestManager checksum_type;
	typedef encryption::Keychain                   keychain_type;
	typedef core::MetaEntry                        meta_entry_type;
	typedef VariantReaderObjects                   objects_type;
	typedef containers::VariantBlock               block_entry_type;
	typedef containers::MetaContainer              meta_container_type;
	typedef containers::GenotypeContainer          gt_container_type;
	typedef containers::InfoContainerInterface     info_interface_type;
	typedef containers::FormatContainerInterface   format_interface_type;
	typedef containers::GenotypeSummary            genotype_summary_type;

	// Function pointers
	typedef void (self_type::*print_format_function)(buffer_type& buffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	typedef void (self_type::*print_info_function)(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const;
	typedef void (self_type::*print_filter_function)(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;
	typedef buffer_type& (*print_meta_function)(buffer_type& buffer, const char& delimiter, const meta_entry_type& meta_entry, const header_type& header, const SettingsCustomOutput& controller);

public:
	VariantReader();
	VariantReader(const std::string& filename);
	VariantReader(const self_type& other);
	~VariantReader();

	/**<
	 * Retrieve current settings records. Settings object is used
	 * prior to each read of a `block` and can be freely modified
	 * before each call. Modifying the settings after a read block
	 * has been invoked has no effect on the loaded `block` data.
	 * @return A reference instance of the block settings object
	 */
	inline block_settings_type& getBlockSettings(void){ return(this->block_settings); }

	/**<
	 * Retrieve current settings for the variant reader. This settings
	 * object controls the parsing/output of the reader itself. This is
	 * unlike the `block_settings_type` that controls the `DataBlock`
	 * parsing.
	 * @return A reference instance of the settings object
	 */
	inline settings_type& getSettings(void){ return(this->settings); }

	/**<
	 * Checks if a FORMAT `field` is set in the header and then checks
	 * if that field exists in the current block. If it does return
	 * the local key. If the field is not described in the header at
	 * all then return -2. If it exists in the header but not in the
	 * loaded block then return -1.
	 * @param field_name FORMAT field name to search for (e.g. "GL")
	 * @return Returns local key if found in this block. Returns -2 if not found in header, or -1 if found in header but not in block
	 */
	const int has_format_field(const std::string& field_name) const;

	/**<
	 * Checks if a INFO `field` is set in the header and then checks
	 * if that field exists in the current block. If it does return
	 * the local key. If the field is not described in the header at
	 * all then return -2. If it exists in the header but not in the
	 * loaded block then return -1.
	 * @param field_name INFO field name to search for (e.g. "AC")
	 * @return Returns local key if found in this block. Returns -2 if not found in header, or -1 if found in header but not in block
	 */
	const int has_info_field(const std::string& field_name) const;

	/**<
	 * Checks if a FILTER `field` is set in the header and then checks
	 * if that field exists in the current block. If it does return
	 * the local key. If the field is not described in the header at
	 * all then return -2. If it exists in the header but not in the
	 * loaded block then return -1.
	 * @param field_name FILTER field name to search for (e.g. "PASS")
	 * @return Returns local key if found in this block. Returns -2 if not found in header, or -1 if found in header but not in block
	 */
	const int has_filter_field(const std::string& field_name) const;

	/**<
	 * Calculates which INFO pattern matches are found for the given field
	 * name in the current loaded block.
	 * @param field_name INFO field name
	 * @return           Returns a vector of booleans representing pattern matches
	 */
	const std::vector<bool> get_info_field_pattern_matches(const std::string& field_name) const;

	/**<
	 * Calculates which FORMAT pattern matches are found for the given field
	 * name in the current loaded block.
	 * @param field_name FORMAT field name
	 * @return           Returns a vector of booleans representing pattern matches
	 */
	const std::vector<bool> get_format_field_pattern_matches(const std::string& field_name) const;

	/**<
	 * Factory function for FORMAT container given an input `field` name
	 * @param field_name FORMAT field name to create a container for
	 * @return           Returns an instance of a `FormatContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::FormatContainer<T>* get_format_container(const std::string& field_name) const{
		int format_field = this->has_format_field(field_name);
		if(format_field >= 0) return(new containers::FormatContainer<T>(this->block.format_containers[format_field], this->header.getSampleNumber()));
		else return nullptr;
	}

	/**<
	 * Factory function for balanced FORMAT container given an input `field` name.
	 * Balances the container such that variants that do not have the given field
	 * will have an empty container placed instead.
	 * @param field_name     FORMAT field name to create a container for
	 * @param meta_container Container for meta objects used to balance the FORMAT container
	 * @return               Returns an instance of a balanced `FormatContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::FormatContainer<T>* get_balanced_format_container(const std::string& field_name, const containers::MetaContainer& meta_container) const{
		int format_field = this->has_format_field(field_name);
		if(format_field >= 0){
			const std::vector<bool> pattern_matches = this->get_format_field_pattern_matches(field_name);
			U32 matches = 0;
			for(U32 i = 0; i < pattern_matches.size(); ++i)
				matches += pattern_matches[i];

			if(matches == 0)
				return nullptr;

			return(new containers::FormatContainer<T>(this->block.format_containers[format_field], meta_container, pattern_matches, this->header.getSampleNumber()));
		}
		else return nullptr;
	}

	/**<
	 * Factory function for INFO container given an input `field` name
	 * @param field_name INFO field name to create a container for
	 * @return           Returns an instance of a `InfoContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::InfoContainer<T>* get_info_container(const std::string& field_name) const{
		int info_field = this->has_info_field(field_name);
		if(info_field >= 0) return(new containers::InfoContainer<T>(this->block.info_containers[info_field]));
		else return nullptr;
	}

	/**<
	 * Factory function for balanced INFO container given an input `field` name.
	 * Balances the container such that variants that do not have the given field
	 * will have an empty container placed instead.
	 * @param field_name     INFO field name to create a container for
	 * @param meta_container Container for meta objects used to balance the INFO container
	 * @return               Returns an instance of a balanced `InfoContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::InfoContainer<T>* get_balanced_info_container(const std::string& field_name, const containers::MetaContainer& meta_container) const{
		int info_field = this->has_info_field(field_name);
		if(info_field >= 0){
			const std::vector<bool> pattern_matches = this->get_info_field_pattern_matches(field_name);

			U32 matches = 0;
			for(U32 i = 0; i < pattern_matches.size(); ++i)
				matches += pattern_matches[i];

			if(matches == 0)
				return nullptr;

			return(new containers::InfoContainer<T>(this->block.info_containers[info_field], meta_container, pattern_matches));
		}
		else return nullptr;
	}

	/**<
	 * Opens a YON file. Performs all prerequisite
	 * checks and loads all auxiliary data structures
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool open(void);

	/**<
	 * Opens a YON file. Performs all prerequisite
	 * checks and loads all auxiliary data structures
	 * @param filename Target input filename
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	inline bool open(const std::string& filename){
		this->settings.input = filename;
		return(this->open());
	}

	/**<
	 * Overloaded operator for blocks. Useful when
	 * looping over a range of blocks. This happens
	 * frequently in parallel programs.
	 * @param index Block index value in range [0..n_blocks)
	 * @return      Returns TRUE if operation was successful or FALSE otherwise
	 */
	bool operator[](const U32 position);

	/**<
	 *
	 * @param position
	 * @return
	 */
	bool seektoBlock(const U32 position);

	/**<
	 *
	 * @param chromosome_name
	 * @return
	 */
	bool seekToBlockChromosome(const std::string& chromosome_name);

	/**<
	 *
	 * @param chromosome_name
	 * @param from_bp_position
	 * @param to_bp_position
	 * @return
	 */
	bool seekToBlockChromosome(const std::string& chromosome_name, const U32 from_bp_position, const U32 to_bp_position);

	/**<
	 * Get the next YON block in-order
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool nextBlock(void);

	/**<
	 * Get the target YON block
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool getBlock(const index_entry_type& index_entry);

	/**<
	 * Get the current YON block in-order as a copy
	 * @return Returns a YON block. The container has a size of 0 upon fail/empty
	 */
	block_entry_type getBlock(void);


	/**<
	 * Seeks to a specific YON block without loading anything.
	 * This allows the user to seek to a specific block and
	 * change the block_settings (i.e. what fields to load) and
	 * then invoke nextBlock() for example.
	 * @param blockID
	 * @return
	 */
	bool seek_to_block(const U32& blockID);

	/**<
	 *
	 * @return
	 */
	bool parseSettings(void);

	/**<
	 * Primary construction function for generating the appropriate instances of
	 * iterators / containers
	 * @param objects Target objects
	 * @return        Returns reference to input target objects
	 */
	objects_type& loadObjects(objects_type& objects) const;

	/**<
	 * Outputs
	 * @return
	 */
	const U64 outputVCF(void){
		U64 n_variants = 0;

		if(this->block_settings.annotate_extra){
			// fixme
			// if special
			// "FS_A", "AN", "NM", "NPM", "AC", "AC_FW", "AC_REV", "AF", "HWE_P", "VT", "MULTI_ALLELIC", "F_PIC"
			if(this->header.getInfoField("FS_A") == false)          this->header.literals += "\n##INFO=<ID=FS_A,Number=A,Type=Float>";
			if(this->header.getInfoField("AN") == false)            this->header.literals += "\n##INFO=<ID=AN,Number=A,Type=Integer>";
			if(this->header.getInfoField("NM") == false)            this->header.literals += "\n##INFO=<ID=NM,Number=A,Type=Integer>";
			if(this->header.getInfoField("NPM") == false)           this->header.literals += "\n##INFO=<ID=NPM,Number=A,Type=Integer>";
			if(this->header.getInfoField("AC") == false)            this->header.literals += "\n##INFO=<ID=AC,Number=A,Type=Integer>";
			if(this->header.getInfoField("AC_FWD") == false)        this->header.literals += "\n##INFO=<ID=AC_FWD,Number=A,Type=Integer>";
			if(this->header.getInfoField("AC_REV") == false)        this->header.literals += "\n##INFO=<ID=AC_REV,Number=A,Type=Integer>";
			if(this->header.getInfoField("HWE_P") == false)         this->header.literals += "\n##INFO=<ID=HWE_P,Number=A,Type=Float>";
			if(this->header.getInfoField("VT") == false)            this->header.literals += "\n##INFO=<ID=VT,Number=A,Type=String>";
			if(this->header.getInfoField("AF") == false)            this->header.literals += "\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1)\">";
			if(this->header.getInfoField("MULTI_ALLELIC") == false) this->header.literals += "\n##INFO=<ID=MULTI_ALLELIC,Number=0,Type=Flag>";
			if(this->header.getInfoField("F_PIC") == false)         this->header.literals += "\n##INFO=<ID=F_PIC,Number=A,Type=Float,Description=\"Population inbreeding coefficient (F-statistics)\">";
		}

		this->header.literals += "\n##tachyon_viewVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
		this->header.literals += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
				  + SSLeay_version(SSLEAY_VERSION) + "," + "ZSTD-" + ZSTD_versionString() + "; timestamp=" + utility::datetime();

		this->header.literals += "\n##tachyon_viewCommand=" + tachyon::constants::LITERAL_COMMAND_LINE + '\n';
		this->header.literals += this->getSettings().get_settings_string();

		// Output VCF header
		if(this->block_settings.show_vcf_header){
			this->header.writeVCFHeaderString(std::cout, this->block_settings.load_format || this->block_settings.format_list.size());
		}

		// If seek is active for targetted intervals
		if(this->getSettings().interval_strings.size()){
			if(this->getSettings().getBlockList().size()){
				for(U32 i = 0; i < this->getSettings().getBlockList().size(); ++i){
					if(this->getBlock(this->getSettings().getBlockList()[i]) == false){
						return(0);
					}
					n_variants += this->outputBlockVCF();
				}

				return(n_variants);
			} else { // Provided intervals but no matching YON blocks
				return(0);
			}
		}

		// While there are YON blocks
		while(this->nextBlock()) n_variants += this->outputBlockVCF();
		return(n_variants);
	}

	/**<
	 *
	 * @return
	 */
	const U64 outputCustom(void){
		U64 n_variants = 0;

		// While there are YON blocks
		while(this->nextBlock()) n_variants += this->outputBlockCustom();
		return(n_variants);
	}

	/**<
	 *
	 * @return
	 */
	const U32 outputBlockVCF(void) const{
		objects_type objects;
		this->loadObjects(objects);

		// Todo
		VariantReaderFilters filters;

		//filters.filter_af(0.9, YON_CMP_GREATER_EQUAL);
		//filters.filter_n_alts(5, YON_CMP_GREATER_EQUAL);
		//filters.filter_uniform_phase(true);
		filters.filter_missing(false);

		filters.build();

		// Reserve memory for output buffer
		// This is much faster than writing directly to ostream because of syncing
		io::BasicBuffer output_buffer(256000);
		if(this->block_settings.load_format) output_buffer.resize(256000 + this->header.getSampleNumber()*2);

		// Todo: in cases of non-diploid
		std::vector<core::GTObject> genotypes_unpermuted(this->header.getSampleNumber());
		for(U32 i = 0; i < this->header.getSampleNumber(); ++i){
			genotypes_unpermuted[i].alleles = new core::GTObjectAllele;
		}

		print_format_function print_format = &self_type::printFORMATDummy;
		if(this->block_settings.format_ID_list.size()) print_format = &self_type::printFORMATCustom;
		else if(block_settings.load_format) print_format = &self_type::printFORMATVCF;
		print_info_function   print_info   = &self_type::printINFOVCF;
		print_meta_function   print_meta   = &utility::to_vcf_string;
		print_filter_function print_filter = &self_type::printFILTER;

		// Cycling over loaded meta objects
		for(U32 p = 0; p < objects.meta->size(); ++p){
			if(filters.filter(objects, p) == false){
				//std::cerr << "failed filter" << std::endl;
				continue;
			}

			if(this->block_settings.custom_output_format)
				utility::to_vcf_string(output_buffer, '\t', (*objects.meta)[p], this->header, this->block_settings.custom_output_controller);
			else
				utility::to_vcf_string(output_buffer, '\t', (*objects.meta)[p], this->header);

			// Filter options
			(this->*print_filter)(output_buffer, p, objects);
			(this->*print_info)(output_buffer, '\t', p, objects);
			if(this->block_settings.annotate_extra) this->getGenotypeSummary(output_buffer, p, objects); // Todo: fixme
			(this->*print_format)(output_buffer, '\t', p, objects, genotypes_unpermuted);
			output_buffer += '\n';

			if(output_buffer.size() > 65536){
				std::cout.write(output_buffer.data(), output_buffer.size());
				output_buffer.reset();
				std::cout.flush();
			}
		}

		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
		std::cout.flush();

		return(objects.meta->size());
	}

	/**<
	 *
	 * @return
	 */
	const U32 outputBlockCustom(void) const{
		objects_type objects;
		this->loadObjects(objects);

		// Reserve memory for output buffer
		// This is much faster than writing directly to ostream because of syncing
		io::BasicBuffer output_buffer(256000 + this->header.getSampleNumber()*2);
		std::vector<core::GTObject> genotypes_unpermuted(this->header.getSampleNumber());

		// Todo: move to function
		//U32 info_match_limit = 1; // any match
		//info_match_limit = this->block_settings.info_list.size(); // all match

		// Function pointer to use
		print_format_function print_format = &self_type::printFORMATDummy;
		print_info_function   print_info   = &self_type::printINFODummy;
		print_meta_function   print_meta   = &utility::to_vcf_string;
		print_filter_function print_filter = &self_type::printFILTERDummy;

		if(block_settings.output_json) print_meta = &utility::to_json_string;

		if(block_settings.load_format || this->block.n_format_loaded){
			if(block_settings.output_json){
				print_format = &self_type::printFORMATCustomVectorJSON;
			} else {
				if(block_settings.output_format_vector) print_format = &self_type::printFORMATCustomVector;
				else print_format = &self_type::printFORMATCustom;
			}
		}

		if(block_settings.load_info || this->block.n_info_loaded){
			if(block_settings.output_json) print_info = &self_type::printINFOCustomJSON;
			else print_info = &self_type::printINFOCustom;
		}

		if(this->block_settings.custom_output_controller.show_filter){
			if(block_settings.output_json) print_filter = &self_type::printFILTERJSON;
			else print_filter = &self_type::printFILTERCustom;
		}

		U32 n_records_returned = 0;

		if(block_settings.output_json) output_buffer += "\"block\":[";
		for(U32 position = 0; position < objects.meta->size(); ++position){
			//if(info_keep[objects.meta->at(p).getInfoPatternID()] < info_match_limit)
			//	continue;
			if(block_settings.output_json){
				if(position != 0) output_buffer += ",\n";


				output_buffer += "{";
			}
			++n_records_returned;

			(*print_meta)(output_buffer, this->block_settings.custom_delimiter_char, (*objects.meta)[position], this->header, this->block_settings.custom_output_controller);
			(this->*print_filter)(output_buffer, position, objects);
			(this->*print_info)(output_buffer, this->block_settings.custom_delimiter_char, position, objects);
			(this->*print_format)(output_buffer, this->block_settings.custom_delimiter_char, position, objects, genotypes_unpermuted);

			if(block_settings.output_json) output_buffer += "}";
			else output_buffer += '\n';
			//output_buffer += "}";

			// Flush if buffer is large
			if(output_buffer.size() > 65536){
				std::cout.write(output_buffer.data(), output_buffer.size());
				output_buffer.reset();
				std::cout.flush();
			}
		}
		if(block_settings.output_json) output_buffer += "]";

		// Flush buffer
		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
		std::cout.flush();

		return(n_records_returned);
	}

	// Dummy functions as interfaces for function pointers
	inline void printFILTERDummy(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const{}
	inline void printFORMATDummy(buffer_type& buffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const{}
	inline void printINFODummy(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const{}

	// FILTER functions
	void printFILTER(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;
	void printFILTERCustom(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;
	void printFILTERJSON(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;

	// FORMAT functions
	void printFORMATVCF(buffer_type& buffer, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	void printFORMATVCF(buffer_type& buffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	void printFORMATCustom(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	void printFORMATCustomVector(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	void printFORMATCustomVectorJSON(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;

	// INFO functions
	void printINFOVCF(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;
	void printINFOVCF(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const;
	void printINFOCustom(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const;
	void printINFOCustomJSON(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const;

	// Calculations
	TACHYON_VARIANT_CLASSIFICATION_TYPE classifyVariant(const meta_entry_type& meta, const U32& allele) const;

	/**<
	 * Parse interval strings. These strings have to match the regular expression
	 * patterns
	 * YON_REGEX_CONTIG_ONLY, YON_REGEX_CONTIG_POSITION, or YON_REGEX_CONTIG_RANGE
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool parseIntervals(void){
		// Intervals pass expression tests
		if(this->getSettings().validateIntervalStrings() == false)
			return(false);

		// Data reference
		std::vector<std::string>& intervals = this->getSettings().interval_strings;

		// No intervals to parse
		if(intervals.size() == 0)
			return true;

		// Blocks to load
		std::vector<index_entry_type>& blocks_to_load = this->getSettings().getBlockList();

		// Parse each interval
		for(U32 i = 0; i < intervals.size(); ++i){
			intervals[i] = utility::remove_whitespace(intervals[i]);
			core::HeaderContig* contig = nullptr;

			if (std::regex_match (intervals[i], constants::YON_REGEX_CONTIG_ONLY )){
				std::cerr << "chromosome only" << std::endl;
				if(!this->header.getContig(intervals[i],contig)){
					std::cerr << "cant find contig: " << intervals[i] << std::endl;
					return(false);
				}

				std::cerr << "Parsed: " << intervals[i] << std::endl;

			} else if (std::regex_match (intervals[i], constants::YON_REGEX_CONTIG_POSITION )){
				std::cerr << "chromosome pos" << std::endl;
				std::vector<std::string> substrings = utility::split(intervals[i], ':');
				if(substrings[0].size() == 0 || substrings[1].size() == 0){
					std::cerr << "illegal form" << std::endl;
					return false;
				}

				if(!this->header.getContig(substrings[0],contig)){
					std::cerr << "cant find contig: " << substrings[0] << std::endl;
					return(false);
				}

				U64 position = atof(substrings[1].data());
				std::cerr << "Parsed: " << substrings[0] << "," << position << std::endl;

				std::vector<index_entry_type> target_blocks = this->index.findOverlap(contig->contigID, position);
				blocks_to_load.insert( blocks_to_load.end(), target_blocks.begin(), target_blocks.end() );

			} else if (std::regex_match (intervals[i], constants::YON_REGEX_CONTIG_RANGE )){
				std::cerr << "chromosome pos - pos" << std::endl;
				std::vector<std::string> substrings = utility::split(intervals[i], ':');
				if(substrings[0].size() == 0 || substrings[1].size() == 0){
					std::cerr << "illegal form" << std::endl;
					return false;
				}

				if(!this->header.getContig(substrings[0],contig)){
					std::cerr << "cant find contig: " << substrings[0] << std::endl;
					return(false);
				}

				std::vector<std::string> position_strings = utility::split(substrings[1], '-');
				if(position_strings[0].size() == 0 || position_strings[1].size() == 0){
					std::cerr << "illegal form" << std::endl;
					return false;
				}
				U64 position_from = atof(position_strings[0].data());
				U64 position_to   = atof(position_strings[1].data());
				if(position_from > position_to) std::swap(position_from, position_to);

				std::cerr << "Parsed: " << substrings[0] << "," << position_from << "," << position_to << std::endl;

				std::vector<index_entry_type> target_blocks = this->index.findOverlap(contig->contigID, position_from, position_to);
				blocks_to_load.insert( blocks_to_load.end(), target_blocks.begin(), target_blocks.end() );

			} else {
				std::cerr << utility::timestamp("ERROR") << "Uninterpretable interval string: " << intervals[i] << std::endl;
				return false;
			}
		}

		if(blocks_to_load.size() == 0)
			return true;

		// Dedupe
		std::sort(blocks_to_load.begin(), blocks_to_load.end());
		std::vector<index_entry_type> blocks_to_load_deduped;
		blocks_to_load_deduped.push_back(blocks_to_load[0]);
		for(U32 i = 1; i < blocks_to_load.size(); ++i){
			if(blocks_to_load_deduped.back() != blocks_to_load[i]){
				blocks_to_load_deduped.push_back(blocks_to_load[i]);
			}
		}

		// Debug
		//std::cerr << "size: " << blocks_to_load_deduped.size() << std::endl;
		//for(U32 i = 0; i < blocks_to_load_deduped.size(); ++i){
		//	std::cerr << i << "\t";
		//	blocks_to_load_deduped[i].print(std::cerr);
		//	std::cerr << std::endl;
		//}

		// Assign deduped vector to settings reference
		blocks_to_load = blocks_to_load_deduped;

		return true;
	}


	//<----------------- EXAMPLE FUNCTIONS -------------------------->


	U64 timings_meta(){
		containers::MetaContainer meta(this->block);
		buffer_type temp(meta.size() * 1000);

		for(U32 p = 0; p < meta.size(); ++p){
			utility::to_vcf_string(temp, this->block_settings.custom_delimiter_char, meta[p], this->header);
			//utility::to_vcf_string(std::cout, meta[p], this->header);
			temp += '\n';
			//if(temp.size() > 65536){
			//	std::cout.write(temp.data(), temp.size());
			//	temp.reset();
			//}
		}
		std::cout.write(temp.data(), temp.size());
		return(meta.size());
	}


	U64 iterate_genotypes(std::ostream& stream = std::cout){
		containers::MetaContainer meta(this->block);
		containers::GenotypeContainer gt(this->block, meta);

		for(U32 i = 0; i < gt.size(); ++i){
			// All of these functions are in relative terms very expensive!
			// Avoid using them unless you absolutely have to!
			// Vector of literal genotype representations (lower level)
			//std::vector<core::GTObject> objects     = gt[i].getLiteralObjects();
			// Vector of genotype objects (high level permuted)
			//std::vector<core::GTObject> objects_all = gt[i].getObjects(this->header.getSampleNumber());
			// Vector of genotype objects (high level unpermuted - original)
			std::vector<core::GTObject> objects_true = gt[i].getObjects(this->header.getSampleNumber(), this->block.ppa_manager);

			std::cout << (int)objects_true[i].alleles[0].allele << (objects_true[i].alleles[1].phase ? '/' : '|') << (int)objects_true[i].alleles[1].allele;
			for(U32 i = 1; i < objects_true.size(); ++i){
				std::cout << '\t' << (int)objects_true[i].alleles[0].allele << (objects_true[i].alleles[1].phase ? '/' : '|') << (int)objects_true[i].alleles[1].allele;
			}
			std::cout << std::endl;
		}
		return(gt.size());
	}

	U64 calculateIBS(math::SquareMatrix<double>& square, math::SquareMatrix<double>& square_temporary){
		algorithm::Timer timer;
		timer.Start();

		containers::MetaContainer meta(this->block);
		containers::GenotypeContainer gt(this->block, meta);
		for(U32 i = 0; i < gt.size(); ++i)
			gt[i].comparePairwise(square_temporary);

		//square /= (U64)2*this->header.getSampleNumber()*gt.size();
		square.addUpperTriagonal(square_temporary, this->block.ppa_manager);
		square_temporary.clear();

		// 2 * (Upper triagonal + diagonal) * number of variants
		const U64 updates = 2*((this->header.getSampleNumber()*this->header.getSampleNumber() - this->header.getSampleNumber())/2 + this->header.getSampleNumber()) * gt.size();
		std::cerr << utility::timestamp("DEBUG") << "Updates: " << utility::ToPrettyString(updates) << '\t' << timer.ElapsedString() << '\t' << utility::ToPrettyString((U64)((double)updates/timer.Elapsed().count())) << "/s" << std::endl;
		return((U64)2*this->header.getSampleNumber()*gt.size());
	}

	U64 getTiTVRatios(std::ostream& stream, std::vector<core::TsTvObject>& global){
		containers::MetaContainer meta(this->block);
		containers::GenotypeContainer gt(this->block, meta);

		std::vector<core::TsTvObject> objects(this->header.getSampleNumber());
		for(U32 i = 0; i < gt.size(); ++i)
			gt[i].getTsTv(objects);

		for(U32 i = 0; i < objects.size(); ++i)
			global[this->block.ppa_manager[i]] += objects[i];

		return(gt.size());
	}

	std::vector<double> calculateStrandBiasAlleles(const meta_entry_type& meta, const genotype_summary_type& genotype_summary, const bool phred_scale = true) const{
		std::vector<double> strand_bias_p_values(meta.n_alleles);
		double fisher_left_p, fisher_right_p, fisher_twosided_p;

		kt_fisher_exact(
		genotype_summary.vectorA_[2], // A: Allele on forward strand
		genotype_summary.vectorB_[2], // B: Allele on reverse strand
		genotype_summary.alleleCountA() - (genotype_summary.vectorA_[2]), // C: Not allele on forward strand
		genotype_summary.alleleCountB() - (genotype_summary.vectorB_[2]), // D: Not allele on reverse strand
		&fisher_left_p, &fisher_right_p, &fisher_twosided_p);

		if(phred_scale) strand_bias_p_values[0] = std::abs(-10 * log10(fisher_twosided_p));
		else strand_bias_p_values[0] = fisher_twosided_p;

		// If n_alleles = 2 then they are identical because of symmetry
		if(meta.n_alleles > 2){
			for(U32 p = 1; p < meta.n_alleles; ++p){
				kt_fisher_exact(
				genotype_summary.vectorA_[2+p], // A: Allele on forward strand
				genotype_summary.vectorB_[2+p], // B: Allele on reverse strand
				genotype_summary.alleleCountA() - (genotype_summary.vectorA_[2+p]), // C: Not allele on forward strand
				genotype_summary.alleleCountB() - (genotype_summary.vectorB_[2+p]), // D: Not allele on reverse strand
				&fisher_left_p, &fisher_right_p, &fisher_twosided_p);

				if(phred_scale) strand_bias_p_values[p] = std::abs(-10 * log10(fisher_twosided_p));
				else strand_bias_p_values[p] = fisher_twosided_p;
			}
		}
		return(strand_bias_p_values);
	}

	void getGenotypeSummary(buffer_type& buffer, const U32& position, objects_type& objects) const{
		if(this->block_settings.load_alleles == false || this->block_settings.load_genotypes_all == false || this->block_settings.load_controller == false || this->block_settings.load_set_membership == false){
			std::cerr << utility::timestamp("ERROR") << "Cannot run function without loading: SET-MEMBERSHIP, GT, REF or ALT, CONTIG or POSITION..." << std::endl;
			return;
		}

		if(buffer.back() != ';' && buffer.back() != '\t') buffer += ';';

		//objects_type objects;
		//this->loadObjects(objects);
		//U32 n_variants_parsed = 0;

		//for(U32 i = 0; i < objects.genotypes->size(); ++i){
			if(objects.meta->at(position).isDiploid() == false){
				std::cerr << "is not diploid" << std::endl;
				return;
			}

			// If set membership is -1 then calculate all fields
			// Set target FLAG set to all ones; update with actual values if they exist
			U16 target_flag_set = 65535;
			if(objects.meta->at(position).getInfoPatternID() != -1)
				target_flag_set = objects.additional_info_execute_flag_set[objects.meta->at(position).getInfoPatternID()];

			// Get genotype summary data
			objects.genotypes->at(position).getSummary(*objects.genotype_summary);
			std::vector<double> hwe_p = objects.genotype_summary->calculateHardyWeinberg(objects.meta->at(position));
			std::vector<double> af    = objects.genotype_summary->calculateAlleleFrequency(objects.meta->at(position));


			//utility::to_vcf_string(stream, this->block_settings.custom_delimiter_char, meta, this->header);

			if(target_flag_set & 1){
				std::vector<double> allele_bias = this->calculateStrandBiasAlleles(objects.meta->at(position), *objects.genotype_summary, true);
				buffer += "FS_A=";
				buffer.AddReadble(allele_bias[0]);
				for(U32 p = 1; p < allele_bias.size(); ++p){
					buffer += ',';
					buffer.AddReadble(allele_bias[p]);
				}
			}

			if(target_flag_set & 2){
				buffer += ";AN=";
				buffer.AddReadble(objects.genotype_summary->alleleCount());
			}

			if(target_flag_set & 4){
				if(objects.genotype_summary->vectorA_[1] + objects.genotype_summary->vectorB_[1]){
					buffer += ";NM=";
					buffer.AddReadble(objects.genotype_summary->vectorA_[1] + objects.genotype_summary->vectorB_[1]);
				}
			}

			if(target_flag_set & 8){
				if(objects.genotype_summary->vectorA_[0] + objects.genotype_summary->vectorB_[0]){
					buffer += ";NPM=";
					buffer.AddReadble(objects.genotype_summary->vectorA_[0] + objects.genotype_summary->vectorB_[0]);
				}
			}

			if(target_flag_set & 16){
				buffer += ";AC=";
				buffer.AddReadble(objects.genotype_summary->vectorA_[2] + objects.genotype_summary->vectorB_[2]);
				for(U32 p = 1; p < objects.meta->at(position).n_alleles; ++p){
					buffer += ",";
					buffer.AddReadble(objects.genotype_summary->vectorA_[2+p] + objects.genotype_summary->vectorB_[2+p]);
				}
			}

			if(target_flag_set & 32){
				buffer += ";AC_FWD=";
				buffer.AddReadble(objects.genotype_summary->vectorA_[2]);
				for(U32 p = 1; p < objects.meta->at(position).n_alleles; ++p){
					buffer += ",";
					buffer.AddReadble(objects.genotype_summary->vectorA_[2+p]);
				}
			}

			if(target_flag_set & 64){
				buffer += ";AC_REV=";
				buffer.AddReadble(objects.genotype_summary->vectorB_[2]);
				for(U32 p = 1; p < objects.meta->at(position).n_alleles; ++p){
					buffer += ",";
					buffer.AddReadble(objects.genotype_summary->vectorB_[2+p]);
				}
			}

			if(target_flag_set & 128){
				buffer += ";AF=";
				buffer.AddReadble(af[0]);
				for(U32 p = 1; p < af.size(); ++p){
					buffer += ",";
					buffer.AddReadble(af[p]);
				}
			}

			if(target_flag_set & 256){
				buffer += ";HWE_P=";
				buffer.AddReadble(hwe_p[0]);
				for(U32 p = 1; p < hwe_p.size(); ++p){
					buffer += ",";
					buffer.AddReadble(hwe_p[p]);
				}
			}

			if(target_flag_set & 512){
				// Classify
				buffer += ";VT=";
				buffer += TACHYON_VARIANT_CLASSIFICATION_STRING[this->classifyVariant(objects.meta->at(position), 1)];

				for(U32 p = 2; p < objects.meta->at(position).n_alleles; ++p){
					buffer += ',';
					buffer += TACHYON_VARIANT_CLASSIFICATION_STRING[this->classifyVariant(objects.meta->at(position), p)];
				}
			}

			if(target_flag_set & 1024){
				if(objects.meta->at(position).n_alleles != 2) buffer += ";MULTI_ALLELIC";
			}

			// Population inbreeding coefficient: F = (Hexp - Hobs) /Hexp
			if(target_flag_set & 2048){
				// Allele frequency of A
				const double p = ((double)2*objects.genotype_summary->matrix_[2][2] + objects.genotype_summary->matrix_[2][3] + objects.genotype_summary->matrix_[3][2]) / (2*objects.genotype_summary->genotypeCount());
				// Genotype frequency of heterozyotes
				const double pg = ((double)objects.genotype_summary->matrix_[2][3] + objects.genotype_summary->matrix_[3][2]) / objects.genotype_summary->genotypeCount();
				// Expected heterozygosity
				const double exp = 2*p*(1-p);
				// Population inbreeding coefficient: F
				const double f_pic = exp > 0 ? (exp-pg)/exp : 0;
				buffer += ";F_PIC=";
				buffer.AddReadble(f_pic);
			}

			//stream.put('\n');
			objects.genotype_summary->clear();
			//++n_variants_parsed;
		//}
		//return(n_variants_parsed);
	}

	U64 countVariants(std::ostream& stream = std::cout){
		containers::MetaContainer meta(this->block);
		return(meta.size());
	}

	U64 iterateMeta(std::ostream& stream = std::cout){
		containers::MetaContainer meta(this->block);
		containers::GenotypeContainer gt(this->block, meta);
		containers::GenotypeSummary gt_summary;
		for(U32 i = 0; i < gt.size(); ++i){
			// If there's > 5 alleles continue
			if(gt[i].getMeta().getNumberAlleles() >= 5) continue;
			// Calculate summary statistics
			//gt[i].getSummary(gt_summary);

			// Calculate total number of alt-alleles (allele 1, where 0 is ref)
			//std::cerr << gt_summary << '\n';
			gt_summary.clear(); // Recycle summary object
		}
		//std::cerr << std::endl;
		//std::cerr << gt.size() << std::endl;
		return(gt.size());
		//std::cerr << gt[0] << std::endl;;

		//return true;

		core::HeaderMapEntry* entry = nullptr;
		if(this->header.getInfoField("AF", entry)){
			containers::InfoContainer<double> it_i(this->block.info_containers[1]);
			//math::MathSummaryStatistics stats = it_i.getSummaryStatistics();
			//std::cerr << stats.n_total << '\t' << stats.mean << '\t' << stats.standard_deviation << '\t' << stats.min << "-" << stats.max << std::endl;
			for(U32 i = 0; i < it_i.size(); ++i){
				//if(it_i[i].size() < 3) continue;
				//it[i].toVCFString(stream, this->header, this->block.index_entry.contigID, this->block.index_entry.minPosition);

				//stream << (int)it_i[i][0];
				for(U32 j = 0; j < it_i[i].size(); ++j)
					stream << it_i[i][j] << ' ';
			}
			stream << '\n';
		}
		return(0);
	}

public:
	std::ifstream       stream;
	U64                 filesize;

	// Actual data
	block_entry_type    block;

	// Supportive objects
	block_settings_type block_settings;
	settings_type       settings;
	header_type         header;
	footer_type         footer;
	index_type          index;
	checksum_type       checksums;
	codec_manager_type  codec_manager;
	keychain_type       keychain;
};

}

#endif /* CORE_TACHYON_READER_H_ */
