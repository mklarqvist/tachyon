#ifndef CONTAINERS_VARIANT_READER_SETTINGS_H_
#define CONTAINERS_VARIANT_READER_SETTINGS_H_

#include "../index/index_entry.h"

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

}



#endif /* CONTAINERS_VARIANT_READER_SETTINGS_H_ */
