#ifndef CONTAINERS_VARIANT_READER_SETTINGS_H_
#define CONTAINERS_VARIANT_READER_SETTINGS_H_

#include "third_party/intervalTree.h"
#include "index/index_entry.h"

namespace tachyon{

struct VariantReaderSettings{
public:
	typedef VariantReaderSettings          self_type;
	typedef index::IndexEntry              index_entry_type;
	typedef algorithm::Interval<U32, S64>  interval_type;
	typedef algorithm::IntervalTree<interval_type, S64> interval_tree_type;

public:
	VariantReaderSettings() :
		drop_format(false),
		header_only(false),
		show_header(true),
		annotate_genotypes(false),
		output("-"),
		use_htslib(false)
	{}

	~VariantReaderSettings() = default;

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

public:
	bool drop_format; // drop FORMAT fields
	bool header_only; // show only the VCF header
	bool show_header; // show the VCF header
	bool annotate_genotypes;
	bool use_htslib;
	std::string input;
	std::string output;
	std::string keychain_file;
	std::string output_type;
	std::string htslib_output_type;
	std::string sample_names_file;
	std::vector<std::string> sample_names;
};

}



#endif /* CONTAINERS_VARIANT_READER_SETTINGS_H_ */
