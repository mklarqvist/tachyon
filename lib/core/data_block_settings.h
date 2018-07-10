#ifndef CORE_BLOCKENTRYSETTINGS_H_
#define CORE_BLOCKENTRYSETTINGS_H_

#include <regex>
#include "support/helpers.h"
#include "core/header/variant_header.h"
#include "containers/components/data_container_header.h"

namespace tachyon{

struct DataBlockSettingsPair{
public:
	DataBlockSettingsPair() : load(false), display(false){}
	DataBlockSettingsPair(const bool load, const bool display) : load(load), display(display){}
	~DataBlockSettingsPair() = default;

	inline void operator()(const bool& load){ this->load = load; this->display = load; }
	inline void operator()(const bool& load, const bool& display){ this->load = load; this->display = display; }

public:
	bool load;
	bool display;
};

/**<
 * Settings
 */
struct DataBlockSettings{
public:
	typedef DataBlockSettings     self_type;
	typedef DataBlockSettingsPair pair_type;
	typedef core::VariantHeader   header_type;
	typedef core::HeaderMapEntry  header_map_type;

public:
	DataBlockSettings();
	~DataBlockSettings() = default;

	self_type& loadAll(const bool set = true);
	self_type& loadAllMeta(const bool set = true);
	self_type& loadAllFILTER(const bool set = true);
	self_type& loadAllINFO(const bool set = true);
	self_type& loadINFO(const std::string& field_name);
	self_type& loadINFO(const U32 field_id);
	self_type& loadGenotypes(const bool set);
	self_type& loadPermutationArray(const bool set);
	self_type& loadAllFORMAT(const bool set);
	self_type& loadFORMAT(const std::string& field_name);
	self_type& loadFORMAT(const U32 field_id);
	self_type& setCustomDelimiter(const char delimiter);
	bool parse(const header_type& header);
	bool parseCommandString(const std::vector<std::string>& command, const header_type& header, const bool customOutputFormat = false);

	/**<
	 *
	 * @param header
	 * @return
	 */
	bool parseSettings(const header_type& header);

public:
	bool show_vcf_header;

	// Special display
	bool display_ref;
	bool display_alt;
	bool display_filter;

	// Load/display pairs
	pair_type contig;
	pair_type positions;
	pair_type controller;
	pair_type quality;
	pair_type names;
	pair_type alleles;
	pair_type set_membership;
	pair_type genotypes_all;
	pair_type genotypes_rle;
	pair_type genotypes_simple;
	pair_type genotypes_other;
	pair_type genotypes_support;
	pair_type ppa;
	pair_type info_all;
	pair_type format_all;

	bool construct_occ_table;
	bool custom_delimiter;
	bool custom_output_format;
	char custom_delimiter_char;

	bool output_json;
	bool output_format_vector;

	bool annotate_extra;

	//SettingsCustomOutput custom_output_controller;

	std::vector<std::string> info_list;
	std::vector<std::string> format_list;

	std::vector<U32> info_ID_list;
	std::vector<U32> format_ID_list;

	// blocks to load
	std::vector<U32> blocks_numbers;
};

}

#endif /* CORE_BLOCKENTRYSETTINGS_H_ */
