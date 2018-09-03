#ifndef CORE_BLOCKENTRYSETTINGS_H_
#define CORE_BLOCKENTRYSETTINGS_H_

#include <regex>
#include "support/helpers.h"
#include "core/header/variant_header.h"
#include "containers/components/data_container_header.h"

namespace tachyon{

/**<
 * Load and display Settings for the basic variant data block
 */
struct DataBlockSettings{
public:
	typedef DataBlockSettings     self_type;
	typedef VariantHeader         header_type;

public:
	DataBlockSettings();
	~DataBlockSettings() = default;

	self_type& LoadWrapper(bool set, const int field_bv);
	self_type& DisplayWrapper(bool set, const int field_bv);
	self_type& LoadDisplayWrapper(bool set, const int field_bv);
	self_type& LoadCore(const bool set = true);
	self_type& LoadAll(const bool set = true);
	self_type& LoadAllMeta(const bool set = true);
	self_type& LoadAllFilter(const bool set = true);
	self_type& LoadAllInfo(const bool set = true);
	self_type& LoadInfo(const std::string& field_name);
	self_type& LoadInfo(const uint32_t field_id);
	self_type& LoadGenotypes(const bool set);
	self_type& LoadPermutationArray(const bool set);
	self_type& LoadAllFormat(const bool set);
	self_type& LoadFormat(const std::string& field_name);
	self_type& LoadFormat(const uint32_t field_id);
	self_type& LoadMinimumVcf(const bool set = true);

	self_type& DisplayCore(const bool set = true);
	self_type& DisplayAll(const bool set = true);
	self_type& DisplayAllMeta(const bool set = true);
	self_type& DisplayAllFilter(const bool set = true);
	self_type& DisplayAllInfo(const bool set = true);
	self_type& DisplayGenotypes(const bool set);
	self_type& DisplayAllFormat(const bool set);
	self_type& DisplayMinimumVcf(const bool set = true);

	bool Parse(const header_type& header);
	bool ParseCommandString(const std::vector<std::string>& command, const header_type& header);

public:
	bool show_vcf_header;

	// Special display
	bool display_ref;
	bool display_alt;
	bool display_filter;

	// Load/display pairs
	uint32_t load_static;
	uint32_t display_static;

	bool construct_occ_table;

	bool annotate_extra;

	//SettingsCustomOutput custom_output_controller;

	std::vector<std::string> info_list;
	std::vector<std::string> format_list;

	std::vector<uint32_t> info_id_global;
	std::vector<uint32_t> format_id_global;

	// blocks to load
	std::vector<uint32_t> blocks_numbers;
};

}

#endif /* CORE_BLOCKENTRYSETTINGS_H_ */
