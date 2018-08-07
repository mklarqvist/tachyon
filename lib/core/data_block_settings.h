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
 * Load and display Settings for the basic variant data block
 */
struct DataBlockSettings{
public:
	typedef DataBlockSettings     self_type;
	typedef DataBlockSettingsPair pair_type;
	typedef VariantHeader         header_type;

public:
	DataBlockSettings();
	~DataBlockSettings() = default;

	self_type& LoadCore(bool display = true){
		for(U32 i = YON_BLK_CONTIG; i <= YON_BLK_ID_FILTER; ++i){
			this->LoadWrapper(true, i);
			this->DisplayWrapper(display, i);
		}
		return(*this);
	}
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

	bool Parse(const header_type& header);
	bool ParseCommandString(const std::vector<std::string>& command, const header_type& header, const bool customOutputFormat = false);
	bool ParseSettings(const header_type& header);

	inline self_type& LoadWrapper(bool set, const int field_id){
		const U32 offset = 1 << field_id;
		this->load_static &= ~(offset);
		this->load_static |= offset;
		return(*this);
	}

	inline self_type& DisplayWrapper(bool set, const int field_id){
		const U32 offset = 1 << field_id;
		this->display_static &= ~(offset);
		this->display_static |= offset;
		return(*this);
	}

	inline self_type& LoadDisplayWrapper(bool set, const int field_id){
		const U32 offset = 1 << field_id;
		this->LoadWrapper(set, offset);
		this->DisplayWrapper(set, offset);
		return(*this);
	}

public:
	bool show_vcf_header;

	// Special display
	bool display_ref;
	bool display_alt;
	bool display_filter;

	// Load/display pairs
	U32 load_static;
	U32 display_static;


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

	std::vector<U32> info_id_global;
	std::vector<U32> format_id_global;

	// blocks to load
	std::vector<U32> blocks_numbers;
};

}

#endif /* CORE_BLOCKENTRYSETTINGS_H_ */
