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
	bool parseSettings(const header_type& header);

	inline void LoadDisplayStandard(bool set, const int field_id){
		this->LoadStandard(set, field_id);
		this->DisplayStandard(set, field_id);
	}
	inline void LoadStandard(bool set, const int field_id){ this->load_static |= set << field_id; }
	inline void DisplayStandard(bool set, const int field_id){ this->display_static |= set << field_id; }

	self_type& LoadWrapper(const int target){
		this->load_static &= ~(target);
		this->load_static |= target;
		return(*this);
	}
	inline self_type& LoadPPA(void){ return(this->LoadWrapper(YON_BLK_BV_PPA)); }
	inline self_type& LoadContigs(void){ return(this->LoadWrapper(YON_BLK_BV_CONTIG)); }
	inline self_type& LoadPositions(void){ return(this->LoadWrapper(YON_BLK_BV_POSITION)); }
	inline self_type& LoadRefAlt(void){ return(this->LoadWrapper(YON_BLK_BV_REFALT)); }
	inline self_type& LoadController(void){ return(this->LoadWrapper(YON_BLK_BV_CONTROLLER)); }
	inline self_type& LoadQuality(void){ return(this->LoadWrapper(YON_BLK_BV_QUALITY)); }
	inline self_type& LoadNames(void){ return(this->LoadWrapper(YON_BLK_BV_NAMES)); }
	inline self_type& LoadAlleles(void){ return(this->LoadWrapper(YON_BLK_BV_ALLELES)); }
	inline self_type& LoadBasic(void){
		this->LoadContigs();
		this->LoadPositions();
		this->LoadController();
		return(this->LoadSetMemberships());
	}

	inline self_type& LoadInfoId(void){ return(this->LoadWrapper(YON_BLK_BV_ID_INFO)); }
	inline self_type& LoadFormatId(void){ return(this->LoadWrapper(YON_BLK_BV_ID_FORMAT)); }
	inline self_type& LoadFilterId(void){ return(this->LoadWrapper(YON_BLK_BV_ID_FILTER)); }
	inline self_type& LoadSetMemberships(void){
		this->LoadInfoId();
		this->LoadFormatId();
		return(this->LoadFilterId());
	}

	inline self_type& LoadGt8(void){ return(this->LoadWrapper(YON_BLK_BV_GT_INT8)); }
	inline self_type& LoadGt16(void){ return(this->LoadWrapper(YON_BLK_BV_GT_INT16)); }
	inline self_type& LoadGt32(void){ return(this->LoadWrapper(YON_BLK_BV_GT_INT32)); }
	inline self_type& LoadGt64(void){ return(this->LoadWrapper(YON_BLK_BV_GT_INT64)); }
	inline self_type& LoadGtS8(void){ return(this->LoadWrapper(YON_BLK_BV_GT_S_INT8)); }
	inline self_type& LoadGtS16(void){ return(this->LoadWrapper(YON_BLK_BV_GT_S_INT16)); }
	inline self_type& LoadGtS32(void){ return(this->LoadWrapper(YON_BLK_BV_GT_S_INT32)); }
	inline self_type& LoadGtS64(void){ return(this->LoadWrapper(YON_BLK_BV_GT_S_INT64)); }
	inline self_type& LoadGtSupport(void){ return(this->LoadWrapper(YON_BLK_BV_GT_SUPPORT)); }
	inline self_type& LoadGenotypes(void){
		this->LoadGt8(); this->LoadGt16(); this->LoadGt32(); this->LoadGt64();
		this->LoadGtS8(); this->LoadGtS16(); this->LoadGtS32(); this->LoadGtS64();
		return(this->LoadGtSupport());
	}

	self_type& LoadDisplayWrapper(const int target){
		this->load_static &= ~(target);
		this->load_static |= target;
		this->display_static &= ~(target);
		this->display_static |= target;
		return(*this);
	}
	inline self_type& LoadDisplayPPA(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_PPA)); }
	inline self_type& LoadDisplayContigs(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_CONTIG)); }
	inline self_type& LoadDisplayPositions(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_POSITION)); }
	inline self_type& LoadDisplayRefAlt(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_REFALT)); }
	inline self_type& LoadDisplayController(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_CONTROLLER)); }
	inline self_type& LoadDisplayQuality(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_QUALITY)); }
	inline self_type& LoadDisplayNames(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_NAMES)); }
	inline self_type& LoadDisplayAlleles(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_ALLELES)); }
	inline self_type& LoadDisplayBasic(void){
		this->LoadDisplayContigs();
		this->LoadDisplayPositions();
		this->LoadDisplayController();
		return(this->LoadDisplaySetMemberships());
	}

	inline self_type& LoadDisplayInfoId(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_ID_INFO)); }
	inline self_type& LoadDisplayFormatId(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_ID_FORMAT)); }
	inline self_type& LoadDisplayFilterId(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_ID_FILTER)); }
	inline self_type& LoadDisplaySetMemberships(void){
		this->LoadDisplayInfoId();
		this->LoadDisplayFormatId();
		return(this->LoadDisplayFilterId());
	}

	inline self_type& LoadDisplayGt8(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_GT_INT8)); }
	inline self_type& LoadDisplayGt16(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_GT_INT16)); }
	inline self_type& LoadDisplayGt32(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_GT_INT32)); }
	inline self_type& LoadDisplayGt64(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_GT_INT64)); }
	inline self_type& LoadDisplayGtS8(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_GT_S_INT8)); }
	inline self_type& LoadDisplayGtS16(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_GT_S_INT16)); }
	inline self_type& LoadDisplayGtS32(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_GT_S_INT32)); }
	inline self_type& LoadDisplayGtS64(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_GT_S_INT64)); }
	inline self_type& LoadDisplayGtSupport(void){ return(this->LoadDisplayWrapper(YON_BLK_BV_GT_SUPPORT)); }
	inline self_type& LoadDisplayGenotypes(void){
		this->LoadDisplayGt8(); this->LoadDisplayGt16(); this->LoadDisplayGt32(); this->LoadDisplayGt64();
		this->LoadDisplayGtS8(); this->LoadDisplayGtS16(); this->LoadDisplayGtS32(); this->LoadDisplayGtS64();
		return(this->LoadDisplayGtSupport());
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

	std::vector<U32> info_ID_list;
	std::vector<U32> format_ID_list;

	// blocks to load
	std::vector<U32> blocks_numbers;
};

}

#endif /* CORE_BLOCKENTRYSETTINGS_H_ */
