#ifndef CORE_BLOCKENTRYSETTINGS_H_
#define CORE_BLOCKENTRYSETTINGS_H_

#include <regex>
#include "support/helpers.h"
#include "core/header/variant_header.h"

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
 * Supportive structure for Block
 */
struct SettingsMap{
	typedef SettingsMap                     self_type;
	typedef containers::DataContainerHeader header_type;

	SettingsMap() : load_order_index(-1), stream_id_local(-1), stream_id_global(-1), offset(nullptr){}
	SettingsMap(const U32 load_order_index, const S32 target_stream_disk, const header_type* offset) :
		load_order_index(load_order_index),
		stream_id_local(target_stream_disk),
		stream_id_global(offset->data_header.global_key),
		offset(offset)
	{}

	~SettingsMap(){}

	SettingsMap(const SettingsMap& other) :
		load_order_index(other.load_order_index),
		stream_id_local(other.stream_id_local),
		stream_id_global(other.stream_id_global),
		offset(other.offset)
	{}

	SettingsMap(SettingsMap&& other) :
		load_order_index(other.load_order_index),
		stream_id_local(other.stream_id_local),
		stream_id_global(other.stream_id_global),
		offset(other.offset)
	{}

	SettingsMap& operator=(const SettingsMap& other){
		this->load_order_index = other.load_order_index;
		this->stream_id_local = other.stream_id_local;
		this->stream_id_global = other.stream_id_global;
		this->offset = other.offset;
		return *this;
	}

	SettingsMap& operator=(SettingsMap&& other){
		if(this!=&other) // prevent self-move
		{
			this->load_order_index = other.load_order_index;
			this->stream_id_local = other.stream_id_local;
			this->stream_id_global = other.stream_id_global;
			this->offset = other.offset;
		}
		return *this;
	}

	inline bool operator<(const self_type& other) const{ return(this->offset->data_header.offset < other.offset->data_header.offset); }
	inline bool operator>(const self_type& other) const{ return(!((*this) < other)); }

	S32 load_order_index;  // Incrementor index
	S32 stream_id_local;   // Local target index
	S32 stream_id_global;  // Global target index
	const header_type* offset; // Header object of target data container
};

/**<
 * Settings
 */
struct DataBlockSettings{
public:
	typedef DataBlockSettings     self_type;
	typedef SettingsMap           map_type;
	typedef DataBlockSettingsPair pair_type;
	typedef core::VariantHeader   header_type;
	typedef core::HeaderMapEntry  header_map_type;

public:
	DataBlockSettings() :
		show_vcf_header(true),
		display_ref(false),
		display_alt(false),
		display_filter(false),
		construct_occ_table(false),
		custom_delimiter(false),
		custom_output_format(false),
		custom_delimiter_char('\t'),
		output_json(false),
		output_format_vector(false),
		annotate_extra(false)
	{}

	self_type& loadAll(const bool set = true){
		this->loadAllMeta(true);
		this->set_membership(set, set);
		this->genotypes_all(set, set);
		this->genotypes_rle(set, set);
		this->genotypes_simple(set, set);
		this->genotypes_other(set, set);
		this->genotypes_support(set, set);
		this->info_all(set, set);
		this->format_all(set, set);
		this->ppa(set, set);
		this->display_alt = true;
		this->display_ref = true;
		this->display_filter = true;
		return(*this);
	}

	self_type& loadAllMeta(const bool set = true){
		this->contig(set, set);
		this->positions(set, set);
		this->controller(set, set);
		this->quality(set, set);
		this->names(set, set);
		this->alleles(set, set);
		this->display_alt = true;
		this->display_ref = true;
		this->display_filter = true;
		return(*this);
	}

	self_type& loadAllFILTER(const bool set = true){
		this->set_membership(set, set);
		this->display_filter = true;
		return(*this);
	}

	self_type& loadAllINFO(const bool set = true){
		this->info_all(set, set);
		this->contig.load = set;
		this->positions.load = set;
		this->set_membership.load = set;
		return(*this);
	}

	self_type& loadINFO(const std::string& field_name){
		if(field_name.size() == 0) return(*this);
		this->contig.load = true;
		this->positions.load = true;
		this->set_membership.load = true;
		this->controller(true, false);
		this->info_list.push_back(field_name);
		return(*this);
	}

	self_type& loadINFO(const U32 field_id){
		this->info_ID_list.push_back(field_id);
		this->contig.load = true;
		this->positions.load = true;
		this->set_membership.load = true;
		this->controller.load = true;
		return(*this);
	}

	self_type& loadGenotypes(const bool set){
		this->genotypes_all(set, set);
		this->genotypes_rle(set, set);
		this->genotypes_simple(set, set);
		this->genotypes_other(set, set);
		this->genotypes_support.load = set;
		return(*this);
	}

	self_type& loadPermutationArray(const bool set){
		this->ppa.load = set;
		return(*this);
	}

	self_type& loadAllFORMAT(const bool set){
		this->ppa(set, set);
		this->loadGenotypes(set);
		this->format_all.load = set;
		this->contig.load = set;
		this->positions.load = set;
		this->set_membership.load = set;
		return(*this);
	}

	self_type& loadFORMAT(const std::string& field_name){
		if(field_name.size() == 0) return(*this);
		this->contig.load = true;
		this->positions.load = true;
		this->set_membership.load = true;
		if(field_name == "GT") this->loadGenotypes(true);
		this->format_list.push_back(field_name);
		return(*this);
	}

	self_type& loadFORMAT(const U32 field_id){
		this->contig.load = true;
		this->positions.load = true;
		this->set_membership.load = true;
		this->format_ID_list.push_back(field_id);
		return(*this);
	}

	self_type& setCustomDelimiter(const char delimiter){
		this->custom_delimiter = true;
		this->custom_delimiter_char = delimiter;
		return(*this);
	}

	bool parse(const header_type& header){
		std::regex field_identifier_regex("^[A-Za-z_0-9]{1,}$");

		for(U32 i = 0; i < this->info_list.size(); ++i){
			std::vector<std::string> ind = utility::split(this->info_list[i], ',');
			for(U32 j = 0; j < ind.size(); ++j){
				ind[j] = utility::remove_excess_whitespace(ind[j]);
				if(std::regex_match(ind[j], field_identifier_regex)){
					const header_map_type* map = header.getInfoField(ind[j]);
					if(map == nullptr){
						std::cerr << utility::timestamp("ERROR") << "Cannot find INFO field: " << ind[j] << " in string " << this->info_list[i] << std::endl;
						continue;
					}
					this->loadINFO(ind[j]);
				} else {
					std::cerr << utility::timestamp("ERROR") << "Illegal field name: " << ind[j] << ". Must match \"[A-Za-z_0-9]\"..." << std::endl;
					return(false);
				}
				this->loadINFO(ind[j]);
			}
		}

		for(U32 i = 0; i < this->format_list.size(); ++i){

		}

		return true;
	}

	bool parseCommandString(const std::vector<std::string>& command, const header_type& header, const bool customOutputFormat = false){
		this->custom_output_format = customOutputFormat; // Todo
		bool allGood = true;

		std::regex field_identifier_regex("^[A-Za-z_0-9]{1,}$");
		for(U32 i = 0; i < command.size(); ++i){
			std::vector<std::string> partitions = utility::split(command[i], ';');
			for(U32 p = 0; p < partitions.size(); ++p){
				partitions[p].erase(std::remove(partitions[p].begin(), partitions[p].end(), ' '), partitions[p].end()); // remove all spaces
				if(strncasecmp(partitions[p].data(), "INFO=", 5) == 0){
					std::vector<std::string> ind = utility::split(partitions[p].substr(5,command.size()-5), ',');
					for(U32 j = 0; j < ind.size(); ++j){
						ind[j] = utility::remove_excess_whitespace(ind[j]);
						if(std::regex_match(ind[j], field_identifier_regex)){
							const header_map_type* map = header.getInfoField(ind[j]);
							if(map == nullptr){
								std::cerr << utility::timestamp("ERROR") << "Cannot find INFO field: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							this->loadINFO(ind[j]);
						} else {
							std::cerr << utility::timestamp("ERROR") << "Illegal field name: " << ind[j] << ". Must match \"[A-Za-z_0-9]\"..." << std::endl;
							allGood = false;
						}
					}
				} else if(strncasecmp(partitions[p].data(), "INFO", 4) == 0 && partitions[p].size() == 4){
					this->loadAllINFO(true);
				} else if(strncasecmp(partitions[p].data(), "FORMAT", 6) == 0 && partitions[p].size() == 6){
					this->format_all(true, true);
					this->set_membership(true, true);
					this->loadGenotypes(true);
					this->set_membership.load = true;
					this->positions.load = true;
					this->contig.load = true;
					this->ppa.load = true;
					this->controller.load = true;

				} else if(strncasecmp(partitions[p].data(), "FORMAT=", 7) == 0){
					std::vector<std::string> ind = utility::split(partitions[p].substr(7,command.size()-7), ',');
					for(U32 j = 0; j < ind.size(); ++j){
						std::transform(ind[j].begin(), ind[j].end(), ind[j].begin(), ::toupper); // transform to UPPERCASE
						if(std::regex_match(ind[j], field_identifier_regex)){
							// Special case for genotypes
							if(strncasecmp(ind[j].data(), "GT", 2) == 0 && ind[j].size() == 2){
								this->contig.load = true;
								this->positions.load = true;
								this->controller.load = true;
								this->loadGenotypes(true);
								this->set_membership.load = true;
							} else if(strncasecmp(ind[j].data(), "GENOTYPES", 9) == 0 && ind[j].size() == 9){
								this->contig.load = true;
								this->positions.load = true;
								this->controller.load = true;
								this->loadGenotypes(true);
								this->set_membership.load = true;
							}
							// Any other FORMAT
							else {
								const header_map_type* map = header.getFormatField(ind[j]);
								if(map == nullptr){
									std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
									allGood = false;
									continue;
								}
								this->loadFORMAT(ind[j]);
							}
						} else {
							std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
							allGood = false;
						}
					}

				} else if((strncasecmp(partitions[p].data(), "CONTIG", 6) == 0 && partitions[p].length() == 6) ||
						  (strncasecmp(partitions[p].data(), "CHROM", 5) == 0 && partitions[p].length() == 5)  ||
						  (strncasecmp(partitions[p].data(), "CHROMOSOME", 10) == 0 && partitions[p].length() == 10)){
					this->contig(true, true);
				} else if((strncasecmp(partitions[p].data(), "POSITION", 8) == 0 && partitions[p].length() == 8) ||
						  (strncasecmp(partitions[p].data(), "POS", 3) == 0 && partitions[p].length() == 3)){
					this->positions(true, true);
				} else if((strncasecmp(partitions[p].data(), "REF", 3) == 0 && partitions[p].length() == 3) ||
						  (strncasecmp(partitions[p].data(), "REFERENCE", 9) == 0 && partitions[p].length() == 9)){
					this->alleles(true, true);
					this->controller(true, true);
					this->display_ref = true;
				} else if((strncasecmp(partitions[p].data(), "ALT", 3) == 0 && partitions[p].length() == 3) ||
						  (strncasecmp(partitions[p].data(), "ALTERNATE", 9) == 0 && partitions[p].length() == 9)){
					this->alleles(true, true);
					this->controller(true, true);
					this->display_alt = true;
				} else if((strncasecmp(partitions[p].data(), "QUALITY", 7) == 0 && partitions[p].length() == 7) ||
						  (strncasecmp(partitions[p].data(), "QUAL", 4) == 0 && partitions[p].length() == 4)){
					this->quality(true, true);
				} else if((strncasecmp(partitions[p].data(), "NAMES", 5) == 0 && partitions[p].length() == 5) ||
						  (strncasecmp(partitions[p].data(), "NAME", 4) == 0 && partitions[p].length() == 4)){
					this->names(true, true);
				} else if((strncasecmp(partitions[p].data(), "FILTERS", 7) == 0 && partitions[p].length() == 7) ||
						  (strncasecmp(partitions[p].data(), "FILTER", 6) == 0 && partitions[p].length() == 6)){
					this->set_membership(true, true);
					this->display_filter = true;
				} else {
					std::cerr << utility::timestamp("ERROR") << "Unknown pattern: " << partitions[p] << std::endl;
					allGood = false;
				}
			}
		}

		if(allGood == false) return false;
		return true;
	}

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

	//
	std::vector<map_type> info_map;
	std::vector<map_type> format_map;
	std::vector<map_type> load_info_ID_loaded;
	std::vector<map_type> load_format_ID_loaded;

	// blocks to load
	std::vector<U32> blocks_numbers;
};

}

#endif /* CORE_BLOCKENTRYSETTINGS_H_ */
