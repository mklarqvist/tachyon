#ifndef CORE_BLOCKENTRYSETTINGS_H_
#define CORE_BLOCKENTRYSETTINGS_H_

#include <regex>
#include "../support/helpers.h"
#include "../core/header/variant_header.h"

namespace tachyon{

/**<
 * Supportive structure for Block
 */
struct SettingsMap{
	typedef SettingsMap                     self_type;
	typedef containers::DataContainerHeader header_type;

	SettingsMap() : iterator_index(0), target_stream_local(-1), offset(nullptr){}
	SettingsMap(const U32 iterator_index, const S32 target_stream_disk, const header_type* offset) : iterator_index(iterator_index), target_stream_local(target_stream_disk), offset(offset){}
	~SettingsMap(){}

	SettingsMap(const SettingsMap& other) : iterator_index(other.iterator_index), target_stream_local(other.target_stream_local), offset(other.offset){}
	SettingsMap(SettingsMap&& other) : iterator_index(other.iterator_index), target_stream_local(other.target_stream_local), offset(other.offset){}
	SettingsMap& operator=(const SettingsMap& other){
		this->iterator_index = other.iterator_index;
		this->target_stream_local = other.target_stream_local;
		this->offset = other.offset;
		return *this;
	}
	SettingsMap& operator=(SettingsMap&& other){
		if(this!=&other) // prevent self-move
		{
			this->iterator_index = other.iterator_index;
			this->target_stream_local = other.target_stream_local;
			this->offset = other.offset;
		}
		return *this;
	}

	inline bool operator<(const self_type& other) const{ return(this->offset->data_header.offset < other.offset->data_header.offset); }
	inline bool operator>(const self_type& other) const{ return(!((*this) < other)); }

	U32 iterator_index;        // Incrementor index
	S32 target_stream_local;   // Local target index
	const header_type* offset; // Header object of target data container
};

struct SettingsCustomOutput{
	SettingsCustomOutput() :
		show_contig(false),
		show_position(false),
		show_quality(false),
		show_names(false),
		show_ref(false),
		show_alt(false),
		show_filter(false),
		show_format_map(false)
	{}
	~SettingsCustomOutput(){}

	bool show_contig;
	bool show_position;
	bool show_quality;
	bool show_names;
	bool show_ref;
	bool show_alt;
	bool show_filter;
	bool show_format_map;
};

/**<
 * Settings
 */
struct DataBlockSettings{
	typedef DataBlockSettings   self_type;
	typedef SettingsMap         map_type;
	typedef core::VariantHeader header_type;

	DataBlockSettings() :
		show_vcf_header(true),
		load_contig(false),
		load_positons(false),
		load_controller(false),
		load_quality(false),
		load_names(false),
		load_alleles(false),
		load_set_membership(false),
		load_genotypes_all(false),
		load_genotypes_rle(false),
		load_genotypes_simple(false),
		load_genotypes_other(false),
		load_genotypes_support(false),
		load_ppa(false),
		load_info(false),
		load_format(false),
		construct_occ_table(false),
		custom_delimiter(false),
		custom_output_format(false),
		custom_delimiter_char('\t'),
		output_json(false),
		output_format_vector(false),
		annotate_extra(false)
	{}

	self_type& loadAll(const bool set = true){
		load_contig = set;
		load_positons = set;
		load_controller = set;
		load_quality = set;
		load_names = set;
		load_alleles = set;
		load_set_membership = set;
		load_genotypes_all = set;
		load_genotypes_rle = set;
		load_genotypes_simple = set;
		load_genotypes_other = set;
		load_genotypes_support = set;
		load_info = set;
		load_format = set;
		load_ppa = set;
		return(*this);
	}

	self_type& loadAllMeta(const bool set = true){
		load_contig = set;
		load_positons = set;
		load_controller = set;
		load_quality = set;
		load_names = set;
		load_alleles = set;
		return(*this);
	}

	self_type& loadAllFILTER(const bool set = true){
		load_set_membership = set;
		return(*this);
	}

	self_type& loadAllINFO(const bool set = true){
		load_info = set;
		return(*this);
	}

	self_type& loadINFO(const std::string& field_name){
		if(field_name.size() == 0) return(*this);
		this->load_contig = true;
		this->load_positons = true;
		this->load_set_membership = true;
		this->info_list.push_back(field_name);
		return(*this);
	}

	self_type& loadINFO(const U32 field_id){
		this->info_ID_list.push_back(field_id);
		this->load_contig = true;
		this->load_positons = true;
		this->load_set_membership = true;
		return(*this);
	}

	self_type& loadGenotypes(const bool set){
		load_genotypes_all = set;
		load_genotypes_rle = set;
		load_genotypes_simple = set;
		load_genotypes_other = set;
		load_genotypes_support = set;
		return(*this);
	}

	self_type& loadAllFORMAT(const bool set){
		this->load_ppa = set;
		this->loadGenotypes(set);
		this->load_format = set;
		this->load_contig = true;
		this->load_positons = true;
		this->load_set_membership = true;
		return(*this);
	}

	self_type& loadFORMAT(const std::string& field_name){
		if(field_name.size() == 0) return(*this);
		this->load_contig = true;
		this->load_positons = true;
		this->load_set_membership = true;
		if(field_name == "GT") this->loadGenotypes(true);
		this->format_list.push_back(field_name);
		return(*this);
	}

	self_type& loadFORMAT(const U32 field_id){
		this->load_contig = true;
		this->load_positons = true;
		this->load_set_membership = true;
		this->format_ID_list.push_back(field_id);
		return(*this);
	}

	self_type& setCustomDelimiter(const char delimiter){
		this->custom_delimiter = true;
		this->custom_delimiter_char = delimiter;
		return(*this);
	}

	bool parseCommandString(const std::vector<std::string>& command, const header_type& header, const bool customOutputFormat = false){
		this->custom_output_format = customOutputFormat;
		bool allGood = true;

		std::regex field_identifier_regex("^[A-Z_0-9]{1,}$");
		for(U32 i = 0; i < command.size(); ++i){
			std::vector<std::string> partitions = utility::split(command[i], ';');
			for(U32 p = 0; p < partitions.size(); ++p){
				partitions[p].erase(std::remove(partitions[p].begin(), partitions[p].end(), ' '), partitions[p].end()); // remove all spaces
				if(strncasecmp(partitions[p].data(), "INFO=", 5) == 0){
					std::vector<std::string> ind = utility::split(partitions[p].substr(5,command.size()-5), ',');
					for(U32 j = 0; j < ind.size(); ++j){
						ind[j] = utility::remove_excess_whitespace(ind[j]);
						//ind[j] = std::regex_replace(ind[j], std::regex("^ +| +$|( ) +"), std::string("$1")); // remove excess white space
						//std::transform(ind[j].begin(), ind[j].end(), ind[j].begin(), ::toupper); // transform to UPPERCASE
						if(std::regex_match(ind[j], field_identifier_regex)){
							const core::HeaderMapEntry* map = header.getInfoField(ind[j]);
							if(map == false){
								std::cerr << utility::timestamp("ERROR") << "Cannot find INFO field: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							this->loadINFO(ind[j]);
						} else {
							std::cerr << utility::timestamp("ERROR") << "Cannot find INFO field: " << ind[j] << " in string " << partitions[p] << std::endl;
							allGood = false;
						}
					}
				} else if(strncasecmp(partitions[p].data(), "INFO", 4) == 0){
					this->load_info = true;
					this->load_set_membership = true;
				} else if(strncasecmp(partitions[p].data(), "FORMAT=", 7) == 0){
					std::vector<std::string> ind = utility::split(partitions[p].substr(7,command.size()-7), ',');
					for(U32 j = 0; j < ind.size(); ++j){
						//ind[j] = std::regex_replace(ind[j], std::regex("^ +| +$|( ) +"), "$1"); // remove excess white space
						std::transform(ind[j].begin(), ind[j].end(), ind[j].begin(), ::toupper); // transform to UPPERCASE
						if(std::regex_match(ind[j], field_identifier_regex)){
							// Special case for genotypes
							if(strncasecmp(ind[j].data(), "GT", 2) == 0 && ind[j].size() == 2){
								this->load_contig = true;
								this->load_positons = true;
								this->load_controller = true;
								this->loadGenotypes(true);
								this->load_set_membership = true;
							} else if(strncasecmp(ind[j].data(), "GENOTYPES", 9) == 0 && ind[j].size() == 9){
								this->load_contig = true;
								this->load_positons = true;
								this->load_controller = true;
								this->loadGenotypes(true);
								this->load_set_membership = true;
							}
							// Any other FORMAT
							else {
								const core::HeaderMapEntry* map = header.getFormatField(ind[j]);
								if(map == false){
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
					this->custom_output_controller.show_contig = true;
					this->load_contig = true;
				} else if((strncasecmp(partitions[p].data(), "POSITION", 8) == 0 && partitions[p].length() == 8) ||
						  (strncasecmp(partitions[p].data(), "POS", 3) == 0 && partitions[p].length() == 3)){
					this->custom_output_controller.show_position = true;
					this->load_positons = true;
				} else if((strncasecmp(partitions[p].data(), "REF", 3) == 0 && partitions[p].length() == 3) ||
						  (strncasecmp(partitions[p].data(), "REFERENCE", 9) == 0 && partitions[p].length() == 9)){
					this->custom_output_controller.show_ref = true;
					this->load_alleles = true;
					this->load_controller = true;
				} else if((strncasecmp(partitions[p].data(), "ALT", 3) == 0 && partitions[p].length() == 3) ||
						  (strncasecmp(partitions[p].data(), "ALTERNATE", 9) == 0 && partitions[p].length() == 9)){
					this->custom_output_controller.show_alt = true;
					this->load_alleles = true;
					this->load_controller = true;
				} else if((strncasecmp(partitions[p].data(), "QUALITY", 7) == 0 && partitions[p].length() == 7) ||
						  (strncasecmp(partitions[p].data(), "QUAL", 4) == 0 && partitions[p].length() == 4)){
					this->custom_output_controller.show_quality = true;
					this->load_quality = true;
				} else if((strncasecmp(partitions[p].data(), "NAMES", 5) == 0 && partitions[p].length() == 5) ||
						  (strncasecmp(partitions[p].data(), "NAME", 4) == 0 && partitions[p].length() == 4)){
					this->custom_output_controller.show_names = true;
					this->load_names = true;
				} else if((strncasecmp(partitions[p].data(), "FILTERS", 7) == 0 && partitions[p].length() == 7) ||
						  (strncasecmp(partitions[p].data(), "FILTER", 6) == 0 && partitions[p].length() == 6)){
					this->custom_output_controller.show_filter = true;
					this->load_set_membership = true;
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
	bool load_contig;
	bool load_positons;
	bool load_controller;
	bool load_quality;
	bool load_names;
	bool load_alleles;
	bool load_set_membership;
	bool load_genotypes_all;
	bool load_genotypes_rle;
	bool load_genotypes_simple;
	bool load_genotypes_other;
	bool load_genotypes_support;
	bool load_ppa;
	bool load_info;
	bool load_format;
	bool construct_occ_table;
	bool custom_delimiter;
	bool custom_output_format;
	char custom_delimiter_char;

	bool output_json;
	bool output_format_vector;

	bool annotate_extra;

	SettingsCustomOutput custom_output_controller;

	std::vector<std::string> info_list;
	std::vector<std::string> format_list;
	std::vector<U32> info_ID_list;
	std::vector<U32> format_ID_list;

	//
	std::vector<map_type> load_info_ID_loaded;
	std::vector<map_type> load_format_ID_loaded;

	// blocks to load
	std::vector<U32> blocks_numbers;
};

}

#endif /* CORE_BLOCKENTRYSETTINGS_H_ */
