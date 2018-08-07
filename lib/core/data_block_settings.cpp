#include "containers/components/variant_block_footer.h"
#include "data_block_settings.h"

namespace tachyon{

DataBlockSettings::DataBlockSettings() :
	show_vcf_header(true),
	display_ref(false),
	display_alt(false),
	display_filter(false),
	load_static(std::numeric_limits<U32>::max()),
	display_static(std::numeric_limits<U32>::max()),
	construct_occ_table(false),
	custom_delimiter(false),
	custom_output_format(false),
	custom_delimiter_char('\t'),
	output_json(false),
	output_format_vector(false),
	annotate_extra(false)
{}

DataBlockSettings& DataBlockSettings::loadAll(const bool set){
	this->loadAllMeta(true);
	for(U32 i = 0; i < YON_BLK_N_STATIC; ++i)
		this->LoadWrapper(set, i);

	this->LoadWrapper(set, YON_BLK_N_STATIC); // all info
	this->LoadWrapper(set, YON_BLK_N_STATIC + 1); // all format

	this->display_alt = true;
	this->display_ref = true;
	this->display_filter = true;
	return(*this);
}

DataBlockSettings& DataBlockSettings::loadAllMeta(const bool set){
	this->LoadWrapper(set, YON_BLK_CONTIG);
	this->LoadWrapper(set, YON_BLK_POSITION);
	this->LoadWrapper(set, YON_BLK_CONTROLLER);
	this->LoadWrapper(set, YON_BLK_QUALITY);
	this->LoadWrapper(set, YON_BLK_NAMES);
	this->LoadWrapper(set, YON_BLK_ALLELES);

	this->display_alt = true;
	this->display_ref = true;
	this->display_filter = true;
	return(*this);
}

DataBlockSettings& DataBlockSettings::loadAllFILTER(const bool set){
	this->LoadWrapper(set, YON_BLK_ID_INFO);
	this->LoadWrapper(set, YON_BLK_ID_FORMAT);
	this->LoadWrapper(set, YON_BLK_ID_FILTER);
	this->LoadWrapper(set, YON_BLK_CONTROLLER);

	this->display_filter  = true;
	return(*this);
}

DataBlockSettings& DataBlockSettings::loadAllINFO(const bool set){
	this->LoadWrapper(set, YON_BLK_N_STATIC); // all info
	this->LoadCore(true);

	return(*this);
}

DataBlockSettings& DataBlockSettings::loadINFO(const std::string& field_name){
	if(field_name.size() == 0) return(*this);
	this->info_list.push_back(field_name);
	this->LoadWrapper(true, YON_BLK_CONTIG);
	this->LoadWrapper(true, YON_BLK_POSITION);
	this->LoadWrapper(true, YON_BLK_CONTROLLER);
	this->LoadWrapper(true, YON_BLK_ID_INFO);
	this->LoadWrapper(true, YON_BLK_ID_FORMAT);
	this->LoadWrapper(true, YON_BLK_ID_FILTER);
	return(*this);
}

DataBlockSettings& DataBlockSettings::loadINFO(const U32 field_id){
	this->info_id_global.push_back(field_id);
	this->LoadWrapper(true, YON_BLK_CONTIG);
	this->LoadWrapper(true, YON_BLK_POSITION);
	this->LoadWrapper(true, YON_BLK_CONTROLLER);
	this->LoadWrapper(true, YON_BLK_ID_INFO);
	this->LoadWrapper(true, YON_BLK_ID_FORMAT);
	this->LoadWrapper(true, YON_BLK_ID_FILTER);
	return(*this);
}

DataBlockSettings& DataBlockSettings::loadGenotypes(const bool set){
	this->LoadWrapper(true, YON_BLK_GT_INT8);
	this->LoadWrapper(true, YON_BLK_GT_INT16);
	this->LoadWrapper(true, YON_BLK_GT_INT32);
	this->LoadWrapper(true, YON_BLK_GT_INT64);
	this->LoadWrapper(true, YON_BLK_GT_S_INT8);
	this->LoadWrapper(true, YON_BLK_GT_S_INT16);
	this->LoadWrapper(true, YON_BLK_GT_S_INT32);
	this->LoadWrapper(true, YON_BLK_GT_S_INT64);
	this->LoadWrapper(true, YON_BLK_GT_N_INT8);
	this->LoadWrapper(true, YON_BLK_GT_N_INT16);
	this->LoadWrapper(true, YON_BLK_GT_N_INT32);
	this->LoadWrapper(true, YON_BLK_GT_N_INT64);
	this->LoadWrapper(true, YON_BLK_GT_SUPPORT);
	this->LoadWrapper(true, YON_BLK_BV_GT_PLOIDY);
	this->LoadWrapper(true, YON_BLK_PPA);
	return(*this);
}

DataBlockSettings& DataBlockSettings::loadPermutationArray(const bool set){
	this->LoadWrapper(true, YON_BLK_PPA);
	return(*this);
}

DataBlockSettings& DataBlockSettings::loadAllFORMAT(const bool set){
	this->loadGenotypes(set);
	this->LoadWrapper(set, YON_BLK_N_STATIC + 1); // all format
	this->LoadCore(true);
	return(*this);
}

DataBlockSettings& DataBlockSettings::loadFORMAT(const std::string& field_name){
	if(field_name.size() == 0) return(*this);
	this->LoadWrapper(true, YON_BLK_CONTIG);
	this->LoadWrapper(true, YON_BLK_POSITION);
	this->LoadWrapper(true, YON_BLK_CONTROLLER);
	this->LoadWrapper(true, YON_BLK_ID_INFO);
	this->LoadWrapper(true, YON_BLK_ID_FORMAT);
	this->LoadWrapper(true, YON_BLK_ID_FILTER);
	if(field_name == "GT") this->loadGenotypes(true);
	this->format_list.push_back(field_name);
	return(*this);
}

DataBlockSettings& DataBlockSettings::loadFORMAT(const U32 field_id){
	this->LoadWrapper(true, YON_BLK_CONTIG);
	this->LoadWrapper(true, YON_BLK_POSITION);
	this->LoadWrapper(true, YON_BLK_CONTROLLER);
	this->LoadWrapper(true, YON_BLK_ID_INFO);
	this->LoadWrapper(true, YON_BLK_ID_FORMAT);
	this->LoadWrapper(true, YON_BLK_ID_FILTER);
	this->format_id_global.push_back(field_id);
	return(*this);
}

DataBlockSettings& DataBlockSettings::setCustomDelimiter(const char delimiter){
	this->custom_delimiter = true;
	this->custom_delimiter_char = delimiter;
	return(*this);
}

bool DataBlockSettings::Parse(const header_type& header){
	std::regex field_identifier_regex("^[A-Za-z_0-9]{1,}$");

	for(U32 i = 0; i < this->info_list.size(); ++i){
		std::vector<std::string> ind = utility::split(this->info_list[i], ',');
		for(U32 j = 0; j < ind.size(); ++j){
			ind[j] = utility::remove_excess_whitespace(ind[j]);
			if(std::regex_match(ind[j], field_identifier_regex)){
				const io::VcfInfo* info = header.GetInfo(ind[j]);
				if(info == nullptr){
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

bool DataBlockSettings::ParseCommandString(const std::vector<std::string>& command, const header_type& header, const bool customOutputFormat){
	this->custom_output_format = customOutputFormat; // Todo
	bool allGood = true;

	this->display_static = 0;
	this->load_static = 0;

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
						const io::VcfInfo* info = header.GetInfo(ind[j]);
						if(info == nullptr){
							std::cerr << utility::timestamp("ERROR") << "Cannot find INFO field: " << ind[j] << " in string " << partitions[p] << std::endl;
							allGood = false;
							continue;
						}
						this->loadINFO(info->idx);
					} else {
						std::cerr << utility::timestamp("ERROR") << "Illegal field name: " << ind[j] << ". Must match \"[A-Za-z_0-9]\"..." << std::endl;
						allGood = false;
					}
				}
			} else if(strncasecmp(partitions[p].data(), "INFO", 4) == 0 && partitions[p].size() == 4){
				this->loadAllINFO(true);
			} else if(strncasecmp(partitions[p].data(), "FORMAT", 6) == 0 && partitions[p].size() == 6){
				this->loadAllFORMAT(true);

			} else if(strncasecmp(partitions[p].data(), "FORMAT=", 7) == 0){
				std::vector<std::string> ind = utility::split(partitions[p].substr(7,command.size()-7), ',');
				for(U32 j = 0; j < ind.size(); ++j){
					std::transform(ind[j].begin(), ind[j].end(), ind[j].begin(), ::toupper); // transform to UPPERCASE
					if(std::regex_match(ind[j], field_identifier_regex)){
						// Special case for genotypes
						if(strncasecmp(ind[j].data(), "GT", 2) == 0 && ind[j].size() == 2){
							this->LoadCore(true);
							this->loadGenotypes(true);

							const io::VcfFormat* fmt = header.GetFormat(ind[j]);
							if(fmt == nullptr){
								std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							this->loadFORMAT(fmt->idx);

						} else if(strncasecmp(ind[j].data(), "GENOTYPES", 9) == 0 && ind[j].size() == 9){
							this->LoadCore(true);
							this->loadGenotypes(true);

							const io::VcfFormat* fmt = header.GetFormat(ind[j]);
							if(fmt == nullptr){
								std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							this->loadFORMAT(fmt->idx);
						}
						// Any other FORMAT
						else {
							const io::VcfFormat* fmt = header.GetFormat(ind[j]);
							if(fmt == nullptr){
								std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							this->loadFORMAT(fmt->idx);
						}
					} else {
						std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
						allGood = false;
					}
				}

			} else if((strncasecmp(partitions[p].data(), "CONTIG", 6) == 0 && partitions[p].length() == 6) ||
					  (strncasecmp(partitions[p].data(), "CHROM", 5) == 0 && partitions[p].length() == 5)  ||
					  (strncasecmp(partitions[p].data(), "CHROMOSOME", 10) == 0 && partitions[p].length() == 10)){
				this->LoadWrapper(true, YON_BLK_CONTIG);
			} else if((strncasecmp(partitions[p].data(), "POSITION", 8) == 0 && partitions[p].length() == 8) ||
					  (strncasecmp(partitions[p].data(), "POS", 3) == 0 && partitions[p].length() == 3)){
				this->LoadWrapper(true, YON_BLK_POSITION);
			} else if((strncasecmp(partitions[p].data(), "REF", 3) == 0 && partitions[p].length() == 3) ||
					  (strncasecmp(partitions[p].data(), "REFERENCE", 9) == 0 && partitions[p].length() == 9)){
				this->LoadWrapper(true, YON_BLK_ALLELES);
				this->LoadWrapper(true, YON_BLK_REFALT);
				this->LoadWrapper(true, YON_BLK_CONTROLLER);
				this->display_ref = true;
			} else if((strncasecmp(partitions[p].data(), "ALT", 3) == 0 && partitions[p].length() == 3) ||
					  (strncasecmp(partitions[p].data(), "ALTERNATE", 9) == 0 && partitions[p].length() == 9)){
				this->LoadWrapper(true, YON_BLK_ALLELES);
				this->LoadWrapper(true, YON_BLK_REFALT);
				this->LoadWrapper(true, YON_BLK_CONTROLLER);
				this->display_alt = true;
			} else if((strncasecmp(partitions[p].data(), "QUALITY", 7) == 0 && partitions[p].length() == 7) ||
					  (strncasecmp(partitions[p].data(), "QUAL", 4) == 0 && partitions[p].length() == 4)){
				this->LoadWrapper(true, YON_BLK_QUALITY);
			} else if((strncasecmp(partitions[p].data(), "NAMES", 5) == 0 && partitions[p].length() == 5) ||
					  (strncasecmp(partitions[p].data(), "NAME", 4) == 0 && partitions[p].length() == 4)){
				this->LoadWrapper(true, YON_BLK_NAMES);
			} else if((strncasecmp(partitions[p].data(), "FILTERS", 7) == 0 && partitions[p].length() == 7) ||
					  (strncasecmp(partitions[p].data(), "FILTER", 6) == 0 && partitions[p].length() == 6)){
				this->LoadWrapper(true, YON_BLK_CONTROLLER);
				this->LoadWrapper(true, YON_BLK_ID_FILTER);
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

}
