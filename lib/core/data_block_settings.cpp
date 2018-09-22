#include "data_block_settings.h"

namespace tachyon{

DataBlockSettings::DataBlockSettings() :
	show_vcf_header(true),
	display_ref(true),
	display_alt(true),
	display_filter(true),
	load_static(0),
	display_static(std::numeric_limits<uint32_t>::max()),
	construct_occ_table(false),
	annotate_extra(false)
{}

DataBlockSettings& DataBlockSettings::LoadWrapper(bool set, const int field_bv){
	this->load_static &= ~(field_bv);
	if(set) this->load_static |= field_bv;
	return(*this);
}

DataBlockSettings& DataBlockSettings::DisplayWrapper(bool set, const int field_bv){
	this->display_static &= ~(field_bv);
	if(set) this->display_static |= field_bv;
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadDisplayWrapper(bool set, const int field_bv){
	this->LoadWrapper(set, field_bv);
	this->DisplayWrapper(set, field_bv);
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadCore(const bool set){
	for(uint32_t i = YON_BLK_CONTIG; i <= YON_BLK_ID_FILTER; ++i){
		const uint32_t bv = 1 << i;
		this->LoadWrapper(set, bv);
	}
	return(*this);
}

DataBlockSettings& DataBlockSettings::DisplayCore(const bool set){
	for(uint32_t i = YON_BLK_CONTIG; i <= YON_BLK_ID_FILTER; ++i){
		const uint32_t bv = 1 << i;
		this->DisplayWrapper(set, bv);
	}
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadAll(const bool set){
	if(set){
		this->load_static = std::numeric_limits<uint32_t>::max();
	} else {
		this->load_static = 0;
	}
	return(*this);
}

DataBlockSettings& DataBlockSettings::DisplayAll(const bool set){
	if(set){
		this->display_static = std::numeric_limits<uint32_t>::max();
	} else {
		this->display_static = 0;
	}
	this->display_alt = set;
	this->display_ref = set;
	this->display_filter = set;
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadAllMeta(const bool set){
	for(uint32_t i = YON_BLK_CONTIG; i <= YON_BLK_ID_FILTER; ++i){
		const uint32_t bv = 1 << i;
		this->LoadWrapper(set, i);
	}
	return(*this);
}

DataBlockSettings& DataBlockSettings::DisplayAllMeta(const bool set){
	for(uint32_t i = YON_BLK_CONTIG; i <= YON_BLK_ID_FILTER; ++i){
		const uint32_t bv = 1 << i;
		this->DisplayWrapper(set, i);
	}

	this->display_alt = set;
	this->display_ref = set;
	this->display_filter = set;
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadAllFilter(const bool set){
	this->LoadWrapper(set, YON_BLK_BV_ID_INFO);
	this->LoadWrapper(set, YON_BLK_BV_ID_FORMAT);
	this->LoadWrapper(set, YON_BLK_BV_ID_FILTER);
	this->LoadWrapper(set, YON_BLK_BV_CONTROLLER);

	return(*this);
}

DataBlockSettings& DataBlockSettings::DisplayAllFilter(const bool set){
	this->DisplayWrapper(set, YON_BLK_BV_ID_INFO);
	this->DisplayWrapper(set, YON_BLK_BV_ID_FORMAT);
	this->DisplayWrapper(set, YON_BLK_BV_ID_FILTER);
	this->DisplayWrapper(set, YON_BLK_BV_CONTROLLER);

	this->display_filter  = true;
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadAllInfo(const bool set){
	this->LoadWrapper(set, YON_BLK_BV_INFO); // all info
	if(set) this->LoadMinimumVcf(true);

	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadInfo(const std::string& field_name){
	if(field_name.size() == 0) return(*this);
	this->info_list.push_back(field_name);
	this->LoadMinimumVcf(true);
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadInfo(const uint32_t field_id){
	this->info_id_global.push_back(field_id);
	this->LoadMinimumVcf(true);
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadGenotypes(const bool set){
	this->LoadWrapper(set, YON_BLK_BV_GT);
	this->LoadWrapper(set, YON_BLK_BV_GT_SUPPORT);
	this->LoadWrapper(set, YON_BLK_BV_GT_PLOIDY);
	this->LoadWrapper(set, YON_BLK_BV_PPA);
	return(*this);
}

DataBlockSettings& DataBlockSettings::DisplayGenotypes(const bool set){
	this->DisplayWrapper(set, YON_BLK_BV_GT);
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadPermutationArray(const bool set){
	this->LoadWrapper(set, YON_BLK_BV_PPA);
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadAllFormat(const bool set){
	this->LoadGenotypes(set);
	this->LoadWrapper(set, YON_BLK_BV_FORMAT); // all format
	if(set) this->LoadCore(set);
	return(*this);
}

DataBlockSettings& DataBlockSettings::DisplayAllFormat(const bool set){
	this->DisplayGenotypes(set);
	this->DisplayWrapper(set, YON_BLK_BV_FORMAT); // all format
	return(*this);
}

DataBlockSettings& DataBlockSettings::DisplayAllInfo(const bool set){
	this->DisplayWrapper(set, YON_BLK_BV_INFO);
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadFormat(const std::string& field_name){
	if(field_name.size() == 0) return(*this);
	this->LoadMinimumVcf(true);
	if(field_name == "GT") this->LoadGenotypes(true);
	this->format_list.push_back(field_name);
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadFormat(const uint32_t field_id){
	this->LoadMinimumVcf(true);
	this->format_id_global.push_back(field_id);
	return(*this);
}

DataBlockSettings& DataBlockSettings::LoadMinimumVcf(const bool set){
	this->LoadWrapper(set, YON_BLK_BV_CONTIG);
	this->LoadWrapper(set, YON_BLK_BV_POSITION);
	this->LoadWrapper(set, YON_BLK_BV_CONTROLLER);
	this->LoadWrapper(set, YON_BLK_BV_ID_INFO);
	this->LoadWrapper(set, YON_BLK_BV_ID_FORMAT);
	this->LoadWrapper(set, YON_BLK_BV_ID_FILTER);
	this->LoadWrapper(set, YON_BLK_BV_ALLELES);
	this->LoadWrapper(set, YON_BLK_BV_REFALT);
	return(*this);
}

DataBlockSettings& DataBlockSettings::DisplayMinimumVcf(const bool set){
	this->DisplayWrapper(set, YON_BLK_BV_CONTIG);
	this->DisplayWrapper(set, YON_BLK_BV_POSITION);
	return(*this);
}

bool DataBlockSettings::Parse(const header_type& header){
	std::regex field_identifier_regex("^[A-Za-z_0-9]{1,}$");

	for(uint32_t i = 0; i < this->info_list.size(); ++i){
		std::vector<std::string> ind = utility::split(this->info_list[i], ',');
		for(uint32_t j = 0; j < ind.size(); ++j){
			ind[j] = utility::remove_excess_whitespace(ind[j]);
			if(std::regex_match(ind[j], field_identifier_regex)){
				const VcfInfo* info = header.GetInfo(ind[j]);
				if(info == nullptr){
					std::cerr << utility::timestamp("ERROR") << "Cannot find INFO field: " << ind[j] << " in string " << this->info_list[i] << std::endl;
					continue;
				}
				this->LoadInfo(ind[j]);
			} else {
				std::cerr << utility::timestamp("ERROR") << "Illegal field name: " << ind[j] << ". Must match \"[A-Za-z_0-9]\"..." << std::endl;
				return(false);
			}
			this->LoadInfo(ind[j]);
		}
	}

	for(uint32_t i = 0; i < this->format_list.size(); ++i){

	}

	return true;
}

bool DataBlockSettings::ParseCommandString(const std::vector<std::string>& command, const header_type& header){
	bool allGood = true;

	this->display_static = 0;
	this->load_static = 0;

	std::regex field_identifier_regex("^[A-Za-z_0-9]{1,}$");
	for(uint32_t i = 0; i < command.size(); ++i){
		std::vector<std::string> partitions = utility::split(command[i], ';');
		for(uint32_t p = 0; p < partitions.size(); ++p){
			partitions[p].erase(std::remove(partitions[p].begin(), partitions[p].end(), ' '), partitions[p].end()); // remove all spaces
			if(strncasecmp(partitions[p].data(), "INFO=", 5) == 0){
				std::vector<std::string> ind = utility::split(partitions[p].substr(5,command.size()-5), ',');
				for(uint32_t j = 0; j < ind.size(); ++j){
					ind[j] = utility::remove_excess_whitespace(ind[j]);
					if(std::regex_match(ind[j], field_identifier_regex)){
						const VcfInfo* info = header.GetInfo(ind[j]);
						if(info == nullptr){
							std::cerr << utility::timestamp("ERROR") << "Cannot find INFO field: " << ind[j] << " in string " << partitions[p] << std::endl;
							allGood = false;
							continue;
						}
						this->LoadInfo(info->idx);
						this->DisplayAllInfo(true);
					} else {
						std::cerr << utility::timestamp("ERROR") << "Illegal field name: " << ind[j] << ". Must match \"[A-Za-z_0-9]\"..." << std::endl;
						allGood = false;
					}
				}
			} else if(strncasecmp(partitions[p].data(), "INFO", 4) == 0 && partitions[p].size() == 4){
				this->LoadAllInfo(true);
				this->DisplayAllInfo(true);
			} else if(strncasecmp(partitions[p].data(), "FORMAT", 6) == 0 && partitions[p].size() == 6){
				this->LoadAllFormat(true);
				this->DisplayAllFormat(true);

			} else if(strncasecmp(partitions[p].data(), "FORMAT=", 7) == 0){
				std::vector<std::string> ind = utility::split(partitions[p].substr(7,command.size()-7), ',');
				for(uint32_t j = 0; j < ind.size(); ++j){
					std::transform(ind[j].begin(), ind[j].end(), ind[j].begin(), ::toupper); // transform to UPPERCASE
					if(std::regex_match(ind[j], field_identifier_regex)){
						// Special case for genotypes
						if(strncasecmp(ind[j].data(), "GT", 2) == 0 && ind[j].size() == 2){
							this->LoadMinimumVcf(true);
							this->LoadGenotypes(true);

							const VcfFormat* fmt = header.GetFormat(ind[j]);
							if(fmt == nullptr){
								std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							this->LoadFormat(fmt->idx);
							this->DisplayAllFormat(true);

						} else if(strncasecmp(ind[j].data(), "GENOTYPES", 9) == 0 && ind[j].size() == 9){
							this->LoadMinimumVcf(true);
							this->LoadGenotypes(true);
							this->DisplayAllFormat(true);

							const VcfFormat* fmt = header.GetFormat(ind[j]);
							if(fmt == nullptr){
								std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							this->LoadFormat(fmt->idx);
							this->DisplayAllFormat(true);
						}
						// Any other FORMAT
						else {
							const VcfFormat* fmt = header.GetFormat(ind[j]);
							if(fmt == nullptr){
								std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							this->LoadFormat(fmt->idx);
							this->DisplayAllFormat(true);
						}
					} else {
						std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
						allGood = false;
					}
				}

			} else if((strncasecmp(partitions[p].data(), "CONTIG", 6) == 0 && partitions[p].length() == 6) ||
					  (strncasecmp(partitions[p].data(), "CHROM", 5) == 0 && partitions[p].length() == 5)  ||
					  (strncasecmp(partitions[p].data(), "CHROMOSOME", 10) == 0 && partitions[p].length() == 10)){
				this->LoadWrapper(true, YON_BLK_BV_CONTIG);
				this->DisplayWrapper(true, YON_BLK_BV_CONTIG);
			} else if((strncasecmp(partitions[p].data(), "POSITION", 8) == 0 && partitions[p].length() == 8) ||
					  (strncasecmp(partitions[p].data(), "POS", 3) == 0 && partitions[p].length() == 3)){
				this->LoadWrapper(true, YON_BLK_BV_POSITION);
			} else if((strncasecmp(partitions[p].data(), "REF", 3) == 0 && partitions[p].length() == 3) ||
					  (strncasecmp(partitions[p].data(), "REFERENCE", 9) == 0 && partitions[p].length() == 9)){
				this->LoadWrapper(true, YON_BLK_BV_ALLELES);
				this->LoadWrapper(true, YON_BLK_BV_REFALT);
				this->LoadWrapper(true, YON_BLK_BV_CONTROLLER);
				this->DisplayWrapper(true, YON_BLK_BV_ALLELES);
				this->DisplayWrapper(true, YON_BLK_BV_REFALT);
				this->DisplayWrapper(true, YON_BLK_BV_CONTROLLER);
				this->display_ref = true;
			} else if((strncasecmp(partitions[p].data(), "ALT", 3) == 0 && partitions[p].length() == 3) ||
					  (strncasecmp(partitions[p].data(), "ALTERNATE", 9) == 0 && partitions[p].length() == 9)){
				this->LoadWrapper(true, YON_BLK_BV_ALLELES);
				this->LoadWrapper(true, YON_BLK_BV_REFALT);
				this->LoadWrapper(true, YON_BLK_BV_CONTROLLER);
				this->DisplayWrapper(true, YON_BLK_BV_ALLELES);
				this->DisplayWrapper(true, YON_BLK_BV_REFALT);
				this->DisplayWrapper(true, YON_BLK_BV_CONTROLLER);
				this->display_alt = true;
			} else if((strncasecmp(partitions[p].data(), "QUALITY", 7) == 0 && partitions[p].length() == 7) ||
					  (strncasecmp(partitions[p].data(), "QUAL", 4) == 0 && partitions[p].length() == 4)){
				this->LoadWrapper(true, YON_BLK_BV_QUALITY);
				this->DisplayWrapper(true, YON_BLK_BV_QUALITY);
			} else if((strncasecmp(partitions[p].data(), "NAMES", 5) == 0 && partitions[p].length() == 5) ||
					  (strncasecmp(partitions[p].data(), "NAME", 4) == 0 && partitions[p].length() == 4)){
				this->LoadWrapper(true, YON_BLK_BV_NAMES);
				this->DisplayWrapper(true, YON_BLK_BV_NAMES);
			} else if((strncasecmp(partitions[p].data(), "FILTERS", 7) == 0 && partitions[p].length() == 7) ||
					  (strncasecmp(partitions[p].data(), "FILTER", 6) == 0 && partitions[p].length() == 6)){
				this->LoadWrapper(true, YON_BLK_BV_CONTROLLER);
				this->LoadWrapper(true, YON_BLK_BV_ID_FILTER);
				this->DisplayWrapper(true, YON_BLK_BV_CONTROLLER);
				this->DisplayWrapper(true, YON_BLK_BV_ID_FILTER);
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
