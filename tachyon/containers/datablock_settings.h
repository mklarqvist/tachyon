#ifndef CORE_BLOCKENTRYSETTINGS_H_
#define CORE_BLOCKENTRYSETTINGS_H_

namespace tachyon{
namespace core{

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
	typedef DataBlockSettings self_type;
	typedef SettingsMap       map_type;

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
		custom_delimiter_char('\t')
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

	self_type& loadFORMAT(const bool set){
		this->load_ppa = set;
		this->loadGenotypes(set);
		this->load_format = set;
		return(*this);
	}

	self_type& loadFORMAT(const std::string& field_name){
		if(field_name.size() == 0) return(*this);
		this->format_list.push_back(field_name);
		this->loadFORMAT(true);
		return(*this);
	}

	self_type& loadFORMAT(const U32 field_id){
		this->format_ID_list.push_back(field_id);
		this->loadFORMAT(true);
		return(*this);
	}

	self_type& setCustomDelimiter(const char delimiter){
		this->custom_delimiter = true;
		this->custom_delimiter_char = delimiter;
		return(*this);
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

	SettingsCustomOutput custom_output_controller;

	std::vector<std::string> samples_list;
	std::vector<std::string> info_list;
	std::vector<std::string> format_list;
	std::vector<U32> info_ID_list;
	std::vector<U32> format_ID_list;

	//
	std::vector<map_type> load_info_ID_loaded;
	std::vector<map_type> load_format_ID_loaded;
};

}
}

#endif /* CORE_BLOCKENTRYSETTINGS_H_ */
