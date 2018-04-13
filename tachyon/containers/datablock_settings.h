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

/**<
 * Settings
 */
struct DataBlockSettings{
	typedef DataBlockSettings self_type;
	typedef SettingsMap       map_type;

	DataBlockSettings() :
		loadContig_(false),
		loadPositons_(false),
		loadController_(false),
		loadQuality_(false),
		loadNames_(false),
		loadAlleles_(false),
		loadSetMembership_(false),
		loadGenotypesAll_(false),
		loadGenotypesRLE_(false),
		loadGenotypesSimple_(false),
		loadGenotypesOther_(false),
		loadGenotypesSupport_(false),
		loadPPA_(false),
		loadINFO_(false),
		loadFORMAT_(false),
		constructOccTable_(false)
	{}

	self_type& loadAll(const bool set = true){
		loadContig_ = set;
		loadPositons_ = set;
		loadController_ = set;
		loadQuality_ = set;
		loadNames_ = set;
		loadAlleles_ = set;
		loadSetMembership_ = set;
		loadGenotypesAll_ = set;
		loadGenotypesRLE_ = set;
		loadGenotypesSimple_ = set;
		loadGenotypesOther_ = set;
		loadGenotypesSupport_ = set;
		loadINFO_ = set;
		loadFORMAT_ = set;
		loadPPA_ = set;
		return(*this);
	}

	self_type& loadAllMeta(const bool set = true){
		loadContig_ = set;
		loadPositons_ = set;
		loadController_ = set;
		loadQuality_ = set;
		loadNames_ = set;
		loadAlleles_ = set;
		return(*this);
	}

	self_type& loadAllFILTER(const bool set = true){
		loadSetMembership_ = set;
		return(*this);
	}

	self_type& loadAllINFO(const bool set = true){
		//loadContig_ = set;
		//loadPositons_ = set;
		//loadSetMembership_ = set;
		loadINFO_ = set;
		return(*this);
	}

	self_type& loadINFO(const std::string& field_name){
		if(field_name.size() == 0) return(*this);
		this->loadContig_ = true;
		this->loadPositons_ = true;
		this->loadSetMembership_ = true;
		this->info_list.push_back(field_name);
		return(*this);
	}

	self_type& loadINFO(const U32 field_id){
		this->info_ID_list.push_back(field_id);
		this->loadContig_ = true;
		this->loadPositons_ = true;
		this->loadSetMembership_ = true;
		return(*this);
	}

	self_type& loadGenotypes(const bool set){
		//loadContig_ = set;
		//loadPositons_ = set;
		//loadController_ = set;
		//loadSetMembership_ = set;
		loadGenotypesAll_ = set;
		loadGenotypesRLE_ = set;
		loadGenotypesSimple_ = set;
		loadGenotypesOther_ = set;
		loadGenotypesSupport_ = set;
		return(*this);
	}

	self_type& loadFORMAT(const bool set){
		this->loadPPA_ = set;
		this->loadGenotypes(set);
		this->loadFORMAT_ = set;
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


public:
	bool loadContig_;
	bool loadPositons_;
	bool loadController_;
	bool loadQuality_;
	bool loadNames_;
	bool loadAlleles_;
	bool loadSetMembership_;
	bool loadGenotypesAll_;
	bool loadGenotypesRLE_;
	bool loadGenotypesSimple_;
	bool loadGenotypesOther_;
	bool loadGenotypesSupport_;
	bool loadPPA_;
	bool loadINFO_;
	bool loadFORMAT_;
	bool constructOccTable_;

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
