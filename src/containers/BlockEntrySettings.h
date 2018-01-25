#ifndef CORE_BLOCKENTRYSETTINGS_H_
#define CORE_BLOCKENTRYSETTINGS_H_

namespace tachyon{
namespace core{

/**<
 * Supportive structure for Block
 */
struct SettingsMap{
	typedef SettingsMap self_type;
	typedef index::BlockIndexOffsetsHeader offset_minimal_type;

	SettingsMap() : iterator_index(0), target_stream_local(-1), offset(nullptr){}
	SettingsMap(const U32 iterator_index, const S32 target_stream_disk, const offset_minimal_type* offset) : iterator_index(iterator_index), target_stream_local(target_stream_disk), offset(offset){}
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

	inline bool operator<(const self_type& other) const{ return(this->offset->offset < other.offset->offset); }
	inline bool operator>(const self_type& other) const{ return(!((*this) < other)); }

	U32 iterator_index;
	S32 target_stream_local;
	const offset_minimal_type* offset;
};

/**<
 * Settings
 */
struct BlockEntrySettings{
	typedef BlockEntrySettings    self_type;
	typedef SettingsMap map_type;

	BlockEntrySettings() :
		importPPA(false),
		importGT(false),
		importGTSimple(false),
		importMetaHot(false),
		importMetaCold(false),
		importInfoAll(false),
		importFormatAll(false),
		constructOccTable(false)
	{}

	self_type& loadAll(const bool set = true){
		this->importPPA = set;
		this->importGT = set;
		this->importGTSimple = set;
		this->importMetaCold = set;
		this->importMetaHot = set;
		this->importInfoAll = set;
		this->importFormatAll = set;
		return(*this);
	}

	self_type& loadMeta(const bool set = true){
		this->importMetaHot = set;
		this->importMetaCold = set;
		return(*this);
	}

	self_type& loadMetaHot(const bool set = true){
		this->importMetaHot = true;
		return(*this);
	}

	self_type& loadINFO(const bool set = true){
		this->importInfoAll = set;
		return(*this);
	}

	self_type& loadINFO(const std::string& field_name){
		if(field_name.size() == 0) return(*this);
		this->info_list.push_back(field_name);
		return(*this);
	}

	self_type& loadINFO(const U32 field_id){
		this->info_ID_list.push_back(field_id);
		return(*this);
	}

	self_type& loadGenotypes(const bool set = true){
		this->loadMeta(true);
		this->importGT = true;
		this->importGTSimple = true;
		return(*this);
	}

	self_type& loadFORMAT(const bool set = true){
		this->importPPA = set;
		this->importGT = set;
		this->importGTSimple = set;
		this->importFormatAll = set;
		return(*this);
	}

	self_type& loadFORMAT(const std::string& field_name){
		if(field_name.size() == 0) return(*this);
		this->format_list.push_back(field_name);
		return(*this);
	}

	self_type& loadFORMAT(const U32 field_id){
		this->format_ID_list.push_back(field_id);
		return(*this);
	}


public:
	bool importPPA;
	bool importGT;
	bool importGTSimple;
	bool importMetaHot;
	bool importMetaCold;
	bool importInfoAll;
	bool importFormatAll;
	bool constructOccTable;

	std::vector<std::string> samples_list;
	std::vector<std::string> info_list;   // has to be in order
	std::vector<std::string> format_list; // has to be in order
	std::vector<U32> info_ID_list;
	std::vector<U32> format_ID_list;

	//
	std::vector<map_type> load_info_ID_loaded;
	std::vector<map_type> load_format_ID_loaded;
};

}
}

#endif /* CORE_BLOCKENTRYSETTINGS_H_ */
