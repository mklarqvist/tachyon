#ifndef CORE_BLOCKENTRYSETTINGS_H_
#define CORE_BLOCKENTRYSETTINGS_H_

namespace Tachyon{
namespace Core{

/**<
 * Supportive structure for BlockEntry
 */
struct BlockEntrySettingsMap{
	typedef BlockEntrySettingsMap self_type;
	typedef Index::BlockIndexHeaderOffsets offset_minimal_type;

	BlockEntrySettingsMap() : iterator_index(0), target_stream_local(-1), offset(nullptr){}
	BlockEntrySettingsMap(const U32 iterator_index, const S32 target_stream_disk, const offset_minimal_type* offset) : iterator_index(iterator_index), target_stream_local(target_stream_disk), offset(offset){}
	~BlockEntrySettingsMap(){}

	BlockEntrySettingsMap(const BlockEntrySettingsMap& other) : iterator_index(other.iterator_index), target_stream_local(other.target_stream_local), offset(other.offset){}
	BlockEntrySettingsMap(BlockEntrySettingsMap&& other) : iterator_index(other.iterator_index), target_stream_local(other.target_stream_local), offset(other.offset){}
	BlockEntrySettingsMap& operator=(const BlockEntrySettingsMap& other){
		this->iterator_index = other.iterator_index;
		this->target_stream_local = other.target_stream_local;
		this->offset = other.offset;
		return *this;
	}
	BlockEntrySettingsMap& operator=(BlockEntrySettingsMap&& other){
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
	typedef BlockEntrySettingsMap map_type;

	BlockEntrySettings() :
		loadPPA(false),
		loadGT(false),
		loadGTSimple(false),
		loadMetaHot(false),
		loadMetaCold(false),
		loadInfoAll(false),
		loadFormatAll(false),
		constructOccTable(false)
	{}

	void loadAll(void){
		this->loadPPA = true;
		this->loadGT = true;
		this->loadGTSimple = true;
		this->loadMetaCold = true;
		this->loadMetaHot = true;
		this->loadInfoAll = true;
		this->loadFormatAll = true;
	}

	void unsetFormat(void){
		this->loadPPA = false;
		this->loadGT = false;
		this->loadGTSimple = false;
		this->loadFormatAll = false;
		this->constructOccTable = false;
	}

	void unsetInfo(void){
		this->loadInfoAll = false;
	}

	bool loadPPA;
	bool loadGT;
	bool loadGTSimple;
	bool loadMetaHot;
	bool loadMetaCold;
	bool loadInfoAll;
	bool loadFormatAll;
	bool constructOccTable;

	std::vector<std::string> load_samples;
	std::vector<std::string> load_info_names;   // has to be in order
	std::vector<std::string> load_format_names; // has to be in order
	std::vector<U32> load_info_ID;
	std::vector<U32> load_format_ID;


	std::vector<map_type> load_info_ID_loaded;
	std::vector<map_type> load_format_ID_loaded;
};

}
}

#endif /* CORE_BLOCKENTRYSETTINGS_H_ */
