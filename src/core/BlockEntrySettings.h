#ifndef CORE_BLOCKENTRYSETTINGS_H_
#define CORE_BLOCKENTRYSETTINGS_H_

namespace Tachyon{
namespace Core{

struct BlockEntrySettingsMap{
	typedef BlockEntrySettingsMap self_type;

	BlockEntrySettingsMap() : key(0), target_stream(-1), target_stream_local(-1), offset(0){}
	BlockEntrySettingsMap(const U32& key, const S32 target_stream, const S32 target_stream_disk, const U32& offset) : key(key), target_stream(target_stream), target_stream_local(target_stream_disk), offset(offset){}

	inline bool operator<(const self_type& other) const{ return(this->offset < other.offset); }
	inline bool operator>(const self_type& other) const{ return(!((*this) < other)); }

	U32 key;
	S32 target_stream;
	S32 target_stream_local;
	U32 offset;
};

struct BlockEntrySettings{
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


	std::vector<BlockEntrySettingsMap> load_info_ID_loaded;
	std::vector<BlockEntrySettingsMap> load_format_ID_loaded;
};

}
}

#endif /* CORE_BLOCKENTRYSETTINGS_H_ */
