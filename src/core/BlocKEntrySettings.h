#ifndef CORE_BLOCKENTRYSETTINGS_H_
#define CORE_BLOCKENTRYSETTINGS_H_

namespace Tachyon{
namespace Core{

struct BlockEntrySettings{

	bool loadPPA;
	bool loadGT;
	bool loadGTSimple;
	bool loadMetaHot;
	bool loadMetaCold;
	bool loadInfoAll;
	bool loadFormatAll;
	bool constructOccTable;

	std::vector<std::string> load_samples;
	std::vector<std::string> load_info_names;
	std::vector<std::string> load_format_names;
};

}
}

#endif /* CORE_BLOCKENTRYSETTINGS_H_ */
