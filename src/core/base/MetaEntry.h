#ifndef CORE_BASE_METAENTRY_H_
#define CORE_BASE_METAENTRY_H_

namespace Tachyon{
namespace Core{

struct MetaEntry{
	typedef MetaHot* hot_entry;
	typedef MetaCold* cold_entry;


	hot_entry* hot;
	cold_entry* cold;
};

}
}



#endif /* CORE_BASE_METAENTRY_H_ */
