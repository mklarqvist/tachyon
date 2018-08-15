#ifndef CONTAINERS_VARIANT_READER_OBJECTS_H_
#define CONTAINERS_VARIANT_READER_OBJECTS_H_

#include <unordered_map>

#include "containers/meta_container.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/info_container_string.h"
#include "containers/format_container.h"
#include "containers/format_container_string.h"
#include "occ.h"

namespace tachyon{

/**<
 * The function of this struct is to keep all loaded containers and
 * field identifiers in a single object. Loading the members of this object
 * occurs OUTSIDE this definition.
 */
struct VariantReaderObjects{
public:
	typedef VariantReaderObjects                 self_type;
	typedef containers::MetaContainer            meta_container_type;
	typedef containers::GenotypeContainer        gt_container_type;
	typedef containers::InfoContainerInterface   info_interface_type;
	typedef containers::FormatContainerInterface format_interface_type;
	typedef yon_gt_summary                       genotype_summary_type;

public:
	VariantReaderObjects() :
		loaded_genotypes(false),
		loaded_meta(false),
		n_loaded_info(0),
		n_loaded_format(0),
		meta_container(nullptr),
		genotype_container(nullptr),
		genotype_summary(nullptr),
		info_containers(nullptr),
		format_containers(nullptr),
		occ(nullptr)
	{
	}

	~VariantReaderObjects(){
		delete this->meta_container;
		delete this->genotype_container;
		delete this->genotype_summary;

		for(U32 i = 0; i < this->info_id_loaded.size(); ++i)
			delete this->info_containers[this->info_id_loaded[i]];
		delete [] this->info_containers;

		for(U32 i = 0; i < this->format_id_loaded.size(); ++i)
			delete this->format_containers[this->format_id_loaded[i]];
		delete [] this->format_containers;
		delete this->occ;
	}

	bool EvaluateOcc(){
		if(occ == nullptr) return false;
		return(this->occ->BuildTable());
	}

	bool EvaluateOcc(yon_gt_ppa* ppa){
		if(occ == nullptr) return(this->EvaluateOcc());
		return(this->occ->BuildTable(*ppa));
	}

public:
	bool loaded_genotypes;
	bool loaded_meta;
	size_t n_loaded_info;
	size_t n_loaded_format;
	std::vector<int> info_id_loaded;
	std::vector<int> format_id_loaded;

	std::unordered_map<std::string, info_interface_type*>   info_container_map;
	std::unordered_map<std::string, format_interface_type*> format_container_map;
	meta_container_type*     meta_container;
	gt_container_type*       genotype_container;
	genotype_summary_type*   genotype_summary;
	info_interface_type**    info_containers;
	format_interface_type**  format_containers;
	yon_occ* occ;
};

}



#endif /* CONTAINERS_VARIANT_READER_OBJECTS_H_ */
