#ifndef CONTAINERS_VARIANT_READER_OBJECTS_H_
#define CONTAINERS_VARIANT_READER_OBJECTS_H_

#include <unordered_map>

#include "containers/meta_container.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/info_container_string.h"
#include "containers/format_container.h"
#include "containers/format_container_string.h"

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
	typedef containers::GenotypeSummary          genotype_summary_type;

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
		format_containers(nullptr)
	{}

	~VariantReaderObjects(){
		delete this->meta_container;
		delete this->genotype_container;
		delete this->genotype_summary;

		for(U32 i = 0; i < this->n_loaded_info; ++i) delete this->info_containers[i];
		delete [] this->info_containers;

		for(U32 i = 0; i < this->n_loaded_format; ++i) delete this->format_containers[i];
		delete [] this->format_containers;
	}

public:
	bool loaded_genotypes;
	bool loaded_meta;
	size_t n_loaded_info;
	size_t n_loaded_format;
	std::vector<U32> info_id_fields_keep;
	std::vector<U32> format_id_fields_keep;
	std::vector<U16> additional_info_execute_flag_set;
	std::vector< std::vector<U32> > local_match_keychain_info;
	std::vector< std::vector<U32> > local_match_keychain_format;

	std::vector<std::string> info_field_names;
	std::vector<std::string> format_field_names;

	std::unordered_map<std::string, info_interface_type*>   info_container_map;
	std::unordered_map<std::string, format_interface_type*> format_container_map;

	meta_container_type*     meta_container;
	gt_container_type*       genotype_container;
	genotype_summary_type*   genotype_summary;
	info_interface_type**    info_containers;
	format_interface_type**  format_containers;
};

}



#endif /* CONTAINERS_VARIANT_READER_OBJECTS_H_ */
