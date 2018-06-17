#ifndef CONTAINERS_VARIANT_READER_OBJECTS_H_
#define CONTAINERS_VARIANT_READER_OBJECTS_H_

#include "meta_container.h"
#include "genotype_container.h"
#include "info_container.h"
#include "info_container_string.h"
#include "format_container.h"
#include "format_container_string.h"

namespace tachyon{

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
		meta(nullptr),
		genotypes(nullptr),
		genotype_summary(nullptr),
		info_fields(nullptr),
		format_fields(nullptr)
	{}

	~VariantReaderObjects(){
		delete this->meta;
		delete this->genotypes;
		delete this->genotype_summary;

		for(U32 i = 0; i < this->n_loaded_info; ++i) delete this->info_fields[i];
		delete [] this->info_fields;

		for(U32 i = 0; i < this->n_loaded_format; ++i) delete this->format_fields[i];
		delete [] this->format_fields;
	}

public:
	bool loaded_genotypes;
	bool loaded_meta;
	size_t n_loaded_info;
	size_t n_loaded_format;
	std::vector<U32> info_keep;
	std::vector<U32> format_keep;
	std::vector< U16 > additional_info_execute_flag_set;
	std::vector< std::vector<U32> > local_match_keychain_info;
	std::vector< std::vector<U32> > local_match_keychain_format;

	std::vector<std::string> info_field_names;
	std::vector<std::string> format_field_names;

	meta_container_type*     meta;
	gt_container_type*       genotypes;
	genotype_summary_type*   genotype_summary;
	info_interface_type**    info_fields;
	format_interface_type**  format_fields;
};

}



#endif /* CONTAINERS_VARIANT_READER_OBJECTS_H_ */
