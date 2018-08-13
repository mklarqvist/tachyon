#ifndef CORE_OCC_H_
#define CORE_OCC_H_

#include <cstdint>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <unordered_map>

#include "support/helpers.h"
#include "header/variant_header.h"
#include "genotypes.h"

namespace tachyon{

struct yon_occ {
	typedef std::unordered_map<std::string, uint32_t> map_type;

	yon_occ() = default;
	~yon_occ() = default;

	bool ReadTable(const std::string file_name, const VariantHeader& header, const char delimiter = '\t');
	bool BuildTable(void);
	bool BuildTable(const yon_gt_ppa& ppa);

	// Map from group name to row offset in the table.
	map_type map;
	// Unique names of grouping factors.
	std::vector<std::string> row_names;

	// a matrix with proportions samples times groupings
	// rows corresponds to the cumulative sum of a grouping
	// over the samples. The table corresponds to the set
	// membership (presence or absence) and the occ table
	// corresponds to the cumsum.
	std::vector< std::vector<uint32_t> > table;
	std::vector< std::vector<uint32_t> > occ;
};

}

#endif /* CORE_OCC_H_ */
