#include "variant_block_footer.h"

#include <cassert>


namespace tachyon{
namespace containers{

VariantBlockFooter::VariantBlockFooter():
	l_info_bitvector(0),
	l_format_bitvector(0),
	l_filter_bitvector(0),
	n_info_streams(0),
	n_format_streams(0),
	n_filter_streams(0),
	n_info_patterns(0),
	n_format_patterns(0),
	n_filter_patterns(0),
	offsets(new header_type[YON_BLK_N_STATIC]),
	info_offsets(nullptr),
	format_offsets(nullptr),
	filter_offsets(nullptr),
	n_info_patterns_allocated(0),
	n_format_patterns_allocated(0),
	n_filter_patterns_allocated(0),
	info_patterns(nullptr),
	format_patterns(nullptr),
	filter_patterns(nullptr),
	info_map(nullptr),
	format_map(nullptr),
	filter_map(nullptr),
	info_pattern_map(nullptr),
	format_pattern_map(nullptr),
	filter_pattern_map(nullptr)
{}

VariantBlockFooter::~VariantBlockFooter(){
	delete [] this->offsets;
	delete [] this->info_offsets;
	delete [] this->format_offsets;
	delete [] this->filter_offsets;
	delete [] this->info_patterns;
	delete [] this->format_patterns;
	delete [] this->filter_patterns;
	delete this->info_map;
	delete this->format_map;
	delete this->filter_map;
	delete this->info_pattern_map;
	delete this->format_pattern_map;
	delete this->filter_pattern_map;
}

void VariantBlockFooter::reset(void){
	// Headers of the various containers
	for(U32 i = 0; i < YON_BLK_N_STATIC; ++i) this->offsets[i].reset();

	for(U32 i = 0; i < this->n_info_streams; ++i)   this->info_offsets[i].reset();
	for(U32 i = 0; i < this->n_format_streams; ++i) this->format_offsets[i].reset();
	for(U32 i = 0; i < this->n_filter_streams; ++i) this->filter_offsets[i].reset();

	for(U32 i = 0; i < this->n_info_patterns; ++i)   this->info_patterns[i].clear();
	for(U32 i = 0; i < this->n_format_patterns; ++i) this->format_patterns[i].clear();
	for(U32 i = 0; i < this->n_filter_patterns; ++i) this->filter_patterns[i].clear();

	this->n_info_streams    = 0;
	this->n_format_streams  = 0;
	this->n_filter_streams  = 0;
	this->n_info_patterns   = 0;
	this->n_format_patterns = 0;
	this->n_filter_patterns = 0;

	this->resetTables();
}

void VariantBlockFooter::resetTables(){
	if(this->info_map != nullptr)   this->info_map->clear();
	if(this->format_map != nullptr) this->format_map->clear();
	if(this->filter_map != nullptr) this->filter_map->clear();
	if(this->info_pattern_map != nullptr)   this->info_pattern_map->clear();
	if(this->format_pattern_map != nullptr) this->format_pattern_map->clear();
	if(this->filter_pattern_map != nullptr) this->filter_pattern_map->clear();
}

io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VariantBlockFooter& entry){
	buffer += (U16)entry.n_info_streams;
	buffer += (U16)entry.n_format_streams;
	buffer += (U16)entry.n_filter_streams;
	buffer += (U16)entry.n_info_patterns;
	buffer += (U16)entry.n_format_patterns;
	buffer += (U16)entry.n_filter_patterns;

	for(U32 i = 0; i < YON_BLK_N_STATIC; ++i)       buffer << entry.offsets[i];
	for(U32 i = 0; i < entry.n_info_streams; ++i)   buffer << entry.info_offsets[i];
	for(U32 i = 0; i < entry.n_format_streams; ++i) buffer << entry.format_offsets[i];
	for(U32 i = 0; i < entry.n_filter_streams; ++i) buffer << entry.filter_offsets[i];

	for(U32 i = 0; i < entry.n_info_patterns; ++i)   buffer << entry.info_patterns[i];
	for(U32 i = 0; i < entry.n_format_patterns; ++i) buffer << entry.format_patterns[i];
	for(U32 i = 0; i < entry.n_filter_patterns; ++i) buffer << entry.filter_patterns[i];

	return(buffer);
}


io::BasicBuffer& operator>>(io::BasicBuffer& buffer, VariantBlockFooter& entry){
	entry.reset();

	buffer >> entry.n_info_streams;
	buffer >> entry.n_format_streams;
	buffer >> entry.n_filter_streams;
	buffer >> entry.n_info_patterns;
	buffer >> entry.n_format_patterns;
	buffer >> entry.n_filter_patterns;

	entry.l_info_bitvector   = ceil((float)entry.n_info_streams    / 8);
	entry.l_format_bitvector = ceil((float)entry.n_format_streams  / 8);
	entry.l_filter_bitvector = ceil((float)entry.n_filter_streams  / 8);

	entry.BuildMaps(); // Construct new maps.
	entry.BuildPatternMaps(); // Construct new pattern maps.
	entry.offsets        = new DataContainerHeader[YON_BLK_N_STATIC];
	entry.info_offsets   = new DataContainerHeader[entry.n_info_streams];
	entry.format_offsets = new DataContainerHeader[entry.n_format_streams];
	entry.filter_offsets = new DataContainerHeader[entry.n_filter_streams];
	entry.n_info_patterns_allocated   = entry.n_info_streams;
	entry.n_format_patterns_allocated = entry.n_format_streams;
	entry.n_filter_patterns_allocated = entry.n_filter_streams;

	for(U32 i = 0; i < YON_BLK_N_STATIC; ++i)
		buffer >> entry.offsets[i];

	for(U32 i = 0; i < entry.n_info_streams; ++i){
		buffer >> entry.info_offsets[i];
		entry.UpdateInfo(entry.info_offsets[i], i);
	}

	for(U32 i = 0; i < entry.n_format_streams; ++i){
		buffer >> entry.format_offsets[i];
		entry.UpdateFormat(entry.format_offsets[i], i);
	}

	for(U32 i = 0; i < entry.n_filter_streams; ++i){
		buffer >> entry.filter_offsets[i];
		entry.UpdateFilter(entry.filter_offsets[i], i);
	}

	entry.info_patterns = new yon_blk_bv_pair[entry.n_info_patterns];
	for(U32 i = 0; i < entry.n_info_patterns; ++i){
		buffer >> entry.info_patterns[i];
		entry.UpdateInfoPattern(entry.info_patterns[i].pattern, i);
		entry.info_patterns[i].Build(entry.n_info_streams, entry.info_map);
	}

	entry.format_patterns = new yon_blk_bv_pair[entry.n_format_patterns];
	for(U32 i = 0; i < entry.n_format_patterns; ++i){
		buffer >> entry.format_patterns[i];
		entry.UpdateFormatPattern(entry.format_patterns[i].pattern, i);
		entry.format_patterns[i].Build(entry.n_format_streams, entry.format_map);
	}

	entry.filter_patterns = new yon_blk_bv_pair[entry.n_filter_patterns];
	for(U32 i = 0; i < entry.n_filter_patterns; ++i){
		buffer >> entry.filter_patterns[i];
		entry.UpdateFilterPattern(entry.filter_patterns[i].pattern, i);
		entry.filter_patterns[i].Build(entry.n_filter_streams, entry.filter_map);
	}

	return(buffer);
}

}
}
