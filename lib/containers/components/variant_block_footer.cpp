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

VariantBlockFooter::VariantBlockFooter(const self_type& other) :
	l_info_bitvector(other.l_info_bitvector),
	l_format_bitvector(other.l_format_bitvector),
	l_filter_bitvector(other.l_filter_bitvector),
	n_info_streams(other.n_info_streams),
	n_format_streams(other.n_format_streams),
	n_filter_streams(other.n_filter_streams),
	n_info_patterns(other.n_info_patterns),
	n_format_patterns(other.n_format_patterns),
	n_filter_patterns(other.n_filter_patterns),
	offsets(new header_type[YON_BLK_N_STATIC]),
	info_offsets(new header_type[other.n_info_streams * 2]),
	format_offsets(new header_type[other.n_format_streams * 2]),
	filter_offsets(new header_type[other.n_filter_streams * 2]),
	n_info_patterns_allocated(other.n_info_patterns_allocated),
	n_format_patterns_allocated(other.n_format_patterns_allocated),
	n_filter_patterns_allocated(other.n_filter_patterns_allocated),
	info_patterns(new yon_blk_bv_pair[other.n_info_patterns_allocated]),
	format_patterns(new yon_blk_bv_pair[other.n_format_patterns_allocated]),
	filter_patterns(new yon_blk_bv_pair[other.n_filter_patterns_allocated]),
	info_map(nullptr),
	format_map(nullptr),
	filter_map(nullptr),
	info_pattern_map(nullptr),
	format_pattern_map(nullptr),
	filter_pattern_map(nullptr)
{
	for(uint32_t i = 0; i < YON_BLK_N_STATIC; ++i) this->offsets[i] = other.offsets[i];

	for(uint32_t i = 0; i < this->n_info_streams; ++i)   this->info_offsets[i]   = other.info_offsets[i];
	for(uint32_t i = 0; i < this->n_format_streams; ++i) this->format_offsets[i] = other.format_offsets[i];
	for(uint32_t i = 0; i < this->n_filter_streams; ++i) this->filter_offsets[i] = other.filter_offsets[i];

	for(uint32_t i = 0; i < this->n_info_patterns; ++i)   this->info_patterns[i]   = other.info_patterns[i];
	for(uint32_t i = 0; i < this->n_format_patterns; ++i) this->format_patterns[i] = other.format_patterns[i];
	for(uint32_t i = 0; i < this->n_filter_patterns; ++i) this->filter_patterns[i] = other.filter_patterns[i];

	if(other.info_map != nullptr)   this->info_map   = new map_type(*other.info_map);
	if(other.format_map != nullptr) this->format_map = new map_type(*other.format_map);
	if(other.filter_map != nullptr) this->filter_map = new map_type(*other.filter_map);

	if(other.info_pattern_map != nullptr)   this->info_pattern_map   = new map_pattern_type(*other.info_pattern_map);
	if(other.format_pattern_map != nullptr) this->format_pattern_map = new map_pattern_type(*other.format_pattern_map);
	if(other.filter_pattern_map != nullptr) this->filter_pattern_map = new map_pattern_type(*other.filter_pattern_map);
}

VariantBlockFooter::VariantBlockFooter(self_type&& other) noexcept :
	l_info_bitvector(other.l_info_bitvector),
	l_format_bitvector(other.l_format_bitvector),
	l_filter_bitvector(other.l_filter_bitvector),
	n_info_streams(other.n_info_streams),
	n_format_streams(other.n_format_streams),
	n_filter_streams(other.n_filter_streams),
	n_info_patterns(other.n_info_patterns),
	n_format_patterns(other.n_format_patterns),
	n_filter_patterns(other.n_filter_patterns),
	offsets(nullptr),
	info_offsets(nullptr),
	format_offsets(nullptr),
	filter_offsets(nullptr),
	n_info_patterns_allocated(other.n_info_patterns_allocated),
	n_format_patterns_allocated(other.n_format_patterns_allocated),
	n_filter_patterns_allocated(other.n_filter_patterns_allocated),
	info_patterns(nullptr),
	format_patterns(nullptr),
	filter_patterns(nullptr),
	info_map(nullptr),
	format_map(nullptr),
	filter_map(nullptr),
	info_pattern_map(nullptr),
	format_pattern_map(nullptr),
	filter_pattern_map(nullptr)
{
	std::swap(this->offsets, other.offsets);
	std::swap(this->info_offsets, other.info_offsets);
	std::swap(this->format_offsets, other.format_offsets);
	std::swap(this->filter_offsets, other.filter_offsets);
	std::swap(this->info_patterns, other.info_patterns);
	std::swap(this->format_patterns, other.format_patterns);
	std::swap(this->filter_patterns, other.filter_patterns);
	std::swap(this->info_map, other.info_map);
	std::swap(this->format_map, other.format_map);
	std::swap(this->filter_map, other.filter_map);
	std::swap(this->info_pattern_map, other.info_pattern_map);
	std::swap(this->format_pattern_map, other.format_pattern_map);
	std::swap(this->filter_pattern_map, other.filter_pattern_map);
}

VariantBlockFooter& VariantBlockFooter::operator=(const self_type& other){
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
	*this = VariantBlockFooter(other);
	return(*this);
}

VariantBlockFooter& VariantBlockFooter::operator=(self_type&& other) noexcept{
	if(this == &other){
		// precautions against self-moves
		return *this;
	}

	delete [] this->offsets; this->offsets = nullptr;
	delete [] this->info_offsets; this->info_offsets = nullptr;
	delete [] this->format_offsets; this->format_offsets = nullptr;
	delete [] this->filter_offsets; this->filter_offsets = nullptr;
	delete [] this->info_patterns; this->info_patterns = nullptr;
	delete [] this->format_patterns; this->format_patterns = nullptr;
	delete [] this->filter_patterns; this->filter_patterns = nullptr;
	delete this->info_map; this->info_map = nullptr;
	delete this->format_map; this->format_map = nullptr;
	delete this->filter_map; this->filter_map = nullptr;
	delete this->info_pattern_map; this->info_pattern_map = nullptr;
	delete this->format_pattern_map; this->format_pattern_map = nullptr;
	delete this->filter_pattern_map; this->filter_pattern_map = nullptr;

	std::swap(this->offsets, other.offsets);
	std::swap(this->info_offsets, other.info_offsets);
	std::swap(this->format_offsets, other.format_offsets);
	std::swap(this->filter_offsets, other.filter_offsets);
	std::swap(this->info_patterns, other.info_patterns);
	std::swap(this->format_patterns, other.format_patterns);
	std::swap(this->filter_patterns, other.filter_patterns);
	std::swap(this->info_map, other.info_map);
	std::swap(this->format_map, other.format_map);
	std::swap(this->filter_map, other.filter_map);
	std::swap(this->info_pattern_map, other.info_pattern_map);
	std::swap(this->format_pattern_map, other.format_pattern_map);
	std::swap(this->filter_pattern_map, other.filter_pattern_map);

	this->l_info_bitvector   = other.l_info_bitvector;
	this->l_format_bitvector = other.l_format_bitvector;
	this->l_filter_bitvector = other.l_filter_bitvector;
	this->n_info_streams     = other.n_info_streams;
	this->n_format_streams   = other.n_format_streams;
	this->n_filter_streams   = other.n_filter_streams;
	this->n_info_patterns    = other.n_info_patterns;
	this->n_format_patterns  = other.n_format_patterns;
	this->n_filter_patterns  = other.n_filter_patterns;
	this->n_info_patterns_allocated   = other.n_info_patterns_allocated;
	this->n_format_patterns_allocated = other.n_format_patterns_allocated;
	this->n_filter_patterns_allocated = other.n_filter_patterns_allocated;

	return(*this);
}

void VariantBlockFooter::reset(void){
	// Headers of the various containers
	for(uint32_t i = 0; i < YON_BLK_N_STATIC; ++i) this->offsets[i].reset();

	for(uint32_t i = 0; i < this->n_info_streams; ++i)   this->info_offsets[i].reset();
	for(uint32_t i = 0; i < this->n_format_streams; ++i) this->format_offsets[i].reset();
	for(uint32_t i = 0; i < this->n_filter_streams; ++i) this->filter_offsets[i].reset();

	for(uint32_t i = 0; i < this->n_info_patterns; ++i)   this->info_patterns[i].clear();
	for(uint32_t i = 0; i < this->n_format_patterns; ++i) this->format_patterns[i].clear();
	for(uint32_t i = 0; i < this->n_filter_patterns; ++i) this->filter_patterns[i].clear();

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

// Allocate offset vectors
void VariantBlockFooter::AllocateInfoHeaders(const uint32_t n_info_streams){
	delete [] this->info_offsets;
	if(n_info_streams == 0){
		this->info_offsets = nullptr;
		return;
	}
	this->info_offsets = new header_type[n_info_streams];
}

void VariantBlockFooter::AllocateFormatHeaders(const uint32_t n_format_streams){
	delete [] this->format_offsets;
	if(n_format_streams == 0){
		this->format_offsets = nullptr;
		return;
	}
	this->format_offsets = new header_type[n_format_streams];
}

void VariantBlockFooter::AllocateFilterHeaders(const uint32_t n_filter_streams){
	delete [] this->filter_offsets;
	if(n_filter_streams == 0){
		this->filter_offsets = nullptr;
		return;
	}
	this->filter_offsets = new header_type[n_filter_streams];
}

void VariantBlockFooter::AllocateHeaders(const uint32_t n_info_streams,
							const uint32_t n_format_streams,
							const uint32_t n_filter_streams)
{
	this->AllocateInfoHeaders(n_info_streams);
	this->AllocateFormatHeaders(n_format_streams);
	this->AllocateFilterHeaders(n_filter_streams);
}

bool VariantBlockFooter::ConstructInfoBitVector(std::unordered_map<uint32_t,uint32_t>* pattern_map){
	for(uint32_t i = 0; i < this->n_info_patterns; ++i){
		this->info_patterns[i].Build(this->n_info_streams, pattern_map);
	}
	return true;
}

bool VariantBlockFooter::ConstructFormatBitVector(std::unordered_map<uint32_t,uint32_t>* pattern_map){
	for(uint32_t i = 0; i < this->n_format_patterns; ++i){
		this->format_patterns[i].Build(this->n_format_streams, pattern_map);
	}
	return true;
}

bool VariantBlockFooter::ConstructFilterBitVector(std::unordered_map<uint32_t,uint32_t>* pattern_map){
	for(uint32_t i = 0; i < this->n_filter_patterns; ++i){
		this->filter_patterns[i].Build(this->n_filter_streams, pattern_map);
	}
	return true;
}

uint32_t VariantBlockFooter::AddPatternWrapper(const std::vector<int>& pattern,
					  map_pattern_type* pattern_map,
					  yon_blk_bv_pair* bv_pairs,
					  uint16_t& stream_counter)
{
	uint64_t pattern_hash = VariantBlockFooter::HashIdentifiers(pattern);
	const map_pattern_type::const_iterator it = pattern_map->find(pattern_hash); // search for pattern
	if(it == pattern_map->end()){
		(*pattern_map)[pattern_hash] = stream_counter;
		bv_pairs[stream_counter].pattern = pattern;
		++stream_counter;
	}

	return((*pattern_map)[pattern_hash]);
}

uint32_t VariantBlockFooter::AddInfoPattern(const std::vector<int>& pattern){
	if(this->info_pattern_map == nullptr) this->BuildPatternMaps();
	if(this->n_info_patterns_allocated == 0){
		delete [] this->info_patterns;
		this->info_patterns = new yon_blk_bv_pair[100];
		this->n_info_patterns_allocated = 100;
	}

	// Resize if required.
	if(this->n_info_patterns == this->n_info_patterns_allocated){
		yon_blk_bv_pair* temp = this->info_patterns;

		this->info_patterns = new yon_blk_bv_pair[this->n_info_patterns_allocated*2];
		for(uint32_t i = 0; i < this->n_info_patterns_allocated; ++i){
			this->info_patterns[i] = std::move(temp[i]);
		}
		this->n_info_patterns_allocated *= 2;
		delete [] temp;
	}

	return(this->AddPatternWrapper(pattern,
								   this->info_pattern_map,
								   this->info_patterns,
								   this->n_info_patterns));
}

uint32_t VariantBlockFooter::AddFormatPattern(const std::vector<int>& pattern){
	if(this->format_pattern_map == nullptr) this->BuildPatternMaps();
	if(this->n_format_patterns_allocated == 0){
		delete [] this->format_patterns;
		this->format_patterns = new yon_blk_bv_pair[100];
		this->n_format_patterns_allocated = 100;
	}

	// Resize if required.
	if(this->n_format_patterns == this->n_format_patterns_allocated){
		yon_blk_bv_pair* temp = this->format_patterns;

		this->format_patterns = new yon_blk_bv_pair[this->n_format_patterns_allocated*2];
		for(uint32_t i = 0; i < this->n_format_patterns_allocated; ++i){
			this->format_patterns[i] = std::move(temp[i]);
		}
		this->n_format_patterns_allocated *= 2;
		delete [] temp;
	}

	return(this->AddPatternWrapper(pattern,
								   this->format_pattern_map,
								   this->format_patterns,
								   this->n_format_patterns));
}

uint32_t VariantBlockFooter::AddFilterPattern(const std::vector<int>& pattern){
	if(this->filter_pattern_map == nullptr) this->BuildPatternMaps();
	if(this->n_filter_patterns_allocated == 0){
		delete [] this->filter_patterns;
		this->filter_patterns = new yon_blk_bv_pair[100];
		this->n_filter_patterns_allocated = 100;
	}

	// Resize if required.
	if(this->n_filter_patterns == this->n_filter_patterns_allocated){
		yon_blk_bv_pair* temp = this->filter_patterns;

		this->filter_patterns = new yon_blk_bv_pair[this->n_filter_patterns_allocated*2];
		for(uint32_t i = 0; i < this->n_filter_patterns_allocated; ++i){
			this->filter_patterns[i] = std::move(temp[i]);
		}
		this->n_filter_patterns_allocated *= 2;
		delete [] temp;
	}

	return(this->AddPatternWrapper(pattern,
								   this->filter_pattern_map,
								   this->filter_patterns,
								   this->n_filter_patterns));
}

// This wrapper adds patterns to the hash map when the data has
// already been loaded. This occurs when loading an object from
// disk/buffer.
uint32_t VariantBlockFooter::UpdatePatternWrapper(const std::vector<int>& pattern,
					  map_pattern_type* pattern_map,
					  const uint16_t& stream_counter)
{
	uint64_t pattern_hash = VariantBlockFooter::HashIdentifiers(pattern);
	const map_pattern_type::const_iterator it = pattern_map->find(pattern_hash); // search for pattern
	if(it == pattern_map->end())
		(*pattern_map)[pattern_hash] = stream_counter;

	return((*pattern_map)[pattern_hash]);
}

uint32_t VariantBlockFooter::UpdateInfoPattern(const std::vector<int>& pattern, const uint16_t pattern_id){
	if(this->info_pattern_map == nullptr) this->BuildPatternMaps();
	if(this->n_info_patterns_allocated == 0){
		delete [] this->info_patterns;
		this->info_patterns = new yon_blk_bv_pair[100];
		this->n_info_patterns_allocated = 100;
	}
	return(this->UpdatePatternWrapper(pattern, this->info_pattern_map, pattern_id));
}

uint32_t VariantBlockFooter::UpdateFormatPattern(const std::vector<int>& pattern, const uint16_t pattern_id){
	if(this->format_pattern_map == nullptr) this->BuildPatternMaps();
	if(this->n_format_patterns_allocated == 0){
		delete [] this->format_patterns;
		this->format_patterns = new yon_blk_bv_pair[100];
		this->n_format_patterns_allocated = 100;
	}
	return(this->UpdatePatternWrapper(pattern, this->format_pattern_map, pattern_id));
}

uint32_t VariantBlockFooter::UpdateFilterPattern(const std::vector<int>& pattern, const uint16_t pattern_id){
	if(this->filter_pattern_map == nullptr) this->BuildPatternMaps();
	if(this->n_filter_patterns_allocated == 0){
		delete [] this->filter_patterns;
		this->filter_patterns = new yon_blk_bv_pair[100];
		this->n_filter_patterns_allocated = 100;
	}
	return(this->UpdatePatternWrapper(pattern, this->filter_pattern_map, pattern_id));
}

void VariantBlockFooter::Finalize(void){
	this->ConstructInfoBitVector(this->info_map);
	this->ConstructFormatBitVector(this->format_map);
	this->ConstructFilterBitVector(this->filter_map);
}

bool VariantBlockFooter::BuildMaps(void){
	delete this->info_map;
	delete this->filter_map;
	delete this->format_map;

	this->info_map   = new map_type();
	this->filter_map = new map_type();
	this->format_map = new map_type();

	return true;
}

bool VariantBlockFooter::BuildPatternMaps(void){
	delete this->info_pattern_map;
	this->info_pattern_map = new map_pattern_type();
	if(this->n_info_patterns_allocated == 0){
		this->info_patterns = new yon_blk_bv_pair[100];
		this->n_info_patterns_allocated = 100;
	}

	delete this->filter_pattern_map;
	this->filter_pattern_map = new map_pattern_type();
	if(this->n_filter_patterns_allocated == 0){
		this->filter_patterns = new yon_blk_bv_pair[100];
		this->n_filter_patterns_allocated = 100;
	}

	delete this->format_pattern_map;
	this->format_pattern_map = new map_pattern_type();
	if(this->n_format_patterns_allocated == 0){
		this->format_patterns = new yon_blk_bv_pair[100];
		this->n_format_patterns_allocated = 100;
	}

	return true;
}

uint32_t VariantBlockFooter::UpdateOffsetMapWrapper(const header_type& offset, map_type* map, const uint16_t& stream_counter){
	map_type::const_iterator it = map->find(offset.data_header.global_key);
	if(it == map->end())
		(*map)[offset.data_header.global_key] = stream_counter;

	return((*map)[offset.data_header.global_key]);
}

uint32_t VariantBlockFooter::UpdateInfo(const header_type& offset, const uint16_t position){
	if(this->info_map == nullptr) this->BuildMaps();
	return(this->UpdateOffsetMapWrapper(offset, this->info_map, position));
}

uint32_t VariantBlockFooter::UpdateFormat(const header_type& offset, const uint16_t position){
	if(this->format_map == nullptr) this->BuildMaps();
	return(this->UpdateOffsetMapWrapper(offset, this->format_map, position));
}

uint32_t VariantBlockFooter::UpdateFilter(const header_type& offset, const uint16_t position){
	if(this->filter_map == nullptr) this->BuildMaps();
	return(this->UpdateOffsetMapWrapper(offset, this->filter_map, position));
}

uint32_t VariantBlockFooter::AddStreamWrapper(const uint32_t id, map_type* map, header_type*& offsets, uint16_t& stream_counter){
	map_type::const_iterator it = map->find(id);
	if(it == map->end()){
		(*map)[id] = stream_counter;
		offsets[stream_counter].data_header.global_key = id;
		++stream_counter;
	}

	return((*map)[id]);
}

uint32_t VariantBlockFooter::AddInfo(const uint32_t id){
	if(this->info_map == nullptr) this->BuildMaps();
	return(this->AddStreamWrapper(id, this->info_map, this->info_offsets, this->n_info_streams));
}

uint32_t VariantBlockFooter::AddFormat(const uint32_t id){
	if(this->format_map == nullptr) this->BuildMaps();
	return(this->AddStreamWrapper(id, this->format_map, this->format_offsets, this->n_format_streams));
}

uint32_t VariantBlockFooter::AddFilter(const uint32_t id){
	if(this->filter_map == nullptr) this->BuildMaps();
	return(this->AddStreamWrapper(id, this->filter_map, this->filter_offsets, this->n_filter_streams));
}

io::BasicBuffer& operator<<(io::BasicBuffer& buffer, const VariantBlockFooter& entry){
	buffer += (uint16_t)entry.n_info_streams;
	buffer += (uint16_t)entry.n_format_streams;
	buffer += (uint16_t)entry.n_filter_streams;
	buffer += (uint16_t)entry.n_info_patterns;
	buffer += (uint16_t)entry.n_format_patterns;
	buffer += (uint16_t)entry.n_filter_patterns;

	for(uint32_t i = 0; i < YON_BLK_N_STATIC; ++i)       buffer << entry.offsets[i];
	for(uint32_t i = 0; i < entry.n_info_streams; ++i)   buffer << entry.info_offsets[i];
	for(uint32_t i = 0; i < entry.n_format_streams; ++i) buffer << entry.format_offsets[i];
	for(uint32_t i = 0; i < entry.n_filter_streams; ++i) buffer << entry.filter_offsets[i];

	for(uint32_t i = 0; i < entry.n_info_patterns; ++i)   buffer << entry.info_patterns[i];
	for(uint32_t i = 0; i < entry.n_format_patterns; ++i) buffer << entry.format_patterns[i];
	for(uint32_t i = 0; i < entry.n_filter_patterns; ++i) buffer << entry.filter_patterns[i];

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

	for(uint32_t i = 0; i < YON_BLK_N_STATIC; ++i)
		buffer >> entry.offsets[i];

	for(uint32_t i = 0; i < entry.n_info_streams; ++i){
		buffer >> entry.info_offsets[i];
		entry.UpdateInfo(entry.info_offsets[i], i);
	}

	for(uint32_t i = 0; i < entry.n_format_streams; ++i){
		buffer >> entry.format_offsets[i];
		entry.UpdateFormat(entry.format_offsets[i], i);
	}

	for(uint32_t i = 0; i < entry.n_filter_streams; ++i){
		buffer >> entry.filter_offsets[i];
		entry.UpdateFilter(entry.filter_offsets[i], i);
	}

	entry.info_patterns = new yon_blk_bv_pair[entry.n_info_patterns];
	for(uint32_t i = 0; i < entry.n_info_patterns; ++i){
		buffer >> entry.info_patterns[i];
		entry.UpdateInfoPattern(entry.info_patterns[i].pattern, i);
		entry.info_patterns[i].Build(entry.n_info_streams, entry.info_map);
	}

	entry.format_patterns = new yon_blk_bv_pair[entry.n_format_patterns];
	for(uint32_t i = 0; i < entry.n_format_patterns; ++i){
		buffer >> entry.format_patterns[i];
		entry.UpdateFormatPattern(entry.format_patterns[i].pattern, i);
		entry.format_patterns[i].Build(entry.n_format_streams, entry.format_map);
	}

	entry.filter_patterns = new yon_blk_bv_pair[entry.n_filter_patterns];
	for(uint32_t i = 0; i < entry.n_filter_patterns; ++i){
		buffer >> entry.filter_patterns[i];
		entry.UpdateFilterPattern(entry.filter_patterns[i].pattern, i);
		entry.filter_patterns[i].Build(entry.n_filter_streams, entry.filter_map);
	}

	return(buffer);
}

}
}
