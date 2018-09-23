#include "variant_block.h"

#include "third_party/xxhash/xxhash.h"

namespace tachyon{

// vnt controller
yon_vb_hdr_cont::yon_vb_hdr_cont():
	has_gt(0),
	has_gt_permuted(0),
	any_encrypted(0),
	unused(0)
{}
yon_vb_hdr_cont::~yon_vb_hdr_cont(){}

void yon_vb_hdr_cont::clear(){ memset(this, 0, sizeof(uint16_t)); }

std::ostream& operator<<(std::ostream& stream, const yon_vb_hdr_cont& controller){
	const uint16_t c = controller.has_gt |
				  controller.has_gt_permuted << 1 |
				  controller.any_encrypted   << 2 |
				  controller.unused          << 3;

	stream.write(reinterpret_cast<const char*>(&c), sizeof(uint16_t));
	return(stream);
}

std::istream& operator>>(std::istream& stream, yon_vb_hdr_cont& controller){
	uint16_t* c = reinterpret_cast<uint16_t*>(&controller);
	stream.read(reinterpret_cast<char*>(c), sizeof(uint16_t));
	return(stream);
}

// vbblock header
yon_vb_hdr::yon_vb_hdr() :
	l_offset_footer(0),
	block_hash(0),
	contig_id(-1),
	min_position(0),
	max_position(0),
	n_variants(0)
{}

yon_vb_hdr::~yon_vb_hdr(){}

std::ostream& operator<<(std::ostream& stream, const yon_vb_hdr& entry){
	stream.write(reinterpret_cast<const char*>(&entry.l_offset_footer), sizeof(uint32_t));
	stream.write(reinterpret_cast<const char*>(&entry.block_hash),      sizeof(uint64_t));
	stream << entry.controller;
	stream.write(reinterpret_cast<const char*>(&entry.contig_id),       sizeof(uint32_t));
	stream.write(reinterpret_cast<const char*>(&entry.min_position),    sizeof(int64_t));
	stream.write(reinterpret_cast<const char*>(&entry.max_position),    sizeof(int64_t));
	stream.write(reinterpret_cast<const char*>(&entry.n_variants),      sizeof(uint32_t));

	return(stream);
}

std::ifstream& operator>>(std::ifstream& stream, yon_vb_hdr& entry){
	stream.read(reinterpret_cast<char*>(&entry.l_offset_footer), sizeof(uint32_t));
	stream.read(reinterpret_cast<char*>(&entry.block_hash),      sizeof(uint64_t));
	stream >> entry.controller;
	stream.read(reinterpret_cast<char*>(&entry.contig_id),       sizeof(uint32_t));
	stream.read(reinterpret_cast<char*>(&entry.min_position),    sizeof(int64_t));
	stream.read(reinterpret_cast<char*>(&entry.max_position),    sizeof(int64_t));
	stream.read(reinterpret_cast<char*>(&entry.n_variants),      sizeof(uint32_t));

	return(stream);
}

void yon_vb_hdr::reset(void){
	this->l_offset_footer    = 0;
	this->block_hash         = 0;
	this->controller.clear();
	this->contig_id          = -1;
	this->min_position       = 0;
	this->max_position       = 0;
	this->n_variants         = 0;
}

// vblock footer
yon_vb_ftr::yon_vb_ftr():
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

yon_vb_ftr::~yon_vb_ftr(){
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

yon_vb_ftr::yon_vb_ftr(const self_type& other) :
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

yon_vb_ftr::yon_vb_ftr(self_type&& other) noexcept :
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

yon_vb_ftr& yon_vb_ftr::operator=(const self_type& other){
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
	*this = yon_vb_ftr(other); // invoke standard copy ctor followed by move ctor.
	return(*this);
}

yon_vb_ftr& yon_vb_ftr::operator=(self_type&& other) noexcept{
	if(this == &other){
		// precautions against self-moves
		return *this;
	}

	delete [] this->offsets;         this->offsets         = nullptr;
	delete [] this->info_offsets;    this->info_offsets    = nullptr;
	delete [] this->format_offsets;  this->format_offsets  = nullptr;
	delete [] this->filter_offsets;  this->filter_offsets  = nullptr;
	delete [] this->info_patterns;   this->info_patterns   = nullptr;
	delete [] this->format_patterns; this->format_patterns = nullptr;
	delete [] this->filter_patterns; this->filter_patterns = nullptr;
	delete this->info_map;   this->info_map   = nullptr;
	delete this->format_map; this->format_map = nullptr;
	delete this->filter_map; this->filter_map = nullptr;
	delete this->info_pattern_map;   this->info_pattern_map   = nullptr;
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

void yon_vb_ftr::reset(void){
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

void yon_vb_ftr::resetTables(){
	if(this->info_map != nullptr)   this->info_map->clear();
	if(this->format_map != nullptr) this->format_map->clear();
	if(this->filter_map != nullptr) this->filter_map->clear();
	if(this->info_pattern_map != nullptr)   this->info_pattern_map->clear();
	if(this->format_pattern_map != nullptr) this->format_pattern_map->clear();
	if(this->filter_pattern_map != nullptr) this->filter_pattern_map->clear();
}

void yon_vb_ftr::AllocateInfoHeaders(const uint32_t n_info_streams){
	delete [] this->info_offsets;
	if(n_info_streams == 0){
		this->info_offsets = nullptr;
		return;
	}
	this->info_offsets = new header_type[n_info_streams];
}

void yon_vb_ftr::AllocateFormatHeaders(const uint32_t n_format_streams){
	delete [] this->format_offsets;
	if(n_format_streams == 0){
		this->format_offsets = nullptr;
		return;
	}
	this->format_offsets = new header_type[n_format_streams];
}

void yon_vb_ftr::AllocateFilterHeaders(const uint32_t n_filter_streams){
	delete [] this->filter_offsets;
	if(n_filter_streams == 0){
		this->filter_offsets = nullptr;
		return;
	}
	this->filter_offsets = new header_type[n_filter_streams];
}

void yon_vb_ftr::AllocateHeaders(const uint32_t n_info_streams,
                                         const uint32_t n_format_streams,
                                         const uint32_t n_filter_streams)
{
	this->AllocateInfoHeaders(n_info_streams);
	this->AllocateFormatHeaders(n_format_streams);
	this->AllocateFilterHeaders(n_filter_streams);
}

bool yon_vb_ftr::ConstructInfoBitVector(std::unordered_map<uint32_t,uint32_t>* pattern_map){
	for(uint32_t i = 0; i < this->n_info_patterns; ++i){
		this->info_patterns[i].Build(this->n_info_streams, pattern_map);
	}
	return true;
}

bool yon_vb_ftr::ConstructFormatBitVector(std::unordered_map<uint32_t,uint32_t>* pattern_map){
	for(uint32_t i = 0; i < this->n_format_patterns; ++i){
		this->format_patterns[i].Build(this->n_format_streams, pattern_map);
	}
	return true;
}

bool yon_vb_ftr::ConstructFilterBitVector(std::unordered_map<uint32_t,uint32_t>* pattern_map){
	for(uint32_t i = 0; i < this->n_filter_patterns; ++i){
		this->filter_patterns[i].Build(this->n_filter_streams, pattern_map);
	}
	return true;
}

uint32_t yon_vb_ftr::AddPatternWrapper(const std::vector<int>& pattern,
                                               map_pattern_type* pattern_map,
                                               yon_blk_bv_pair* bv_pairs,
                                               uint16_t& stream_counter)
{
	uint64_t pattern_hash = yon_vb_ftr::HashIdentifiers(pattern);
	const map_pattern_type::const_iterator it = pattern_map->find(pattern_hash); // search for pattern
	if(it == pattern_map->end()){
		(*pattern_map)[pattern_hash] = stream_counter;
		bv_pairs[stream_counter].pattern = pattern;
		++stream_counter;
	}

	return((*pattern_map)[pattern_hash]);
}

uint32_t yon_vb_ftr::AddInfoPattern(const std::vector<int>& pattern){
	if(this->info_pattern_map == nullptr) this->AllocatePatternMaps();
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

uint32_t yon_vb_ftr::AddFormatPattern(const std::vector<int>& pattern){
	if(this->format_pattern_map == nullptr) this->AllocatePatternMaps();
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

uint32_t yon_vb_ftr::AddFilterPattern(const std::vector<int>& pattern){
	if(this->filter_pattern_map == nullptr) this->AllocatePatternMaps();
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
uint32_t yon_vb_ftr::UpdatePatternMapWrapper(const std::vector<int>& pattern,
                                                  map_pattern_type* pattern_map,
                                                  const uint16_t& local_position)
{
	uint64_t pattern_hash = yon_vb_ftr::HashIdentifiers(pattern);
	const map_pattern_type::const_iterator it = pattern_map->find(pattern_hash); // search for pattern
	if(it == pattern_map->end())
		(*pattern_map)[pattern_hash] = local_position;

	return((*pattern_map)[pattern_hash]);
}

uint32_t yon_vb_ftr::UpdateInfoPatternMap(const std::vector<int>& pattern, const uint16_t local_position){
	if(this->info_pattern_map == nullptr) this->AllocatePatternMaps();
	if(this->n_info_patterns_allocated == 0){
		delete [] this->info_patterns;
		this->info_patterns = new yon_blk_bv_pair[100];
		this->n_info_patterns_allocated = 100;
	}
	return(this->UpdatePatternMapWrapper(pattern, this->info_pattern_map, local_position));
}

uint32_t yon_vb_ftr::UpdateFormatPatternMap(const std::vector<int>& pattern, const uint16_t local_position){
	if(this->format_pattern_map == nullptr) this->AllocatePatternMaps();
	if(this->n_format_patterns_allocated == 0){
		delete [] this->format_patterns;
		this->format_patterns = new yon_blk_bv_pair[100];
		this->n_format_patterns_allocated = 100;
	}
	return(this->UpdatePatternMapWrapper(pattern, this->format_pattern_map, local_position));
}

uint32_t yon_vb_ftr::UpdateFilterPatternMap(const std::vector<int>& pattern, const uint16_t local_position){
	if(this->filter_pattern_map == nullptr) this->AllocatePatternMaps();
	if(this->n_filter_patterns_allocated == 0){
		delete [] this->filter_patterns;
		this->filter_patterns = new yon_blk_bv_pair[100];
		this->n_filter_patterns_allocated = 100;
	}
	return(this->UpdatePatternMapWrapper(pattern, this->filter_pattern_map, local_position));
}

void yon_vb_ftr::Finalize(void){
	this->ConstructInfoBitVector(this->info_map);
	this->ConstructFormatBitVector(this->format_map);
	this->ConstructFilterBitVector(this->filter_map);
}

bool yon_vb_ftr::AllocateMaps(void){
	delete this->info_map;
	delete this->filter_map;
	delete this->format_map;

	this->info_map   = new map_type();
	this->filter_map = new map_type();
	this->format_map = new map_type();

	return true;
}

bool yon_vb_ftr::AllocatePatternMaps(void){
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

uint32_t yon_vb_ftr::UpdateOffsetMapWrapper(const header_type& offset,
                                                    map_type* map,
                                                    const uint16_t& local_position)
{
	assert(map != nullptr);
	map_type::const_iterator it = map->find(offset.data_header.global_key);
	if(it == map->end())
		(*map)[offset.data_header.global_key] = local_position;

	return((*map)[offset.data_header.global_key]);
}

uint32_t yon_vb_ftr::UpdateInfoMap(const header_type& offset, const uint16_t local_position){
	if(this->info_map == nullptr) this->AllocateMaps();
	return(this->UpdateOffsetMapWrapper(offset, this->info_map, local_position));
}

uint32_t yon_vb_ftr::UpdateFormatMap(const header_type& offset, const uint16_t local_position){
	if(this->format_map == nullptr) this->AllocateMaps();
	return(this->UpdateOffsetMapWrapper(offset, this->format_map, local_position));
}

uint32_t yon_vb_ftr::UpdateFilterMap(const header_type& offset, const uint16_t local_position){
	if(this->filter_map == nullptr) this->AllocateMaps();
	return(this->UpdateOffsetMapWrapper(offset, this->filter_map, local_position));
}

uint32_t yon_vb_ftr::AddStreamWrapper(const uint32_t global_id,
                                              map_type* map,
                                              header_type*& offsets,
                                              uint16_t& n_streams)
{
	map_type::const_iterator it = map->find(global_id);
	if(it == map->end()){
		(*map)[global_id] = n_streams;
		offsets[n_streams].data_header.global_key = global_id;
		++n_streams;
	}

	return((*map)[global_id]);
}

uint32_t yon_vb_ftr::AddInfo(const uint32_t global_id){
	if(this->info_map == nullptr) this->AllocateMaps();
	return(this->AddStreamWrapper(global_id, this->info_map, this->info_offsets, this->n_info_streams));
}

uint32_t yon_vb_ftr::AddFormat(const uint32_t global_id){
	if(this->format_map == nullptr) this->AllocateMaps();
	return(this->AddStreamWrapper(global_id, this->format_map, this->format_offsets, this->n_format_streams));
}

uint32_t yon_vb_ftr::AddFilter(const uint32_t global_id){
	if(this->filter_map == nullptr) this->AllocateMaps();
	return(this->AddStreamWrapper(global_id, this->filter_map, this->filter_offsets, this->n_filter_streams));
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const yon_vb_ftr& entry){
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


yon_buffer_t& operator>>(yon_buffer_t& buffer, yon_vb_ftr& entry){
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

	entry.AllocateMaps(); // Construct new maps.
	entry.AllocatePatternMaps(); // Construct new pattern maps.
	entry.offsets        = new yon_dc_hdr[YON_BLK_N_STATIC];
	entry.info_offsets   = new yon_dc_hdr[entry.n_info_streams];
	entry.format_offsets = new yon_dc_hdr[entry.n_format_streams];
	entry.filter_offsets = new yon_dc_hdr[entry.n_filter_streams];
	entry.n_info_patterns_allocated   = entry.n_info_streams;
	entry.n_format_patterns_allocated = entry.n_format_streams;
	entry.n_filter_patterns_allocated = entry.n_filter_streams;

	for(uint32_t i = 0; i < YON_BLK_N_STATIC; ++i)
		buffer >> entry.offsets[i];

	for(uint32_t i = 0; i < entry.n_info_streams; ++i){
		buffer >> entry.info_offsets[i];
		entry.UpdateInfoMap(entry.info_offsets[i], i);
	}

	for(uint32_t i = 0; i < entry.n_format_streams; ++i){
		buffer >> entry.format_offsets[i];
		entry.UpdateFormatMap(entry.format_offsets[i], i);
	}

	for(uint32_t i = 0; i < entry.n_filter_streams; ++i){
		buffer >> entry.filter_offsets[i];
		entry.UpdateFilterMap(entry.filter_offsets[i], i);
	}

	entry.info_patterns = new yon_blk_bv_pair[entry.n_info_patterns];
	for(uint32_t i = 0; i < entry.n_info_patterns; ++i){
		buffer >> entry.info_patterns[i];
		entry.UpdateInfoPatternMap(entry.info_patterns[i].pattern, i);
		entry.info_patterns[i].Build(entry.n_info_streams, entry.info_map);
	}

	entry.format_patterns = new yon_blk_bv_pair[entry.n_format_patterns];
	for(uint32_t i = 0; i < entry.n_format_patterns; ++i){
		buffer >> entry.format_patterns[i];
		entry.UpdateFormatPatternMap(entry.format_patterns[i].pattern, i);
		entry.format_patterns[i].Build(entry.n_format_streams, entry.format_map);
	}

	entry.filter_patterns = new yon_blk_bv_pair[entry.n_filter_patterns];
	for(uint32_t i = 0; i < entry.n_filter_patterns; ++i){
		buffer >> entry.filter_patterns[i];
		entry.UpdateFilterPatternMap(entry.filter_patterns[i].pattern, i);
		entry.filter_patterns[i].Build(entry.n_filter_streams, entry.filter_map);
	}

	return(buffer);
}

uint64_t yon_vb_ftr::HashIdentifiers(const std::vector<int>& id_vector){
	XXH64_state_t* const state = XXH64_createState();
	if (state==NULL) abort();

	XXH_errorcode const resetResult = XXH64_reset(state, 71236251);
	if (resetResult == XXH_ERROR) abort();

	for(uint32_t i = 0; i < id_vector.size(); ++i){
		XXH_errorcode const addResult = XXH64_update(state, (const void*)&id_vector[i], sizeof(int));
		if (addResult == XXH_ERROR) abort();
	}

	uint64_t hash = XXH64_digest(state);
	XXH64_freeState(state);

	return hash;
}

// Variant block settings
yon_vb_settings::yon_vb_settings() :
	show_vcf_header(true),
	display_ref(true),
	display_alt(true),
	display_filter(true),
	construct_occ_table(false),
	annotate_extra(false),
	load_static(0),
	display_static(std::numeric_limits<uint32_t>::max())
{}

yon_vb_settings& yon_vb_settings::LoadWrapper(bool set, const int field_bv){
	this->load_static &= ~(field_bv);
	if(set) this->load_static |= field_bv;
	return(*this);
}

yon_vb_settings& yon_vb_settings::DisplayWrapper(bool set, const int field_bv){
	this->display_static &= ~(field_bv);
	if(set) this->display_static |= field_bv;
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadDisplayWrapper(bool set, const int field_bv){
	this->LoadWrapper(set, field_bv);
	this->DisplayWrapper(set, field_bv);
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadCore(const bool set){
	for(uint32_t i = YON_BLK_CONTIG; i <= YON_BLK_ID_FILTER; ++i){
		const uint32_t bv = 1 << i;
		this->LoadWrapper(set, bv);
	}
	return(*this);
}

yon_vb_settings& yon_vb_settings::DisplayCore(const bool set){
	for(uint32_t i = YON_BLK_CONTIG; i <= YON_BLK_ID_FILTER; ++i){
		const uint32_t bv = 1 << i;
		this->DisplayWrapper(set, bv);
	}
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadAll(const bool set){
	if(set){
		this->load_static = std::numeric_limits<uint32_t>::max();
	} else {
		this->load_static = 0;
	}
	return(*this);
}

yon_vb_settings& yon_vb_settings::DisplayAll(const bool set){
	if(set){
		this->display_static = std::numeric_limits<uint32_t>::max();
	} else {
		this->display_static = 0;
	}
	this->display_alt = set;
	this->display_ref = set;
	this->display_filter = set;
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadAllMeta(const bool set){
	for(uint32_t i = YON_BLK_CONTIG; i <= YON_BLK_ID_FILTER; ++i){
		const uint32_t bv = 1 << i;
		this->LoadWrapper(set, i);
	}
	return(*this);
}

yon_vb_settings& yon_vb_settings::DisplayAllMeta(const bool set){
	for(uint32_t i = YON_BLK_CONTIG; i <= YON_BLK_ID_FILTER; ++i){
		const uint32_t bv = 1 << i;
		this->DisplayWrapper(set, i);
	}

	this->display_alt = set;
	this->display_ref = set;
	this->display_filter = set;
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadAllFilter(const bool set){
	this->LoadWrapper(set, YON_BLK_BV_ID_INFO);
	this->LoadWrapper(set, YON_BLK_BV_ID_FORMAT);
	this->LoadWrapper(set, YON_BLK_BV_ID_FILTER);
	this->LoadWrapper(set, YON_BLK_BV_CONTROLLER);

	return(*this);
}

yon_vb_settings& yon_vb_settings::DisplayAllFilter(const bool set){
	this->DisplayWrapper(set, YON_BLK_BV_ID_INFO);
	this->DisplayWrapper(set, YON_BLK_BV_ID_FORMAT);
	this->DisplayWrapper(set, YON_BLK_BV_ID_FILTER);
	this->DisplayWrapper(set, YON_BLK_BV_CONTROLLER);

	this->display_filter  = true;
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadAllInfo(const bool set){
	this->LoadWrapper(set, YON_BLK_BV_INFO); // all info
	if(set) this->LoadMinimumVcf(true);

	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadInfo(const std::string& field_name){
	if(field_name.size() == 0) return(*this);
	this->info_list.push_back(field_name);
	this->LoadMinimumVcf(true);
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadInfo(const uint32_t field_id){
	this->info_id_global.push_back(field_id);
	this->LoadMinimumVcf(true);
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadGenotypes(const bool set){
	this->LoadWrapper(set, YON_BLK_BV_GT);
	this->LoadWrapper(set, YON_BLK_BV_GT_SUPPORT);
	this->LoadWrapper(set, YON_BLK_BV_GT_PLOIDY);
	this->LoadWrapper(set, YON_BLK_BV_PPA);
	return(*this);
}

yon_vb_settings& yon_vb_settings::DisplayGenotypes(const bool set){
	this->DisplayWrapper(set, YON_BLK_BV_GT);
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadPermutationArray(const bool set){
	this->LoadWrapper(set, YON_BLK_BV_PPA);
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadAllFormat(const bool set){
	this->LoadGenotypes(set);
	this->LoadWrapper(set, YON_BLK_BV_FORMAT); // all format
	if(set) this->LoadCore(set);
	return(*this);
}

yon_vb_settings& yon_vb_settings::DisplayAllFormat(const bool set){
	this->DisplayGenotypes(set);
	this->DisplayWrapper(set, YON_BLK_BV_FORMAT); // all format
	return(*this);
}

yon_vb_settings& yon_vb_settings::DisplayAllInfo(const bool set){
	this->DisplayWrapper(set, YON_BLK_BV_INFO);
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadFormat(const std::string& field_name){
	if(field_name.size() == 0) return(*this);
	this->LoadMinimumVcf(true);
	if(field_name == "GT") this->LoadGenotypes(true);
	this->format_list.push_back(field_name);
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadFormat(const uint32_t field_id){
	this->LoadMinimumVcf(true);
	this->format_id_global.push_back(field_id);
	return(*this);
}

yon_vb_settings& yon_vb_settings::LoadMinimumVcf(const bool set){
	this->LoadWrapper(set, YON_BLK_BV_CONTIG);
	this->LoadWrapper(set, YON_BLK_BV_POSITION);
	this->LoadWrapper(set, YON_BLK_BV_CONTROLLER);
	this->LoadWrapper(set, YON_BLK_BV_ID_INFO);
	this->LoadWrapper(set, YON_BLK_BV_ID_FORMAT);
	this->LoadWrapper(set, YON_BLK_BV_ID_FILTER);
	this->LoadWrapper(set, YON_BLK_BV_ALLELES);
	this->LoadWrapper(set, YON_BLK_BV_REFALT);
	return(*this);
}

yon_vb_settings& yon_vb_settings::DisplayMinimumVcf(const bool set){
	this->DisplayWrapper(set, YON_BLK_BV_CONTIG);
	this->DisplayWrapper(set, YON_BLK_BV_POSITION);
	return(*this);
}

bool yon_vb_settings::Parse(const header_type& header){
	std::regex field_identifier_regex("^[A-Za-z_0-9]{1,}$");

	for(uint32_t i = 0; i < this->info_list.size(); ++i){
		std::vector<std::string> ind = utility::split(this->info_list[i], ',');
		for(uint32_t j = 0; j < ind.size(); ++j){
			ind[j] = utility::remove_excess_whitespace(ind[j]);
			if(std::regex_match(ind[j], field_identifier_regex)){
				const VcfInfo* info = header.GetInfo(ind[j]);
				if(info == nullptr){
					std::cerr << utility::timestamp("ERROR") << "Cannot find INFO field: " << ind[j] << " in string " << this->info_list[i] << std::endl;
					continue;
				}
				this->LoadInfo(ind[j]);
			} else {
				std::cerr << utility::timestamp("ERROR") << "Illegal field name: " << ind[j] << ". Must match \"[A-Za-z_0-9]\"..." << std::endl;
				return(false);
			}
			this->LoadInfo(ind[j]);
		}
	}

	// Todo
	//for(uint32_t i = 0; i < this->format_list.size(); ++i){
	//
	//}

	return true;
}

bool yon_vb_settings::ParseCommandString(const std::vector<std::string>& command, const header_type& header){
	bool allGood = true;

	this->display_static = 0;
	this->load_static = 0;

	std::regex field_identifier_regex("^[A-Za-z_0-9]{1,}$");
	for(uint32_t i = 0; i < command.size(); ++i){
		std::vector<std::string> partitions = utility::split(command[i], ';');
		for(uint32_t p = 0; p < partitions.size(); ++p){
			partitions[p].erase(std::remove(partitions[p].begin(), partitions[p].end(), ' '), partitions[p].end()); // remove all spaces
			if(strncasecmp(partitions[p].data(), "INFO=", 5) == 0){
				std::vector<std::string> ind = utility::split(partitions[p].substr(5,command.size()-5), ',');
				for(uint32_t j = 0; j < ind.size(); ++j){
					ind[j] = utility::remove_excess_whitespace(ind[j]);
					if(std::regex_match(ind[j], field_identifier_regex)){
						const VcfInfo* info = header.GetInfo(ind[j]);
						if(info == nullptr){
							std::cerr << utility::timestamp("ERROR") << "Cannot find INFO field: " << ind[j] << " in string " << partitions[p] << std::endl;
							allGood = false;
							continue;
						}
						this->LoadInfo(info->idx);
						this->DisplayAllInfo(true);
					} else {
						std::cerr << utility::timestamp("ERROR") << "Illegal field name: " << ind[j] << ". Must match \"[A-Za-z_0-9]\"..." << std::endl;
						allGood = false;
					}
				}
			} else if(strncasecmp(partitions[p].data(), "INFO", 4) == 0 && partitions[p].size() == 4){
				this->LoadAllInfo(true);
				this->DisplayAllInfo(true);
			} else if(strncasecmp(partitions[p].data(), "FORMAT", 6) == 0 && partitions[p].size() == 6){
				this->LoadAllFormat(true);
				this->DisplayAllFormat(true);

			} else if(strncasecmp(partitions[p].data(), "FORMAT=", 7) == 0){
				std::vector<std::string> ind = utility::split(partitions[p].substr(7,command.size()-7), ',');
				for(uint32_t j = 0; j < ind.size(); ++j){
					std::transform(ind[j].begin(), ind[j].end(), ind[j].begin(), ::toupper); // transform to UPPERCASE
					if(std::regex_match(ind[j], field_identifier_regex)){
						// Special case for genotypes
						if(strncasecmp(ind[j].data(), "GT", 2) == 0 && ind[j].size() == 2){
							this->LoadMinimumVcf(true);
							this->LoadGenotypes(true);

							const VcfFormat* fmt = header.GetFormat(ind[j]);
							if(fmt == nullptr){
								std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							this->LoadFormat(fmt->idx);
							this->DisplayAllFormat(true);

						} else if(strncasecmp(ind[j].data(), "GENOTYPES", 9) == 0 && ind[j].size() == 9){
							this->LoadMinimumVcf(true);
							this->LoadGenotypes(true);
							this->DisplayAllFormat(true);

							const VcfFormat* fmt = header.GetFormat(ind[j]);
							if(fmt == nullptr){
								std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							this->LoadFormat(fmt->idx);
							this->DisplayAllFormat(true);
						}
						// Any other FORMAT
						else {
							const VcfFormat* fmt = header.GetFormat(ind[j]);
							if(fmt == nullptr){
								std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							this->LoadFormat(fmt->idx);
							this->DisplayAllFormat(true);
						}
					} else {
						std::cerr << utility::timestamp("ERROR") << "Cannot find FORMAT field: " << ind[j] << " in string " << partitions[p] << std::endl;
						allGood = false;
					}
				}

			} else if((strncasecmp(partitions[p].data(), "CONTIG", 6) == 0 && partitions[p].length() == 6) ||
					  (strncasecmp(partitions[p].data(), "CHROM", 5) == 0 && partitions[p].length() == 5)  ||
					  (strncasecmp(partitions[p].data(), "CHROMOSOME", 10) == 0 && partitions[p].length() == 10)){
				this->LoadWrapper(true, YON_BLK_BV_CONTIG);
				this->DisplayWrapper(true, YON_BLK_BV_CONTIG);
			} else if((strncasecmp(partitions[p].data(), "POSITION", 8) == 0 && partitions[p].length() == 8) ||
					  (strncasecmp(partitions[p].data(), "POS", 3) == 0 && partitions[p].length() == 3)){
				this->LoadWrapper(true, YON_BLK_BV_POSITION);
			} else if((strncasecmp(partitions[p].data(), "REF", 3) == 0 && partitions[p].length() == 3) ||
					  (strncasecmp(partitions[p].data(), "REFERENCE", 9) == 0 && partitions[p].length() == 9)){
				this->LoadWrapper(true, YON_BLK_BV_ALLELES);
				this->LoadWrapper(true, YON_BLK_BV_REFALT);
				this->LoadWrapper(true, YON_BLK_BV_CONTROLLER);
				this->DisplayWrapper(true, YON_BLK_BV_ALLELES);
				this->DisplayWrapper(true, YON_BLK_BV_REFALT);
				this->DisplayWrapper(true, YON_BLK_BV_CONTROLLER);
				this->display_ref = true;
			} else if((strncasecmp(partitions[p].data(), "ALT", 3) == 0 && partitions[p].length() == 3) ||
					  (strncasecmp(partitions[p].data(), "ALTERNATE", 9) == 0 && partitions[p].length() == 9)){
				this->LoadWrapper(true, YON_BLK_BV_ALLELES);
				this->LoadWrapper(true, YON_BLK_BV_REFALT);
				this->LoadWrapper(true, YON_BLK_BV_CONTROLLER);
				this->DisplayWrapper(true, YON_BLK_BV_ALLELES);
				this->DisplayWrapper(true, YON_BLK_BV_REFALT);
				this->DisplayWrapper(true, YON_BLK_BV_CONTROLLER);
				this->display_alt = true;
			} else if((strncasecmp(partitions[p].data(), "QUALITY", 7) == 0 && partitions[p].length() == 7) ||
					  (strncasecmp(partitions[p].data(), "QUAL", 4) == 0 && partitions[p].length() == 4)){
				this->LoadWrapper(true, YON_BLK_BV_QUALITY);
				this->DisplayWrapper(true, YON_BLK_BV_QUALITY);
			} else if((strncasecmp(partitions[p].data(), "NAMES", 5) == 0 && partitions[p].length() == 5) ||
					  (strncasecmp(partitions[p].data(), "NAME", 4) == 0 && partitions[p].length() == 4)){
				this->LoadWrapper(true, YON_BLK_BV_NAMES);
				this->DisplayWrapper(true, YON_BLK_BV_NAMES);
			} else if((strncasecmp(partitions[p].data(), "FILTERS", 7) == 0 && partitions[p].length() == 7) ||
					  (strncasecmp(partitions[p].data(), "FILTER", 6) == 0 && partitions[p].length() == 6)){
				this->LoadWrapper(true, YON_BLK_BV_CONTROLLER);
				this->LoadWrapper(true, YON_BLK_BV_ID_FILTER);
				this->DisplayWrapper(true, YON_BLK_BV_CONTROLLER);
				this->DisplayWrapper(true, YON_BLK_BV_ID_FILTER);
				this->display_filter = true;
			} else {
				std::cerr << utility::timestamp("ERROR") << "Unknown pattern: " << partitions[p] << std::endl;
				allGood = false;
			}
		}
	}

	if(allGood == false) return false;
	return true;
}

yon1_vb_t::yon1_vb_t() :
	n_info_c_allocated(0),
	n_format_c_allocated(0),
	base_containers(new container_type[YON_BLK_N_STATIC]),
	info_containers(nullptr),
	format_containers(nullptr),
	gt_ppa(nullptr),
	load_settings(nullptr),
	end_block_(0),
	start_compressed_data_(0),
	end_compressed_data_(0)
{
	this->base_containers[YON_BLK_ALLELES].SetType(YON_TYPE_STRUCT);
	this->base_containers[YON_BLK_CONTROLLER].SetType(YON_TYPE_16B);
	this->base_containers[YON_BLK_REFALT].SetType(YON_TYPE_8B);
	this->footer_support.resize(65536);
}

yon1_vb_t::yon1_vb_t(const uint16_t n_info, const uint16_t n_format) :
	n_info_c_allocated(n_info),
	n_format_c_allocated(n_format),
	base_containers(new container_type[YON_BLK_N_STATIC]),
	info_containers(new container_type[n_info]),
	format_containers(new container_type[n_format]),
	gt_ppa(nullptr),
	load_settings(nullptr),
	end_block_(0),
	start_compressed_data_(0),
	end_compressed_data_(0)
{
	this->base_containers[YON_BLK_ALLELES].SetType(YON_TYPE_STRUCT);
	this->base_containers[YON_BLK_CONTROLLER].SetType(YON_TYPE_16B);
	this->base_containers[YON_BLK_REFALT].SetType(YON_TYPE_8B);
	this->footer_support.resize(65536);
}

yon1_vb_t::~yon1_vb_t(){
	delete [] this->base_containers;
	delete [] this->info_containers;
	delete [] this->format_containers;
	delete this->gt_ppa;
	delete this->load_settings;
}

yon1_vb_t::yon1_vb_t(const self_type& other) :
	n_info_c_allocated(other.n_info_c_allocated),
	n_format_c_allocated(other.n_info_c_allocated),
	header(other.header),
	footer(other.footer),
	base_containers(new container_type[YON_BLK_N_STATIC]),
	info_containers(new container_type[other.n_info_c_allocated]),
	format_containers(new container_type[other.n_format_c_allocated]),
	gt_ppa(nullptr),
	load_settings(nullptr),
	end_block_(other.end_block_),
	start_compressed_data_(other.start_compressed_data_),
	end_compressed_data_(other.end_compressed_data_),
	footer_support(other.footer_support)
{
	if(other.gt_ppa != nullptr){
		// Copy ppa data to new object.
		this->gt_ppa = new yon_gt_ppa(*other.gt_ppa);
	}

	if(other.load_settings != nullptr){
		this->load_settings = new yon_blk_load_settings(*other.load_settings);
	}

	for(uint32_t i = 0; i < YON_BLK_N_STATIC; ++i) this->base_containers[i] = other.base_containers[i];
	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i)   this->info_containers[i]   = other.info_containers[i];
	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i) this->format_containers[i] = other.format_containers[i];
}

yon1_vb_t::yon1_vb_t(self_type&& other) noexcept :
	n_info_c_allocated(other.n_info_c_allocated),
	n_format_c_allocated(other.n_format_c_allocated),
	header(std::move(other.header)),
	footer(std::move(other.footer)),
	base_containers(nullptr),
	info_containers(nullptr),
	format_containers(nullptr),
	gt_ppa(nullptr),
	load_settings(nullptr),
	end_block_(other.end_block_),
	start_compressed_data_(other.start_compressed_data_),
	end_compressed_data_(other.end_compressed_data_),
	footer_support(std::move(other.footer_support))
{
	std::swap(this->base_containers, other.base_containers);
	std::swap(this->info_containers, other.info_containers);
	std::swap(this->format_containers, other.format_containers);
	std::swap(this->gt_ppa, other.gt_ppa);
	std::swap(this->load_settings, other.load_settings);
}

yon1_vb_t& yon1_vb_t::operator=(const self_type& other){
	delete [] this->base_containers;
	delete [] this->info_containers;
	delete [] this->format_containers;
	delete this->gt_ppa;
	delete this->load_settings;
	*this = yon1_vb_t(other); // invoke copy ctor
	return(*this);
}

yon1_vb_t& yon1_vb_t::operator=(self_type&& other) noexcept{
	if(this == &other){
		// precautions against self-moves
		return *this;
	}

	this->n_info_c_allocated   = other.n_info_c_allocated;
	this->n_format_c_allocated = other.n_format_c_allocated;
	this->header = std::move(other.header);
	this->footer = std::move(other.footer);
	delete [] this->base_containers; this->base_containers = nullptr;
	std::swap(this->base_containers, other.base_containers);
	delete [] this->info_containers; this->info_containers = nullptr;
	std::swap(this->info_containers, other.info_containers);
	delete [] this->format_containers; this->format_containers = nullptr;
	std::swap(this->format_containers, other.format_containers);
	delete this->gt_ppa; this->gt_ppa = nullptr;
	std::swap(this->gt_ppa, other.gt_ppa);
	this->end_block_ = other.end_block_;
	this->start_compressed_data_ = other.start_compressed_data_;
	this->end_compressed_data_ = other.end_compressed_data_;
	this->footer_support = std::move(other.footer_support);
	delete this->load_settings; this->load_settings = nullptr;
	std::swap(this->load_settings, other.load_settings);
	return(*this);
}

void yon1_vb_t::Allocate(const uint16_t n_info,
                            const uint16_t n_format,
                            const uint16_t n_filter)
{
	// Allocate space for INFO containers.
	delete [] this->info_containers;
	this->info_containers = new container_type[n_info];
	this->n_info_c_allocated = n_info;

	// Allocate space for FORMAT containers.
	delete [] this->format_containers;
	this->format_containers = new container_type[n_format];
	this->n_format_c_allocated = n_format;

	// Alocate space for headers.
	this->footer.AllocateHeaders(n_info, n_format, n_filter);
}

void yon1_vb_t::clear(void){
	for(uint32_t i = 0; i < YON_BLK_N_STATIC; ++i) this->base_containers[i].reset();
	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i)   this->info_containers[i].reset();
	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i) this->format_containers[i].reset();

	this->base_containers[YON_BLK_ALLELES].SetType(YON_TYPE_STRUCT);
	this->base_containers[YON_BLK_CONTROLLER].SetType(YON_TYPE_16B);
	this->base_containers[YON_BLK_REFALT].SetType(YON_TYPE_8B);

	this->end_block_             = 0;
	this->start_compressed_data_ = 0;
	this->end_compressed_data_   = 0;

	this->header.reset();
	this->footer.reset();
	this->footer_support.reset();

	if(this->gt_ppa != nullptr) this->gt_ppa->reset();
}

void yon1_vb_t::resize(const uint32_t s){
	if(s == 0) return;

	for(uint32_t i = 0; i < YON_BLK_N_STATIC; ++i)     this->base_containers[i].resize(s);
	for(uint32_t i = 0; i < n_info_c_allocated; ++i)   this->info_containers[i].resize(s);
	for(uint32_t i = 0; i < n_format_c_allocated; ++i) this->format_containers[i].resize(s);
}

void yon1_vb_t::UpdateContainers(const uint32_t n_samples){
	this->base_containers[YON_BLK_CONTIG].UpdateContainer();
	this->base_containers[YON_BLK_POSITION].UpdateContainer();
	this->base_containers[YON_BLK_REFALT].UpdateContainer(false, true);
	this->base_containers[YON_BLK_QUALITY].UpdateContainer();
	this->base_containers[YON_BLK_NAMES].UpdateContainer(false, true);
	this->base_containers[YON_BLK_ALLELES].UpdateContainer(false, true);
	this->base_containers[YON_BLK_ID_FILTER].UpdateContainer();
	this->base_containers[YON_BLK_ID_FORMAT].UpdateContainer();
	this->base_containers[YON_BLK_ID_INFO].UpdateContainer();
	this->base_containers[YON_BLK_GT_SUPPORT].UpdateContainer();
	this->base_containers[YON_BLK_GT_PLOIDY].UpdateContainer();
	this->base_containers[YON_BLK_CONTROLLER].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_INT8].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_INT16].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_INT32].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_INT64].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_S_INT8].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_S_INT16].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_S_INT32].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_S_INT64].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_N_INT8].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_N_INT16].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_N_INT32].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_N_INT64].UpdateContainer(false, true);

	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i){
		this->info_containers[i].UpdateContainer();
	}

	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i){
		// Illegal to have BOOLEAN fields in the Vcf:Format column.
		// Therefore we assert that this is never the case.
		assert(this->format_containers[i].header.data_header.stride != 0);
		this->format_containers[i].UpdateContainerFormat(true, true, n_samples);
	}
}

bool yon1_vb_t::ReadHeaderFooter(std::ifstream& stream){
	if(!stream.good()){
		std::cerr << utility::timestamp("ERROR") << "File stream is corrupted..." << std::endl;
		return false;
	}

	stream >> this->header; // load header
	this->start_compressed_data_ = (uint64_t)stream.tellg(); // start of compressed data
	stream.seekg(this->start_compressed_data_ + this->header.l_offset_footer); // seek to start of footer
	this->end_compressed_data_   = stream.tellg(); // end of compressed data

	assert(stream.good());

	uint32_t footer_uLength = 0;
	uint32_t footer_cLength = 0;
	uint8_t footer_crc[MD5_DIGEST_LENGTH];
	utility::DeserializePrimitive(footer_uLength, stream);
	utility::DeserializePrimitive(footer_cLength, stream);
	stream.read(reinterpret_cast<char*>(&footer_crc[0]), MD5_DIGEST_LENGTH);

	this->footer_support.resize(footer_cLength);
	stream.read(this->footer_support.data.data(), footer_cLength);
	this->footer_support.data.n_chars_ = footer_cLength;
	this->footer_support.data_uncompressed.resize(footer_uLength);
	this->footer_support.data_uncompressed.n_chars_     = footer_uLength;
	this->footer_support.header.data_header.controller.encoder = YON_ENCODE_ZSTD;
	this->footer_support.header.data_header.cLength            = footer_cLength;
	this->footer_support.header.data_header.uLength            = footer_uLength;
	memcpy(&this->footer_support.header.data_header.crc[0], &footer_crc[0], MD5_DIGEST_LENGTH);

	// Assert end-of-block marker
	uint64_t eof_marker;
	utility::DeserializePrimitive(eof_marker, stream);
	assert(eof_marker == TACHYON_BLOCK_EOF);
	this->end_block_ = stream.tellg(); // end-of-block offset
	stream.seekg(this->start_compressed_data_);
	return(stream.good());
}

bool yon1_vb_t::read(std::ifstream& stream){
	if(this->header.controller.has_gt_permuted && this->header.controller.has_gt){
		stream.seekg(this->start_compressed_data_ + this->footer.offsets[YON_BLK_PPA].data_header.offset);
		stream >> this->base_containers[YON_BLK_PPA];
	}

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i)
		this->LoadContainer(stream, this->footer.offsets[i], this->base_containers[i]);

	// Load all INFO
	delete [] this->info_containers;
	this->info_containers = new container_type[this->footer.n_info_streams];
	this->n_info_c_allocated = this->footer.n_info_streams;
	if(this->footer.n_info_streams){
		stream.seekg(this->start_compressed_data_ + this->footer.info_offsets[0].data_header.offset);
		for(uint32_t i = 0; i < this->footer.n_info_streams; ++i)
			this->LoadContainer(stream, this->footer.info_offsets[i], this->info_containers[i]);

	}

	// Load all FORMAT
	delete [] this->format_containers;
	this->format_containers = new container_type[this->footer.n_format_streams];
	this->n_format_c_allocated = this->footer.n_format_streams;
	if(this->footer.n_format_streams){
		stream.seekg(this->start_compressed_data_ + this->footer.format_offsets[0].data_header.offset);
		for(uint32_t i = 0; i < this->footer.n_format_streams; ++i)
			this->LoadContainer(stream, this->footer.format_offsets[i], this->format_containers[i]);

		// EOF assertion
		assert(this->end_compressed_data_ == (uint64_t)stream.tellg());
	}

	stream.seekg(this->end_block_); // seek to end-of-block
	return(true);
}


bool yon1_vb_t::ParseSettings(yon_vb_settings& settings, const yon_vnt_hdr_t& header){
	// Clear previous information (if any).
	this->load_settings->clear();

	// Construct a black-list of INFO fields that should not be loaded
	// or displayed as they are being re-calculated internally and
	// emitted. If a blacklisted field is available in this block then
	// add that tag to the map.
	/*
	std::unordered_map<uint32_t, std::string> blocked_list;
	if(settings.annotate_extra){
		for(uint32_t i = 0; i < YON_GT_ANNOTATE_FIELDS.size(); ++i){
			const YonInfo* info = header.GetInfo(YON_GT_ANNOTATE_FIELDS[i]);
			if(info != nullptr){
				blocked_list[info->idx] = YON_GT_ANNOTATE_FIELDS[i];
			}
		}
	}
	*/

	// Parse Info. If all Info containers are loaded then we simply copy
	// the order in which they occur. If we are provided with a vector
	// of target global identifiers we have to first map these to the
	// (possible) local identifiers.
	if(settings.load_static & YON_BLK_BV_INFO){
		for(uint32_t i = 0; i < this->footer.n_info_streams; ++i){
			//const std::unordered_map<uint32_t, std::string>::const_iterator it = blocked_list.find(this->footer.info_offsets[i].data_header.global_key);
			//if(it == blocked_list.end()){
				//std::cerr << "adding not blocked" << std::endl;
				this->load_settings->info_id_local_loaded.push_back(i);
				this->load_settings->info_id_global_loaded.push_back(this->footer.info_offsets[i].data_header.global_key);
				this->load_settings->info_map_global[this->load_settings->info_id_global_loaded[i]] = i;
			//} else {
				//std::cerr << "skipping blocked" << std::endl;
			//}
		}
		//std::cerr << this->info_id_local_loaded.size() << "," << this->info_id_global_loaded.size() << std::endl;
	} else {
		std::vector<int> local_ids;
		std::vector<int> global_ids;
		for(uint32_t i = 0; i < settings.info_id_global.size(); ++i){
			// Searches for the global Vcf:INFO idx value in the block. If
			// it is found then return that local idx otherwise -1. If the
			// idx is found store it in the loaded idx vector.
			//const std::unordered_map<uint32_t, std::string>::const_iterator it = blocked_list.find(this->footer.info_offsets[i].data_header.global_key);
			//if(it == blocked_list.end()){
				const int local = this->GetInfoPosition(settings.info_id_global[i]);
				if(local >= 0){
					local_ids.push_back(local);
					global_ids.push_back(settings.info_id_global[i]);
				}
			//}
		}

		if(local_ids.size()){
			// Dedupe vectors. This prevents multiple parsings of the same
			// target data container as this is illegal.
			for(uint32_t i = 0; i < local_ids.size(); ++i){
				std::unordered_map<int, int>::const_iterator it = this->load_settings->info_map_global.find(global_ids[i]);
				if(it == this->load_settings->info_map_global.end()){
					this->load_settings->info_id_local_loaded.push_back(local_ids[i]);
					this->load_settings->info_id_global_loaded.push_back(global_ids[i]);
					this->load_settings->info_map_global[global_ids[i]] = i;
				}
			}
		}
	}

	// Parse Format. If all Format containers are loaded then we simply copy
	// the order in which they occur. If we are provided with a vector
	// of target global identifiers we have to first map these to the
	// (possible) local identifiers.
	if(settings.load_static & YON_BLK_BV_FORMAT){
		for(uint32_t i = 0; i < this->footer.n_format_streams; ++i){
			this->load_settings->format_id_local_loaded.push_back(i);
			this->load_settings->format_id_global_loaded.push_back(this->footer.format_offsets[i].data_header.global_key);
			this->load_settings->format_map_global[this->load_settings->format_id_global_loaded[i]] = i;
		}
	} else {
		std::vector<int> local_ids;
		std::vector<int> global_ids;
		for(uint32_t i = 0; i < settings.format_id_global.size(); ++i){
			// Searches for the global Vcf:FORMAT idx value in the block. If
			// it is found then return that local idx otherwise -1. If the
			// idx is found store it in the loaded idx vector.
			const int local = this->GetFormatPosition(settings.format_id_global[i]);
			if(local >= 0){
				local_ids.push_back(local);
				global_ids.push_back(settings.format_id_global[i]);
			}
		}

		if(local_ids.size()){
			// Dedupe vectors. This prevents multiple parsings of the same
			// target data container as this is illegal.
			for(uint32_t i = 0; i < local_ids.size(); ++i){
				std::unordered_map<int, int>::const_iterator it = this->load_settings->format_map_global.find(global_ids[i]);
				if(it == this->load_settings->format_map_global.end()){
					this->load_settings->format_id_local_loaded.push_back(local_ids[i]);
					this->load_settings->format_id_global_loaded.push_back(global_ids[i]);
					this->load_settings->format_map_global[global_ids[i]] = i;
				}
			}
		}
	}

	return(this->ParseLoadedPatterns(settings));
}

bool yon1_vb_t::ParseLoadedPatterns(yon_vb_settings& settings){
	// Clear previous information (if any).
	this->load_settings->info_patterns_local.clear();
	this->load_settings->format_patterns_local.clear();
	this->load_settings->info_patterns_local.resize(this->footer.n_info_patterns);
	this->load_settings->format_patterns_local.resize(this->footer.n_format_patterns);

	// Iterate over Info patterns.
	if(this->load_settings->info_id_global_loaded.size()){
		// If all Vcf::INFO fields are desired then return them
		// in the stored order to guarantee bit-exactness. Otherwise
		// return in the order requested.
		if((settings.load_static & YON_BLK_BV_INFO) && settings.annotate_extra == false){
			for(uint32_t p = 0; p < this->footer.n_info_patterns; ++p){
				this->load_settings->info_patterns_local[p] = this->footer.info_patterns[p].pattern;
			}
		} else { // Return in requested order.
			for(uint32_t p = 0; p < this->footer.n_info_patterns; ++p){
				this->load_settings->info_patterns_local[p] = this->IntersectInfoPatterns(this->load_settings->info_id_global_loaded, p);
			}
		}
	}

	// Iterate over Format patterns.
	if(this->load_settings->format_id_global_loaded.size()){
		if(settings.load_static & YON_BLK_BV_FORMAT){
			// If all Vcf::FORMAT fields are desired then return them
			// in the stored order to guarantee bit-exactness. Otherwise
			// return in the order requested.
			for(uint32_t p = 0; p < this->footer.n_format_patterns; ++p){
				this->load_settings->format_patterns_local[p] = this->footer.format_patterns[p].pattern;
			}
		} else {
			for(uint32_t p = 0; p < this->footer.n_format_patterns; ++p){
				this->load_settings->format_patterns_local[p] = this->IntersectFormatPatterns(this->load_settings->format_id_global_loaded, p);
			}
		}
	}

	return true;
}

bool yon1_vb_t::read(std::ifstream& stream,
                     block_settings_type& settings,
                     const yon_vnt_hdr_t& header)
{
	if(this->load_settings == nullptr)
		this->load_settings = new yon_blk_load_settings;

	// Allocate enough memory to store all available Format and
	// Info containers.
	delete [] this->info_containers;
	this->info_containers    = new yon1_vb_t::container_type[this->footer.n_info_streams];
	this->n_info_c_allocated = this->footer.n_info_streams;

	delete [] this->format_containers;
	this->format_containers    = new yon1_vb_t::container_type[this->footer.n_format_streams];
	this->n_format_c_allocated = this->footer.n_format_streams;

	// Interpret the user-specified block-settings if any. This step converts
	// global index offset values into local offsets and computes new pattern
	// vectors if required. The ordering of the values are according to the
	// input sequence not according to the actual stored order.
	if(this->ParseSettings(settings, header) == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to interpret block settings..." << std::endl;
		return false;
	}

	// Load the FORMAT:GT (GBPBWT) permutation array.
	if(settings.load_static & YON_BLK_BV_PPA){
		// If there is FORMAT:GT field data available AND that data has
		// been permuted then create a new yon_gt_ppa object to store
		// this data.
		if(this->header.controller.has_gt_permuted && this->header.controller.has_gt){
			stream.seekg(this->start_compressed_data_ + this->footer.offsets[YON_BLK_PPA].data_header.offset);
			this->LoadContainerSeek(stream,
									this->footer.offsets[YON_BLK_PPA],
									this->base_containers[YON_BLK_PPA]);

			this->gt_ppa = new yon_gt_ppa;
			this->gt_ppa->n_s = header.GetNumberSamples();
		}
	}

	// Load base meta containers.
	for(uint32_t i = YON_BLK_CONTIG; i < YON_BLK_GT_INT8; ++i){
		if(settings.load_static & (1 << i)){
			this->LoadContainerSeek(stream,
									this->footer.offsets[i],
									this->base_containers[i]);
		}
	}

	// Load genotype containers. At the moment, genotype containers
	// cannot be loaded individually by using this wrapper routine.
	// If you wish to load these separately you will have to do
	// so manually.
	if((settings.load_static & YON_BLK_BV_GT) || (settings.load_static & YON_BLK_BV_FORMAT)){
		this->load_settings->loaded_genotypes = true;
		this->LoadContainerSeek(stream, this->footer.offsets[YON_BLK_GT_INT8], this->base_containers[YON_BLK_GT_INT8]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_INT16],    this->base_containers[YON_BLK_GT_INT16]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_INT32],    this->base_containers[YON_BLK_GT_INT32]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_INT64],    this->base_containers[YON_BLK_GT_INT64]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_S_INT8],   this->base_containers[YON_BLK_GT_S_INT8]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_S_INT16],  this->base_containers[YON_BLK_GT_S_INT16]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_S_INT32],  this->base_containers[YON_BLK_GT_S_INT32]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_S_INT64],  this->base_containers[YON_BLK_GT_S_INT64]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_N_INT8],   this->base_containers[YON_BLK_GT_N_INT8]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_N_INT16],  this->base_containers[YON_BLK_GT_N_INT16]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_N_INT32],  this->base_containers[YON_BLK_GT_N_INT32]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_N_INT64],  this->base_containers[YON_BLK_GT_N_INT64]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_SUPPORT],  this->base_containers[YON_BLK_GT_SUPPORT]);
		this->LoadContainer(stream, this->footer.offsets[YON_BLK_GT_PLOIDY],   this->base_containers[YON_BLK_GT_PLOIDY]);
	}

	// Load Info containers. Technically there is no difference between the two
	// conditions below in terms of outcome. However, the first case guarantees
	// that data is loaded linearly from disk as this can be guaranteed when loading
	// all available data. There is no such guarntees for the second case.
	if(this->footer.n_info_streams && (settings.load_static & YON_BLK_BV_INFO) && settings.annotate_extra == false){
		stream.seekg(this->start_compressed_data_ + this->footer.info_offsets[0].data_header.offset);

		for(uint32_t i = 0; i < this->footer.n_info_streams; ++i){
			this->LoadContainer(stream,
								this->footer.info_offsets[i],
								this->info_containers[i]);
		}
	}
	// If we have a user-supplied list of identifiers parsed above.
	else {
		for(uint32_t i = 0; i < this->load_settings->info_id_local_loaded.size(); ++i){
			this->LoadContainerSeek(stream,
									this->footer.info_offsets[this->load_settings->info_id_local_loaded[i]],
									this->info_containers[this->load_settings->info_id_local_loaded[i]]);
		}

	}

	// Load Format containers. Technically there is no difference between the two
	// conditions below in terms of outcome. However, the first case guarantees
	// that data is loaded linearly from disk as this can be guaranteed when loading
	// all available data. There is no such guarntees for the second case.
	for(uint32_t i = 0; i < this->load_settings->format_id_local_loaded.size(); ++i){
		this->LoadContainerSeek(stream,
		                        this->footer.format_offsets[this->load_settings->format_id_local_loaded[i]],
		                        this->format_containers[this->load_settings->format_id_local_loaded[i]]);
	}

	// Seek to end-of-block position.
	stream.seekg(this->end_block_);
	return(true);
}

uint64_t yon1_vb_t::GetCompressedSize(void) const{
	uint64_t total = 0;
	if(this->header.controller.has_gt && this->header.controller.has_gt_permuted)
		total += this->base_containers[YON_BLK_PPA].GetObjectSize();

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i)              total += this->base_containers[i].GetObjectSize();
	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i)   total += this->info_containers[i].GetObjectSize();
	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i) total += this->format_containers[i].GetObjectSize();

	return(total);
}

uint64_t yon1_vb_t::GetUncompressedSize(void) const{
	uint64_t total = 0;
	if(this->header.controller.has_gt && this->header.controller.has_gt_permuted)
		total += this->base_containers[YON_BLK_PPA].data_uncompressed.size();

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i)              total += this->base_containers[i].data_uncompressed.size() + this->base_containers[i].strides_uncompressed.size();
	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i)   total += this->info_containers[i].data_uncompressed.size() + this->info_containers[i].strides_uncompressed.size();
	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i) total += this->format_containers[i].data_uncompressed.size() + this->format_containers[i].strides_uncompressed.size();

	return(total);
}

void yon1_vb_t::UpdateOutputStatistics(import_stats_type& stats_basic,
                                          import_stats_type& stats_info,
                                          import_stats_type& stats_format)
{
	if(this->header.controller.has_gt && this->header.controller.has_gt_permuted)
		stats_basic[0] += this->base_containers[YON_BLK_PPA];

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i)
		stats_basic[i] += this->base_containers[i];

	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i){
		stats_info[this->footer.info_offsets[i].data_header.global_key] += this->info_containers[i];
	}

	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i){
		stats_format[this->footer.format_offsets[i].data_header.global_key] += this->format_containers[i];
	}
}

bool yon1_vb_t::write(std::ostream& stream)
{
	if(stream.good() == false){
		return false;
	}

	// Keep track of start offset and other offets in the
	// stream.
	const uint64_t begin_pos = stream.tellp();
	this->header.l_offset_footer = this->GetCompressedSize();
	stream << this->header;
	const uint64_t start_pos = stream.tellp();

	if(this->header.controller.has_gt && this->header.controller.has_gt_permuted)
		this->WriteContainer(stream, this->footer.offsets[YON_BLK_PPA], this->base_containers[YON_BLK_PPA], (uint64_t)stream.tellp() - start_pos);

	// Start at offset 1 because offset 0 (YON_BLK_PPA) is encoding for the
	// genotype permutation array that is handled differently.
	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i)
		this->WriteContainer(stream, this->footer.offsets[i], this->base_containers[i], (uint64_t)stream.tellp() - start_pos);

	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i)
		this->WriteContainer(stream, this->footer.info_offsets[i], this->info_containers[i], (uint64_t)stream.tellp() - start_pos);

	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i)
		this->WriteContainer(stream, this->footer.format_offsets[i], this->format_containers[i], (uint64_t)stream.tellp() - start_pos);

	// Assert that the written amount equals the expected amount.
	assert(this->header.l_offset_footer == (uint64_t)stream.tellp() - start_pos);

	return(stream.good());
}

bool yon1_vb_t::operator+=(yon1_vnt_t& rcd){
	// Meta positions
	this->base_containers[YON_BLK_POSITION].Add((int64_t)rcd.pos);
	++this->base_containers[YON_BLK_POSITION];

	// Contig ID
	this->base_containers[YON_BLK_CONTIG].Add((int32_t)rcd.rid);
	++this->base_containers[YON_BLK_CONTIG];

	// Ref-alt data
	if(rcd.UsePackedRefAlt()){ // Is simple SNV and possible extra case when <NON_REF> in gVCF
		rcd.controller.alleles_packed = true;
		const uint8_t ref_alt = rcd.PackRefAltByte();
		this->base_containers[YON_BLK_REFALT].AddLiteral(ref_alt);
		++this->base_containers[YON_BLK_REFALT];
	}
	// add complex
	else {
		// Special encoding
		for(uint32_t i = 0; i < rcd.n_alleles; ++i){
			// Write out allele
			this->base_containers[YON_BLK_ALLELES].AddLiteral((uint16_t)rcd.alleles[i].l_allele);
			this->base_containers[YON_BLK_ALLELES].AddCharacter(rcd.alleles[i].allele, rcd.alleles[i].l_allele);
		}
		++this->base_containers[YON_BLK_ALLELES]; // update before to not trigger
		this->base_containers[YON_BLK_ALLELES].AddStride(rcd.n_alleles);
	}

	// Quality
	this->base_containers[YON_BLK_QUALITY].Add(rcd.qual);
	++this->base_containers[YON_BLK_QUALITY];

	// Variant name
	this->base_containers[YON_BLK_NAMES].AddStride(rcd.name.size());
	this->base_containers[YON_BLK_NAMES].AddCharacter(rcd.name);
	++this->base_containers[YON_BLK_NAMES];

	// Tachyon pattern identifiers
	this->base_containers[YON_BLK_ID_INFO].Add(rcd.info_pid);
	this->base_containers[YON_BLK_ID_FORMAT].Add(rcd.fmt_pid);
	this->base_containers[YON_BLK_ID_FILTER].Add(rcd.flt_pid);
	++this->base_containers[YON_BLK_ID_INFO];
	++this->base_containers[YON_BLK_ID_FORMAT];
	++this->base_containers[YON_BLK_ID_FILTER];

	// Check if all variants are of length 1 (as in all alleles are SNVs)
	bool all_snv = true;
	for(uint32_t i = 0; i < rcd.n_alleles; ++i){
		if(rcd.alleles[i].size() != 1) all_snv = false;
	}
	rcd.controller.all_snv = all_snv;

	// Controller
	this->base_containers[YON_BLK_CONTROLLER].AddLiteral((uint16_t)rcd.controller.ToValue()); // has been overloaded
	++this->base_containers[YON_BLK_CONTROLLER];

	// Ploidy
	this->base_containers[YON_BLK_GT_PLOIDY].Add(rcd.n_base_ploidy);
	++this->base_containers[YON_BLK_GT_PLOIDY];

	return true;
}

std::vector<int> yon1_vb_t::IntersectInfoKeys(const std::vector<int>& info_ids_global) const{
	std::vector<int> info_ids_found;
	if(info_ids_global.size() == 0) return(info_ids_found);

	for(uint32_t i = 0; i < info_ids_global.size(); ++i){
		for(uint32_t j = 0; j < this->footer.n_info_streams; ++j){
			if(this->footer.info_offsets[j].data_header.global_key == info_ids_global[i])
				info_ids_found.push_back(this->footer.info_offsets[j].data_header.global_key);
		}
	}

	return(info_ids_found);
}

std::vector<int> yon1_vb_t::IntersectFormatKeys(const std::vector<int>& format_ids_global) const{
	std::vector<int> format_ids_found;
	if(format_ids_global.size() == 0) return(format_ids_found);

	for(uint32_t i = 0; i < format_ids_global.size(); ++i){
		for(uint32_t j = 0; j < this->footer.n_format_streams; ++j){
			if(this->footer.format_offsets[j].data_header.global_key == format_ids_global[i])
				format_ids_found.push_back(this->footer.format_offsets[j].data_header.global_key);
		}
	}

	return(format_ids_found);
}

std::vector<int> yon1_vb_t::IntersectFilterKeys(const std::vector<int>& filter_ids_global) const{
	std::vector<int> filter_ids_found;
	if(filter_ids_global.size() == 0) return(filter_ids_found);

	for(uint32_t i = 0; i < filter_ids_global.size(); ++i){
		for(uint32_t j = 0; j < this->footer.n_filter_streams; ++j){
			if(this->footer.filter_offsets[j].data_header.global_key == filter_ids_global[i])
				filter_ids_found.push_back(this->footer.filter_offsets[j].data_header.global_key);
		}
	}

	return(filter_ids_found);
}

std::vector<int> yon1_vb_t::IntersectInfoPatterns(const std::vector<int>& info_ids_global, const uint32_t local_id) const{
	std::vector<int> info_ids_found;
	if(info_ids_global.size() == 0) return(info_ids_found);
	assert(local_id < this->footer.n_info_patterns);

	for(uint32_t i = 0; i < info_ids_global.size(); ++i){
		for(uint32_t k = 0; k < this->footer.info_patterns[local_id].pattern.size(); ++k){
			if(this->footer.info_patterns[local_id].pattern[k] == info_ids_global[i]){
				info_ids_found.push_back(this->footer.info_patterns[local_id].pattern[k]);
			}
		}
	}

	return(info_ids_found);
}

std::vector<int> yon1_vb_t::IntersectFormatPatterns(const std::vector<int>& format_ids_global, const uint32_t local_id) const{
	std::vector<int> format_ids_found;
	if(format_ids_global.size() == 0) return(format_ids_found);
	assert(local_id < this->footer.n_format_patterns);

	for(uint32_t i = 0; i < format_ids_global.size(); ++i){
		for(uint32_t k = 0; k < this->footer.format_patterns[local_id].pattern.size(); ++k){
			if(this->footer.format_patterns[local_id].pattern[k] == format_ids_global[i])
				format_ids_found.push_back(this->footer.format_patterns[local_id].pattern[k]);
		}
	}

	return(format_ids_found);
}

std::vector<int> yon1_vb_t::IntersectFilterPatterns(const std::vector<int>& filter_ids_global, const uint32_t local_id) const{
	std::vector<int> filter_ids_found;
	if(filter_ids_global.size() == 0) return(filter_ids_found);
	assert(local_id < this->footer.n_filter_patterns);

	for(uint32_t i = 0; i < filter_ids_global.size(); ++i){
		for(uint32_t k = 0; k < this->footer.filter_patterns[local_id].pattern.size(); ++k){
			if(this->footer.filter_patterns[local_id].pattern[k] == filter_ids_global[i])
				filter_ids_found.push_back(this->footer.filter_patterns[local_id].pattern[k]);
		}
	}

	return(filter_ids_found);
}

std::vector<uint32_t> yon1_vb_t::GetInfoKeys(void) const{
	std::vector<uint32_t> ret;
	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i)
		ret.push_back(this->footer.info_offsets[i].data_header.global_key);

	return(ret);
}

std::vector<uint32_t> yon1_vb_t::GetFormatKeys(void) const{
	std::vector<uint32_t> ret;
	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i)
		ret.push_back(this->footer.format_offsets[i].data_header.global_key);

	return(ret);
}

std::vector<uint32_t> yon1_vb_t::GetFilterKeys(void) const{
	std::vector<uint32_t> ret;
	for(uint32_t i = 0; i < this->footer.n_filter_streams; ++i)
		ret.push_back(this->footer.filter_offsets[i].data_header.global_key);

	return(ret);
}

void yon1_vb_t::Finalize(void){
	this->footer.Finalize();

	// Pre-calculate the virtual file offsets prior to writing the block
	// to disk. This is required during parallel processing as we do not
	// want to the writer slave to spend time compressing or doing any
	// other compute instead of simply writing at I/O saturated speeds.
	uint64_t b_offset = 0;
	if(this->header.controller.has_gt && this->header.controller.has_gt_permuted){
		this->UpdateHeader(this->footer.offsets[YON_BLK_PPA], this->base_containers[YON_BLK_PPA], 0);
		b_offset += this->base_containers[YON_BLK_PPA].GetObjectSize();
	}

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i){
		this->UpdateHeader(this->footer.offsets[i], this->base_containers[i], b_offset);
		b_offset += this->base_containers[i].GetObjectSize();
	}

	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i){
		this->UpdateHeader(this->footer.info_offsets[i], this->info_containers[i], b_offset);
		b_offset += this->info_containers[i].GetObjectSize();
	}

	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i){
		this->UpdateHeader(this->footer.format_offsets[i], this->format_containers[i], b_offset);
		b_offset += this->format_containers[i].GetObjectSize();
	}
}

int32_t yon1_vb_t::GetInfoPosition(const uint32_t global_id) const{
	if(this->footer.info_map == nullptr) return -1;
	yon_vb_ftr::map_type::const_iterator it = this->footer.info_map->find(global_id);
	if(it == this->footer.info_map->end()) return -1;
	return(it->second);
}

int32_t yon1_vb_t::GetFormatPosition(const uint32_t global_id) const{
	if(this->footer.format_map == nullptr) return -1;
	yon_vb_ftr::map_type::const_iterator it = this->footer.format_map->find(global_id);
	if(it == this->footer.format_map->end()) return -1;
	return(it->second);
}

int32_t yon1_vb_t::GetFilterPosition(const uint32_t global_id) const{
	if(this->footer.filter_map == nullptr) return -1;
	yon_vb_ftr::map_type::const_iterator it = this->footer.filter_map->find(global_id);
	if(it == this->footer.filter_map->end()) return -1;
	return(it->second);
}

bool yon1_vb_t::HasInfo(const uint32_t global_id) const{
	if(this->footer.info_map == nullptr) return false;
	yon_vb_ftr::map_type::const_iterator it = this->footer.info_map->find(global_id);
	if(it == this->footer.info_map->end()) return false;
	return(true);
}

bool yon1_vb_t::HasFormat(const uint32_t global_id) const{
	if(this->footer.format_map == nullptr) return false;
	yon_vb_ftr::map_type::const_iterator it = this->footer.format_map->find(global_id);
	if(it == this->footer.format_map->end()) return false;
	return(true);
}

bool yon1_vb_t::HasFilter(const uint32_t global_id) const{
	if(this->footer.filter_map == nullptr) return false;
	yon_vb_ftr::map_type::const_iterator it = this->footer.filter_map->find(global_id);
	if(it == this->footer.filter_map->end()) return false;
	return(true);
}

yon1_dc_t* yon1_vb_t::GetInfoContainer(const uint32_t global_id) const{
	if(this->HasInfo(global_id))
		return(&this->info_containers[this->footer.info_map->at(global_id)]);
	else
		return nullptr;
}

yon1_dc_t* yon1_vb_t::GetFormatContainer(const uint32_t global_id) const{
	if(this->HasFormat(global_id))
		return(&this->format_containers[this->footer.format_map->at(global_id)]);
	else
		return nullptr;
}

std::vector<bool> yon1_vb_t::InfoPatternSetMembership(const int value) const{
	std::vector<bool> matches(this->footer.n_info_patterns, false);
	for(uint32_t i = 0; i < this->footer.n_info_patterns; ++i){
		for(uint32_t j = 0; j < this->footer.info_patterns[i].pattern.size(); ++j){
			if(this->footer.info_patterns[i].pattern[j] == value){
				matches[i] = true;
				break;
			}
		}
	}
	return(matches);
}

std::vector<bool> yon1_vb_t::FormatPatternSetMembership(const int value) const{
	std::vector<bool> matches(this->footer.n_format_patterns, false);
	for(uint32_t i = 0; i < this->footer.n_format_patterns; ++i){
		for(uint32_t j = 0; j < this->footer.format_patterns[i].pattern.size(); ++j){
			if(this->footer.format_patterns[i].pattern[j] == value){
				matches[i] = true;
				break;
			}
		}
	}
	return(matches);
}

std::vector<bool> yon1_vb_t::FilterPatternSetMembership(const int value) const{
	std::vector<bool> matches(this->footer.n_filter_patterns, false);
	for(uint32_t i = 0; i < this->footer.n_filter_patterns; ++i){
		for(uint32_t j = 0; j < this->footer.filter_patterns[i].pattern.size(); ++j){
			if(this->footer.filter_patterns[i].pattern[j] == value){
				matches[i] = true;
				break;
			}
		}
	}
	return(matches);
}

}
