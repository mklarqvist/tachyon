#include "variant_block.h"
#include "algorithm/compression/compression_container.h"
#include "support/helpers.h"

namespace tachyon{
namespace containers{

VariantBlock::VariantBlock() :
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

VariantBlock::VariantBlock(const uint16_t n_info, const uint16_t n_format) :
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

VariantBlock::~VariantBlock(){
	delete [] this->base_containers;
	delete [] this->info_containers;
	delete [] this->format_containers;
	delete this->gt_ppa;
	delete this->load_settings;
}

VariantBlock::VariantBlock(const self_type& other) :
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

VariantBlock::VariantBlock(self_type&& other) noexcept :
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

VariantBlock& VariantBlock::operator=(const self_type& other){
	delete [] this->base_containers;
	delete [] this->info_containers;
	delete [] this->format_containers;
	delete this->gt_ppa;
	delete this->load_settings;
	*this = VariantBlock(other); // invoke copy ctor
	return(*this);
}

VariantBlock& VariantBlock::operator=(self_type&& other) noexcept{
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

void VariantBlock::Allocate(const uint16_t n_info,
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

void VariantBlock::clear(void){
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

void VariantBlock::resize(const uint32_t s){
	if(s == 0) return;

	for(uint32_t i = 0; i < YON_BLK_N_STATIC; ++i) this->base_containers[i].resize(s);
	for(uint32_t i = 0; i < n_info_c_allocated; ++i) this->info_containers[i].resize(s);
	for(uint32_t i = 0; i < n_format_c_allocated; ++i) this->format_containers[i].resize(s);
}

void VariantBlock::UpdateContainers(void){
	this->base_containers[YON_BLK_CONTIG].UpdateContainer();
	this->base_containers[YON_BLK_POSITION].UpdateContainer();
	this->base_containers[YON_BLK_REFALT].UpdateContainer();
	this->base_containers[YON_BLK_QUALITY].UpdateContainer();
	this->base_containers[YON_BLK_NAMES].UpdateContainer();
	this->base_containers[YON_BLK_ALLELES].UpdateContainer(false, false);
	this->base_containers[YON_BLK_ID_FILTER].UpdateContainer();
	this->base_containers[YON_BLK_ID_FORMAT].UpdateContainer();
	this->base_containers[YON_BLK_ID_INFO].UpdateContainer();
	this->base_containers[YON_BLK_GT_SUPPORT].UpdateContainer();
	this->base_containers[YON_BLK_GT_PLOIDY].UpdateContainer();
	this->base_containers[YON_BLK_CONTROLLER].UpdateContainer(false, false);
	this->base_containers[YON_BLK_GT_INT8].UpdateContainer();
	this->base_containers[YON_BLK_GT_INT16].UpdateContainer();
	this->base_containers[YON_BLK_GT_INT32].UpdateContainer();
	this->base_containers[YON_BLK_GT_INT64].UpdateContainer();
	this->base_containers[YON_BLK_GT_S_INT8].UpdateContainer();
	this->base_containers[YON_BLK_GT_S_INT16].UpdateContainer();
	this->base_containers[YON_BLK_GT_S_INT32].UpdateContainer();
	this->base_containers[YON_BLK_GT_S_INT64].UpdateContainer();
	this->base_containers[YON_BLK_GT_N_INT8].UpdateContainer(false, false);
	this->base_containers[YON_BLK_GT_N_INT16].UpdateContainer(false, false);
	this->base_containers[YON_BLK_GT_N_INT32].UpdateContainer(false, false);
	this->base_containers[YON_BLK_GT_N_INT64].UpdateContainer(false, false);

	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i){
		assert(this->info_containers[i].header.data_header.stride != 0);
		this->info_containers[i].UpdateContainer();
	}

	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i){
		assert(this->format_containers[i].header.data_header.stride != 0);
		this->format_containers[i].UpdateContainer();
	}
}

bool VariantBlock::ReadHeaderFooter(std::ifstream& stream){
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
	assert(eof_marker == constants::TACHYON_BLOCK_EOF);
	this->end_block_ = stream.tellg(); // end-of-block offset
	stream.seekg(this->start_compressed_data_);
	return(stream.good());
}

bool VariantBlock::read(std::ifstream& stream){
	if(this->header.controller.hasGTPermuted && this->header.controller.hasGT){
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


bool VariantBlock::ParseSettings(DataBlockSettings& settings, const VariantHeader& header){
	// Clear previous information (if any).
	this->load_settings->clear();

	// Construct a black-list of INFO fields that should not be loaded
	// or displayed as they are being re-calculated internally and
	// emitted. If a blacklisted field is available in this block then
	// add that tag to the map.
	std::unordered_map<uint32_t, std::string> blocked_list;
	if(settings.annotate_extra){
		for(uint32_t i = 0; i < YON_GT_ANNOTATE_FIELDS.size(); ++i){
			const YonInfo* info = header.GetInfo(YON_GT_ANNOTATE_FIELDS[i]);
			if(info != nullptr){
				blocked_list[info->idx] = YON_GT_ANNOTATE_FIELDS[i];
			}
		}
	}

	// Parse Info. If all Info containers are loaded then we simply copy
	// the order in which they occur. If we are provided with a vector
	// of target global identifiers we have to first map these to the
	// (possible) local identifiers.
	if(settings.load_static & YON_BLK_BV_INFO){
		for(uint32_t i = 0; i < this->footer.n_info_streams; ++i){
			const std::unordered_map<uint32_t, std::string>::const_iterator it = blocked_list.find(this->footer.info_offsets[i].data_header.global_key);
			if(it == blocked_list.end()){
				//std::cerr << "adding not blocked" << std::endl;
				this->load_settings->info_id_local_loaded.push_back(i);
				this->load_settings->info_id_global_loaded.push_back(this->footer.info_offsets[i].data_header.global_key);
				this->load_settings->info_map_global[this->load_settings->info_id_global_loaded[i]] = i;
			} else {
				//std::cerr << "skipping blocked" << std::endl;
			}
		}
		//std::cerr << this->info_id_local_loaded.size() << "," << this->info_id_global_loaded.size() << std::endl;
	} else {
		std::vector<int> local_ids;
		std::vector<int> global_ids;
		for(uint32_t i = 0; i < settings.info_id_global.size(); ++i){
			// Searches for the global Vcf:INFO idx value in the block. If
			// it is found then return that local idx otherwise -1. If the
			// idx is found store it in the loaded idx vector.
			const std::unordered_map<uint32_t, std::string>::const_iterator it = blocked_list.find(this->footer.info_offsets[i].data_header.global_key);
			if(it == blocked_list.end()){
				const int local = this->GetInfoPosition(settings.info_id_global[i]);
				if(local >= 0){
					local_ids.push_back(local);
					global_ids.push_back(settings.info_id_global[i]);
				}
			}
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

bool VariantBlock::ParseLoadedPatterns(DataBlockSettings& settings){
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

bool VariantBlock::read(std::ifstream& stream,
                        block_settings_type& settings,
                        const VariantHeader& header)
{
	if(this->load_settings == nullptr)
		this->load_settings = new yon_blk_load_settings;

	// Allocate enough memory to store all available Format and
	// Info containers.
	delete [] this->info_containers;
	this->info_containers    = new VariantBlock::container_type[this->footer.n_info_streams];
	this->n_info_c_allocated = this->footer.n_info_streams;

	delete [] this->format_containers;
	this->format_containers    = new VariantBlock::container_type[this->footer.n_format_streams];
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
		if(this->header.controller.hasGTPermuted && this->header.controller.hasGT){
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

uint64_t VariantBlock::GetCompressedSize(void) const{
	uint64_t total = 0;
	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted)
		total += this->base_containers[YON_BLK_PPA].GetObjectSize();

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i)              total += this->base_containers[i].GetObjectSize();
	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i)   total += this->info_containers[i].GetObjectSize();
	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i) total += this->format_containers[i].GetObjectSize();

	return(total);
}

void VariantBlock::UpdateOutputStatistics(import_stats_type& stats_basic,
                                          import_stats_type& stats_info,
                                          import_stats_type& stats_format)
{
	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted)
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

bool VariantBlock::write(std::ostream& stream)
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

	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted)
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

bool VariantBlock::operator+=(meta_entry_type& meta_entry){
	// Meta positions
	this->base_containers[YON_BLK_POSITION].Add((int32_t)meta_entry.position);
	++this->base_containers[YON_BLK_POSITION];

	// Contig ID
	this->base_containers[YON_BLK_CONTIG].Add((int32_t)meta_entry.contigID);
	++this->base_containers[YON_BLK_CONTIG];

	// Ref-alt data
	if(meta_entry.UsePackedRefAlt()){ // Is simple SNV and possible extra case when <NON_REF> in gVCF
		meta_entry.controller.alleles_packed = true;
		const uint8_t ref_alt = meta_entry.PackRefAltByte();
		this->base_containers[YON_BLK_REFALT].AddLiteral(ref_alt);
		++this->base_containers[YON_BLK_REFALT];
	}
	// add complex
	else {
		// Special encoding
		for(uint32_t i = 0; i < meta_entry.n_alleles; ++i){
			// Write out allele
			this->base_containers[YON_BLK_ALLELES].AddLiteral((uint16_t)meta_entry.alleles[i].l_allele);
			this->base_containers[YON_BLK_ALLELES].AddCharacter(meta_entry.alleles[i].allele, meta_entry.alleles[i].l_allele);
		}
		++this->base_containers[YON_BLK_ALLELES]; // update before to not trigger
		this->base_containers[YON_BLK_ALLELES].AddStride(meta_entry.n_alleles);
	}

	// Quality
	this->base_containers[YON_BLK_QUALITY].Add(meta_entry.quality);
	++this->base_containers[YON_BLK_QUALITY];

	// Variant name
	this->base_containers[YON_BLK_NAMES].AddStride(meta_entry.name.size());
	this->base_containers[YON_BLK_NAMES].AddCharacter(meta_entry.name);
	++this->base_containers[YON_BLK_NAMES];

	// Tachyon pattern identifiers
	this->base_containers[YON_BLK_ID_INFO].Add(meta_entry.info_pattern_id);
	this->base_containers[YON_BLK_ID_FORMAT].Add(meta_entry.format_pattern_id);
	this->base_containers[YON_BLK_ID_FILTER].Add(meta_entry.filter_pattern_id);
	++this->base_containers[YON_BLK_ID_INFO];
	++this->base_containers[YON_BLK_ID_FORMAT];
	++this->base_containers[YON_BLK_ID_FILTER];

	// Check if all variants are of length 1 (as in all alleles are SNVs)
	bool all_snv = true;
	for(uint32_t i = 0; i < meta_entry.n_alleles; ++i){
		if(meta_entry.alleles[i].size() != 1) all_snv = false;
	}
	meta_entry.controller.all_snv = all_snv;

	// Controller
	this->base_containers[YON_BLK_CONTROLLER].AddLiteral((uint16_t)meta_entry.controller.toValue()); // has been overloaded
	++this->base_containers[YON_BLK_CONTROLLER];

	// Ploidy
	this->base_containers[YON_BLK_GT_PLOIDY].Add(meta_entry.n_base_ploidy);
	++this->base_containers[YON_BLK_GT_PLOIDY];

	return true;
}

std::vector<int> VariantBlock::IntersectInfoKeys(const std::vector<int>& info_ids_global) const{
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

std::vector<int> VariantBlock::IntersectFormatKeys(const std::vector<int>& format_ids_global) const{
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

std::vector<int> VariantBlock::IntersectFilterKeys(const std::vector<int>& filter_ids_global) const{
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

std::vector<int> VariantBlock::IntersectInfoPatterns(const std::vector<int>& info_ids_global, const uint32_t local_id) const{
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

std::vector<int> VariantBlock::IntersectFormatPatterns(const std::vector<int>& format_ids_global, const uint32_t local_id) const{
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

std::vector<int> VariantBlock::IntersectFilterPatterns(const std::vector<int>& filter_ids_global, const uint32_t local_id) const{
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

std::vector<uint32_t> VariantBlock::GetInfoKeys(void) const{
	std::vector<uint32_t> ret;
	for(uint32_t i = 0; i < this->footer.n_info_streams; ++i)
		ret.push_back(this->footer.info_offsets[i].data_header.global_key);

	return(ret);
}

std::vector<uint32_t> VariantBlock::GetFormatKeys(void) const{
	std::vector<uint32_t> ret;
	for(uint32_t i = 0; i < this->footer.n_format_streams; ++i)
		ret.push_back(this->footer.format_offsets[i].data_header.global_key);

	return(ret);
}

std::vector<uint32_t> VariantBlock::GetFilterKeys(void) const{
	std::vector<uint32_t> ret;
	for(uint32_t i = 0; i < this->footer.n_filter_streams; ++i)
		ret.push_back(this->footer.filter_offsets[i].data_header.global_key);

	return(ret);
}

int32_t VariantBlock::GetInfoPosition(const uint32_t global_id) const{
	if(this->footer.info_map == nullptr) return false;
	VariantBlockFooter::map_type::const_iterator it = this->footer.info_map->find(global_id);
	if(it == this->footer.info_map->end()) return -1;
	return(it->second);
}

int32_t VariantBlock::GetFormatPosition(const uint32_t global_id) const{
	if(this->footer.format_map == nullptr) return false;
	VariantBlockFooter::map_type::const_iterator it = this->footer.format_map->find(global_id);
	if(it == this->footer.format_map->end()) return -1;
	return(it->second);
}

int32_t VariantBlock::GetFilterPosition(const uint32_t global_id) const{
	if(this->footer.filter_map == nullptr) return false;
	VariantBlockFooter::map_type::const_iterator it = this->footer.filter_map->find(global_id);
	if(it == this->footer.filter_map->end()) return -1;
	return(it->second);
}

bool VariantBlock::HasInfo(const uint32_t global_id) const{
	if(this->footer.info_map == nullptr) return false;
	VariantBlockFooter::map_type::const_iterator it = this->footer.info_map->find(global_id);
	if(it == this->footer.info_map->end()) return false;
	return(true);
}

bool VariantBlock::HasFormat(const uint32_t global_id) const{
	if(this->footer.format_map == nullptr) return false;
	VariantBlockFooter::map_type::const_iterator it = this->footer.format_map->find(global_id);
	if(it == this->footer.format_map->end()) return false;
	return(true);
}

bool VariantBlock::HasFilter(const uint32_t global_id) const{
	if(this->footer.filter_map == nullptr) return false;
	VariantBlockFooter::map_type::const_iterator it = this->footer.filter_map->find(global_id);
	if(it == this->footer.filter_map->end()) return false;
	return(true);
}

DataContainer* VariantBlock::GetInfoContainer(const uint32_t global_id) const{
	if(this->HasInfo(global_id))
		return(&this->info_containers[this->footer.info_map->at(global_id)]);
	else
		return nullptr;
}

DataContainer* VariantBlock::GetFormatContainer(const uint32_t global_id) const{
	if(this->HasFormat(global_id))
		return(&this->format_containers[this->footer.format_map->at(global_id)]);
	else
		return nullptr;
}

std::vector<bool> VariantBlock::InfoPatternSetMembership(const int value) const{
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

std::vector<bool> VariantBlock::FormatPatternSetMembership(const int value) const{
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

std::vector<bool> VariantBlock::FilterPatternSetMembership(const int value) const{
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
}
