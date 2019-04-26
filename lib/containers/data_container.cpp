#include "data_container.h"
#include "algorithm/digest/variant_digest_manager.h"
#include "third_party/xxhash/xxhash.h"

namespace tachyon {

yon_blk_bv_pair::yon_blk_bv_pair() : l_bytes(0), bit_bytes(nullptr) {}
yon_blk_bv_pair::~yon_blk_bv_pair() { delete [] this->bit_bytes; }

void yon_blk_bv_pair::clear(void) {
	this->pattern.clear();
	this->l_bytes = 0;
	delete [] this->bit_bytes;
	this->bit_bytes = nullptr;
}

yon_blk_bv_pair& yon_blk_bv_pair::operator=(const yon_blk_bv_pair& other) {
	delete [] this->bit_bytes;
	this->pattern   = other.pattern;
	this->l_bytes   = other.l_bytes;
	this->bit_bytes = new uint8_t[this->l_bytes];
	memcpy(this->bit_bytes, other.bit_bytes, this->l_bytes);
	return(*this);
}

yon_blk_bv_pair& yon_blk_bv_pair::operator=(yon_blk_bv_pair&& other) noexcept{
	if (this == &other) {
		// take precautions against self-moves
		return *this;
	}

	delete [] this->bit_bytes; this->bit_bytes = nullptr;
	std::swap(this->bit_bytes, other.bit_bytes);
	this->pattern = std::move(other.pattern);
	other.pattern.clear(); // Clear the src pattern vector.
	this->l_bytes = other.l_bytes;
	other.l_bytes = 0; // Clear the src byte length.
	return(*this);
}

void yon_blk_bv_pair::Build(const uint32_t n_footer_total_fields,
                            const std::unordered_map<uint32_t, uint32_t>* local_map)
{
	if (this->pattern.size() == 0) return;
	assert(local_map != nullptr);

	// Determine the required byte width of the bit-vector.
	uint8_t bitvector_width = ceil((float)(n_footer_total_fields + 1) / 8);

	// Allocate new bit-vectors.
	delete [] this->bit_bytes;
	this->l_bytes = bitvector_width;
	this->bit_bytes = new uint8_t[bitvector_width];
	memset(this->bit_bytes, 0, sizeof(uint8_t)*bitvector_width);

	// Cycle over global idx values in the vector.
	for (uint32_t i = 0; i < this->pattern.size(); ++i) {
		std::unordered_map<uint32_t, uint32_t>::const_iterator it = local_map->find(this->pattern[i]);
		assert(it != local_map->end());

		// Map from absolute key to local key.
		uint32_t local_key = it->second;
		assert(local_key <= n_footer_total_fields);

		// Set the target bit to TRUE at the local key position.
		this->bit_bytes[local_key/8] |= 1 << (local_key % 8);
	}
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const yon_blk_bv_pair& entry) {
	SerializePrimitive(entry.l_bytes, buffer);
	buffer += (uint32_t)entry.pattern.size();
	for (uint32_t i = 0; i < entry.pattern.size(); ++i)
		SerializePrimitive(entry.pattern[i], buffer);

	for (uint32_t i = 0; i < entry.l_bytes; ++i)
		SerializePrimitive(entry.bit_bytes[i], buffer);


	return(buffer);
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, yon_blk_bv_pair& entry) {
	entry.pattern.clear();
	DeserializePrimitive(entry.l_bytes, buffer);
	uint32_t l_vector;
	buffer >> l_vector;
	//entry.pattern.resize(l_vector);
	for (uint32_t i = 0; i < l_vector; ++i) {
		int temp;
		DeserializePrimitive(temp, buffer);
		//entry.pattern[i] = temp;
		entry.pattern.push_back(temp);
	}

	entry.bit_bytes = new uint8_t[entry.l_bytes];
	for (uint32_t i = 0; i < entry.l_bytes; ++i)
		DeserializePrimitive(entry.bit_bytes[i], buffer);

	return(buffer);
}

// Header
yon_dc_hdr_cont::yon_dc_hdr_cont() :
	signedness(0), mixedStride(0), type(0), encoder(0),
	uniform(0), encryption(0), preprocessor(0)
{}

yon_dc_hdr_cont::~yon_dc_hdr_cont() {}

inline void yon_dc_hdr_cont::clear() {
	this->signedness   = 0;
	this->mixedStride  = 0;
	this->type         = 0;
	this->encoder      = 0;
	this->uniform      = 0;
	this->encryption   = 0;
	this->preprocessor = 0;
}

yon_dc_hdr_cont& yon_dc_hdr_cont::operator=(const self_type& other) {
	this->signedness   = other.signedness;
	this->mixedStride  = other.mixedStride;
	this->type         = other.type;
	this->encoder      = other.encoder;
	this->uniform      = other.uniform;
	this->encryption   = other.encryption;
	this->preprocessor = other.preprocessor;
	return(*this);
}

bool yon_dc_hdr_cont::operator==(const self_type& other) const{
	if (this->signedness   != other.signedness)   return false;
	if (this->mixedStride  != other.mixedStride)  return false;
	if (this->type         != other.type)         return false;
	if (this->encoder      != other.encoder)      return false;
	if (this->uniform      != other.uniform)      return false;
	if (this->encryption   != other.encryption)   return false;
	if (this->preprocessor != other.preprocessor) return false;
	return true;
}

yon_buffer_t& operator<<(yon_buffer_t& buffer,const yon_dc_hdr_cont& controller) {
	const uint32_t c =
				controller.signedness   << 0  |
				controller.mixedStride  << 1  |
				controller.type         << 2  |
				controller.encoder      << 8  |
				controller.uniform      << 13 |
				controller.encryption   << 14 |
				controller.preprocessor << 16;

	//const uint16_t* c = reinterpret_cast<const uint16_t* const>(&controller);
	buffer += c;
	return(buffer);
}

std::ostream& operator<<(std::ostream& stream, const yon_dc_hdr_cont& controller) {
	const uint32_t c =
				controller.signedness   << 0  |
				controller.mixedStride  << 1  |
				controller.type         << 2  |
				controller.encoder      << 8  |
				controller.uniform      << 13 |
				controller.encryption   << 14 |
				controller.preprocessor << 16;

	//assert(*reinterpret_cast<const uint16_t* const>(&controller) == c);

	stream.write(reinterpret_cast<const char*>(&c), sizeof(uint32_t));
	return(stream);
}

std::istream& operator>>(std::istream& stream, yon_dc_hdr_cont& controller) {
	stream.read(reinterpret_cast<char*>(&controller), sizeof(uint32_t));
	return(stream);
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, yon_dc_hdr_cont& controller) {
	uint32_t* c = reinterpret_cast<uint32_t*>(&controller);
	buffer >> *c;
	return(buffer);
}

yon_dc_hdr_obj::yon_dc_hdr_obj() :
	stride(1),
	offset(0),
	cLength(0),
	uLength(0),
	eLength(0),
	global_key(-1)
{
	memset(&this->crc[0], 0, MD5_DIGEST_LENGTH);
}

yon_dc_hdr_obj::yon_dc_hdr_obj(const yon_dc_hdr_obj& other) :
	controller(other.controller),
	stride(other.stride),
	offset(other.offset),
	cLength(other.cLength),
	uLength(other.uLength),
	eLength(other.eLength),
	global_key(other.global_key)
{
	memcpy(&this->crc[0], &other.crc[0], MD5_DIGEST_LENGTH);
}

yon_dc_hdr_obj::yon_dc_hdr_obj(yon_dc_hdr_obj&& other) noexcept :
	controller(other.controller),
	stride(other.stride),
	offset(other.offset),
	cLength(other.cLength),
	uLength(other.uLength),
	eLength(other.eLength),
	global_key(other.global_key)
{
	memcpy(&this->crc[0], &other.crc[0], MD5_DIGEST_LENGTH);
}

yon_dc_hdr_obj& yon_dc_hdr_obj::operator=(const yon_dc_hdr_obj& other) {
	this->controller = other.controller;
	this->stride     = other.stride;
	this->offset     = other.offset;
	this->cLength    = other.cLength;
	this->uLength    = other.uLength;
	this->eLength    = other.eLength;
	memcpy(&this->crc[0], &other.crc[0], MD5_DIGEST_LENGTH);
	this->global_key = other.global_key;
	return *this;
}

yon_dc_hdr_obj& yon_dc_hdr_obj::operator=(yon_dc_hdr_obj&& other) noexcept{
	this->controller = other.controller;
	this->stride     = other.stride;
	this->offset     = other.offset;
	this->cLength    = other.cLength;
	this->uLength    = other.uLength;
	this->eLength    = other.eLength;
	memcpy(&this->crc[0], &other.crc[0], MD5_DIGEST_LENGTH);
	this->global_key = other.global_key;
	return *this;
}

yon_dc_hdr_obj::~yon_dc_hdr_obj() { }

void yon_dc_hdr_obj::reset(void) {
	this->controller.clear();
	this->stride     = 1;
	this->offset     = 0;
	this->cLength    = 0;
	this->uLength    = 0;
	memset(&this->crc[0], 0, MD5_DIGEST_LENGTH);
	this->global_key = -1;
}

bool yon_dc_hdr_obj::operator==(const self_type& other) const{
	if (this->stride     != other.stride)     return false;
	if (this->offset     != other.offset)     return false;
	if (this->cLength    != other.cLength)    return false;
	if (this->uLength    != other.uLength)    return false;
	if (this->eLength    != other.eLength)    return false;
	if (this->global_key != other.global_key) return false;
	if (this->controller != other.controller) return false;
	for (uint32_t i = 0; i < MD5_DIGEST_LENGTH; ++i)
		if (this->crc[i] != other.crc[i]) return false;

	return true;
}

int8_t yon_dc_hdr_obj::GetPrimitiveWidth(void) const{
	// We do not care about signedness here
	switch(this->controller.type) {
	case(YON_TYPE_UNKNOWN):
	case(YON_TYPE_STRUCT): return(-1);
	case(YON_TYPE_BOOLEAN):
	case(YON_TYPE_CHAR):   return(sizeof(char));
	case(YON_TYPE_8B):     return(sizeof(uint8_t));
	case(YON_TYPE_16B):    return(sizeof(uint16_t));
	case(YON_TYPE_32B):    return(sizeof(uint32_t));
	case(YON_TYPE_64B):    return(sizeof(uint64_t));
	case(YON_TYPE_FLOAT):  return(sizeof(float));
	case(YON_TYPE_DOUBLE): return(sizeof(double));
	}
	return 0;
}

bool yon_dc_hdr_obj::CheckChecksum(const uint8_t* compare) const{
	for (uint32_t i = 0; i < MD5_DIGEST_LENGTH; ++i) {
		if (compare[i] != this->crc[i])
			return false;
	}
	return true;
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const yon_dc_hdr_obj& entry) {
	buffer << entry.controller;
	buffer += entry.stride;
	buffer += entry.offset;
	buffer += entry.cLength;
	buffer += entry.uLength;
	buffer += entry.eLength;
	for (uint32_t i = 0; i < MD5_DIGEST_LENGTH; ++i) buffer += entry.crc[i];
	buffer += entry.global_key;
	return(buffer);
}

std::ostream& operator<<(std::ostream& stream, const yon_dc_hdr_obj& entry) {
	stream << entry.controller;
	stream.write(reinterpret_cast<const char*>(&entry.stride),    sizeof(int32_t));
	stream.write(reinterpret_cast<const char*>(&entry.offset),    sizeof(uint32_t));
	stream.write(reinterpret_cast<const char*>(&entry.cLength),   sizeof(uint32_t));
	stream.write(reinterpret_cast<const char*>(&entry.uLength),   sizeof(uint32_t));
	stream.write(reinterpret_cast<const char*>(&entry.eLength),   sizeof(uint32_t));
	stream.write(reinterpret_cast<const char*>(&entry.crc[0]),    sizeof(uint8_t)*MD5_DIGEST_LENGTH);
	stream.write(reinterpret_cast<const char*>(&entry.global_key),sizeof(int32_t));
	return(stream);
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, yon_dc_hdr_obj& entry) {
	buffer >> entry.controller;
	buffer >> entry.stride;
	buffer >> entry.offset;
	buffer >> entry.cLength;
	buffer >> entry.uLength;
	buffer >> entry.eLength;
	for (uint32_t i = 0; i < MD5_DIGEST_LENGTH; ++i) buffer >> entry.crc[i];
	buffer >> entry.global_key;
	return(buffer);
}

std::ifstream& operator>>(std::ifstream& stream, yon_dc_hdr_obj& entry) {
	stream >> entry.controller;
	stream.read(reinterpret_cast<char*>(&entry.stride),     sizeof(int32_t));
	stream.read(reinterpret_cast<char*>(&entry.offset),     sizeof(uint32_t));
	stream.read(reinterpret_cast<char*>(&entry.cLength),    sizeof(uint32_t));
	stream.read(reinterpret_cast<char*>(&entry.uLength),    sizeof(uint32_t));
	stream.read(reinterpret_cast<char*>(&entry.eLength),    sizeof(uint32_t));
	stream.read(reinterpret_cast<char*>(&entry.crc[0]),     sizeof(uint8_t)*MD5_DIGEST_LENGTH);
	stream.read(reinterpret_cast<char*>(&entry.global_key), sizeof(int32_t));

	return(stream);
}

// header
yon_dc_hdr::yon_dc_hdr() :
	identifier(0),
	n_entries(0),
	n_additions(0),
	n_strides(0)
{
}

yon_dc_hdr::yon_dc_hdr(const self_type& other) :
	identifier(other.identifier),
	n_entries(other.n_entries),
	n_additions(other.n_additions),
	n_strides(other.n_strides),
	data_header(other.data_header),
	stride_header(other.stride_header)
{

}

yon_dc_hdr::~yon_dc_hdr() {}

void yon_dc_hdr::reset(void) {
	this->identifier  = 0;
	this->n_entries   = 0;
	this->n_additions = 0;
	this->n_strides   = 0;
	this->data_header.reset();
	this->stride_header.reset();
}

yon_dc_hdr& yon_dc_hdr::operator=(const self_type& other) {
	this->identifier    = other.identifier;
	this->n_entries     = other.n_entries;
	this->n_additions   = other.n_additions;
	this->n_strides     = other.n_strides;
	this->data_header   = other.data_header;
	this->stride_header = other.stride_header;
	return(*this);
}

yon_dc_hdr& yon_dc_hdr::operator=(self_type&& other) noexcept{
	this->identifier    = other.identifier;
	this->n_entries     = other.n_entries;
	this->n_additions   = other.n_additions;
	this->n_strides     = other.n_strides;
	this->data_header   = std::move(other.data_header);
	this->stride_header = std::move(other.stride_header);
	return(*this);
}

// Comparators
bool yon_dc_hdr::operator==(const self_type& other) const {
	if (this->identifier    != other.identifier)    return false;
	if (this->n_entries     != other.n_entries)     return false;
	if (this->n_additions   != other.n_additions)   return false;
	if (this->n_strides     != other.n_strides)     return false;
	if (this->data_header   != other.data_header)   return false;
	if (this->stride_header != other.stride_header) return false;
	return true;
}

yon_dc_hdr& yon_dc_hdr::operator+=(const self_type& other) {
	this->n_entries     += other.n_entries;
	this->n_additions   += other.n_additions;
	this->n_strides     += other.n_strides;
	if (data_header.stride != other.data_header.stride || other.data_header.controller.mixedStride) {
		//if (data_header.HasMixedStride() == false)
		//	std::cerr << "triggering mixed stride: " << data_header.stride << "!=" << other.data_header.stride << std::endl;
		this->data_header.SetMixedStride(true);
	}
	return(*this);
}

yon_buffer_t& operator<<(yon_buffer_t& buffer, const yon_dc_hdr& entry) {
	buffer += entry.identifier;
	buffer += entry.n_entries;
	buffer += entry.n_additions;
	buffer += entry.n_strides;
	buffer << entry.data_header;

	if (entry.data_header.HasMixedStride())
		buffer << entry.stride_header;

	return(buffer);
}

yon_buffer_t& operator>>(yon_buffer_t& buffer, yon_dc_hdr& entry) {
	buffer >> entry.identifier;
	buffer >> entry.n_entries;
	buffer >> entry.n_additions;
	buffer >> entry.n_strides;
	buffer >> entry.data_header;

	if (entry.data_header.HasMixedStride())
		buffer >> entry.stride_header;

	return(buffer);
}

std::ostream& operator<<(std::ostream& stream, const yon_dc_hdr& entry) {
	stream.write(reinterpret_cast<const char*>(&entry.identifier),  sizeof(uint64_t));
	stream.write(reinterpret_cast<const char*>(&entry.n_entries),   sizeof(uint32_t));
	stream.write(reinterpret_cast<const char*>(&entry.n_additions), sizeof(uint32_t));
	stream.write(reinterpret_cast<const char*>(&entry.n_strides),   sizeof(uint32_t));
	stream << entry.data_header;
	if (entry.data_header.HasMixedStride())
		stream << entry.stride_header;

	return(stream);
}

std::ifstream& operator>>(std::ifstream& stream, yon_dc_hdr& entry) {
	stream.read(reinterpret_cast<char*>(&entry.identifier),  sizeof(uint64_t));
	stream.read(reinterpret_cast<char*>(&entry.n_entries),   sizeof(uint32_t));
	stream.read(reinterpret_cast<char*>(&entry.n_additions), sizeof(uint32_t));
	stream.read(reinterpret_cast<char*>(&entry.n_strides),   sizeof(uint32_t));
	stream >> entry.data_header;
	if (entry.data_header.HasMixedStride())
		stream >> entry.stride_header;

	return(stream);
}

// data container
yon1_dc_t::yon1_dc_t()
{}

yon1_dc_t::yon1_dc_t(const uint32_t start_size) :
	data(start_size),
	strides(start_size),
	data_uncompressed(start_size),
	strides_uncompressed(start_size)
{}

yon1_dc_t::~yon1_dc_t() { }

yon1_dc_t::yon1_dc_t(self_type&& other) noexcept :
	header(std::move(other.header)),
	data(std::move(other.data)),
	strides(std::move(other.strides)),
	data_uncompressed(std::move(other.data_uncompressed)),
	strides_uncompressed(std::move(other.strides_uncompressed))
{

}

yon1_dc_t::yon1_dc_t(const self_type& other) :
	header(other.header),
	data(other.data),
	strides(other.strides),
	data_uncompressed(other.data_uncompressed),
	strides_uncompressed(other.strides_uncompressed)
{}

yon1_dc_t& yon1_dc_t::operator=(const self_type& other) {
	this->data = other.data;
	this->data_uncompressed = other.data_uncompressed;
	this->strides = other.strides;
	this->strides_uncompressed = other.strides_uncompressed;
	this->header = other.header;
	return(*this);
}

yon1_dc_t& yon1_dc_t::operator=(self_type&& other) noexcept{
	this->data = std::move(other.data);
	this->data_uncompressed = std::move(other.data_uncompressed);
	this->strides = std::move(other.strides);
	this->strides_uncompressed = std::move(other.strides_uncompressed);
	this->header = std::move(other.header);
	return(*this);
}

void yon1_dc_t::reset(void) {
	this->data.reset();
	this->data_uncompressed.reset();
	this->strides.reset();
	this->strides_uncompressed.reset();
	this->header.reset();
}

void yon1_dc_t::resize(const uint32_t size) {
	this->data.resize(size);
	this->data_uncompressed.resize(size);
	this->strides.resize(size);
	this->strides_uncompressed.resize(size);
}

void yon1_dc_t::GenerateMd5(void) {
	algorithm::VariantDigestManager::GenerateMd5(this->data_uncompressed.data(), this->data_uncompressed.size(), &this->header.data_header.crc[0]);
	algorithm::VariantDigestManager::GenerateMd5(this->strides_uncompressed.data(), this->strides_uncompressed.size(), &this->header.stride_header.crc[0]);
}

bool yon1_dc_t::CheckMd5(int target) {
	if (target == 0) {
		if (this->data_uncompressed.size() == 0)
			return true;

		// Checksum for main buffer
		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->data_uncompressed.data(), this->data_uncompressed.size(), &md5_compare[0]);
		return(this->header.data_header.CheckChecksum(md5_compare));

	} else if (target == 1) {
		if (this->strides_uncompressed.size() == 0)
			return true;

		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->strides_uncompressed.data(), this->strides_uncompressed.size(), &md5_compare[0]);
		return(this->header.stride_header.CheckChecksum(md5_compare));

	} else if (target == 3) {
		if (this->data.size() == 0)
			return true;

		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->data.data(), this->data.size(), &md5_compare[0]);
		return(this->header.data_header.CheckChecksum(md5_compare));

	} else if (target == 4) {
		if (this->strides.size() == 0)
			return true;

		uint8_t md5_compare[MD5_DIGEST_LENGTH];
		algorithm::VariantDigestManager::GenerateMd5(this->strides.data(), this->strides.size(), &md5_compare[0]);
		return(this->header.stride_header.CheckChecksum(md5_compare));
	}
	return true;
}

bool yon1_dc_t::CheckUniformity(void) {
	if (data_uncompressed.size() == 0)
		return false;

	if (this->header.n_entries == 0)
		return false;

	if (this->header.data_header.controller.type == YON_TYPE_CHAR)
		return false;

	if (this->header.data_header.HasMixedStride())
		return false;

	// We know the stride cannot be uniform if
	// the stride size is uneven
	const int16_t& stride_size = this->header.data_header.stride;
	if (stride_size == -1)
		return false;

	uint32_t stride_update = stride_size;

	uint8_t word_width = sizeof(char);
	switch(this->header.data_header.controller.type) {
	case YON_TYPE_DOUBLE: stride_update *= sizeof(double);   word_width = sizeof(double);  break;
	case YON_TYPE_FLOAT:  stride_update *= sizeof(float);    word_width = sizeof(float);   break;
	case YON_TYPE_8B:     stride_update *= sizeof(uint8_t);  word_width = sizeof(uint8_t); break;
	case YON_TYPE_16B:    stride_update *= sizeof(uint16_t); word_width = sizeof(uint16_t);break;
	case YON_TYPE_32B:    stride_update *= sizeof(int32_t);  word_width = sizeof(int32_t); break;
	case YON_TYPE_64B:    stride_update *= sizeof(uint64_t); word_width = sizeof(uint64_t);break;
	case YON_TYPE_CHAR:   stride_update *= sizeof(char);     word_width = sizeof(char);    break;
	default: return false; break;
	}

	assert(data_uncompressed.size() % word_width == 0);
	assert((data_uncompressed.size() / word_width) % header.n_additions == 0);
	const uint32_t n_e = data_uncompressed.size() / stride_update;

	const uint64_t first_hash = XXH64(this->data_uncompressed.data(), stride_update, 2147483647);

	uint64_t cumulative_position = stride_update;
	for (uint32_t i = 1; i < n_e; ++i) {
		if (XXH64(&this->data_uncompressed.buffer_[cumulative_position], stride_update, 2147483647) != first_hash) {
			return(false);
		}
		cumulative_position += stride_update;
	}
	// std::cerr << "n_entries: " << this->header.n_entries << "/" << n_e << " -> " << cumulative_position << "/" << this->data_uncompressed.size() << " stride: " << stride_update << " word " << (int)word_width << std::endl;
	// std::cerr << "controller=" << this->header.data_header.controller.type << "," << this->header.data_header.controller.uniform << "," << this->header.data_header.controller.signedness << "," << this->header.data_header.controller.mixedStride << std::endl;
	// std::cerr << "Data size=" << data_uncompressed.size() << " with stride " << stride_size << " and stride update=" << stride_update << std::endl;

	assert(cumulative_position == this->data_uncompressed.size());

	this->header.n_entries   = 1;
	this->header.n_strides   = 0;
	this->data_uncompressed.n_chars_                = stride_size * word_width;
	this->header.data_header.uLength                = stride_size * word_width;
	this->header.data_header.cLength                = stride_size * word_width;
	this->header.data_header.controller.uniform     = true;
	this->header.data_header.controller.mixedStride = false;
	this->header.data_header.controller.encoder     = YON_ENCODE_NONE;
	return(true);
}

bool yon1_dc_t::CheckUniformity(const uint32_t n_samples) {
	if (data_uncompressed.size() == 0)
		return false;

	if (this->header.n_entries == 0)
		return false;

	if (this->header.data_header.controller.type == YON_TYPE_CHAR)
		return false;

	if (this->header.data_header.HasMixedStride())
		return false;

	// We know the stride cannot be uniform if
	// the stride size is uneven
	const int16_t& stride_size = this->header.data_header.stride;
	if (stride_size == -1)
		return false;

	uint32_t stride_update = stride_size;

	uint8_t word_width = sizeof(char);
	switch(this->header.data_header.controller.type) {
	case YON_TYPE_DOUBLE: stride_update *= sizeof(double);   word_width = sizeof(double);  break;
	case YON_TYPE_FLOAT:  stride_update *= sizeof(float);    word_width = sizeof(float);   break;
	case YON_TYPE_8B:     stride_update *= sizeof(uint8_t);  word_width = sizeof(uint8_t); break;
	case YON_TYPE_16B:    stride_update *= sizeof(uint16_t); word_width = sizeof(uint16_t);break;
	case YON_TYPE_32B:    stride_update *= sizeof(int32_t);  word_width = sizeof(int32_t); break;
	case YON_TYPE_64B:    stride_update *= sizeof(uint64_t); word_width = sizeof(uint64_t);break;
	case YON_TYPE_CHAR:   stride_update *= sizeof(char);     word_width = sizeof(char);    break;
	default: return false; break;
	}
	stride_update *= n_samples;

	// Hash the first word.
	const uint64_t first_hash = XXH64(this->data_uncompressed.data(), stride_update, 2147483647);

	// Iterate over the remaining words and break if the current hash does not 
	// match the first one.
	uint64_t cumulative_position = stride_update;
	for (uint32_t i = 1; i < this->header.n_entries; ++i) {
		if (XXH64(&this->data_uncompressed.buffer_[cumulative_position], stride_update, 2147483647) != first_hash) {
			return(false);
		}
		cumulative_position += stride_update;
	}
	//std::cerr << cumulative_position << "/" << this->data_uncompressed.size() << " stride: " << stride_update << " word " << (int)word_width << std::endl;
	assert(cumulative_position == this->data_uncompressed.size());

	this->header.n_entries   = 1;
	this->header.n_strides   = 0;
	this->data_uncompressed.n_chars_                = stride_size * word_width * n_samples;
	this->header.data_header.uLength                = stride_size * word_width * n_samples;
	this->header.data_header.cLength                = stride_size * word_width * n_samples;
	this->header.data_header.controller.uniform     = true;
	this->header.data_header.controller.mixedStride = false;
	this->header.data_header.controller.encoder     = YON_ENCODE_NONE;
	return(true);
}

void yon1_dc_t::ReformatInteger() {
	if (data_uncompressed.size() == 0)
		return;

	// Recode integer types only.
	if (!(this->header.data_header.controller.type == YON_TYPE_32B &&
	     this->header.data_header.controller.signedness == true))
	{
		return;
	}

	// Do not recode if the data is uniform.
	if (this->header.data_header.IsUniform())
		return;

	// At this point all integers are signed 32-bit integers.
	assert(this->data_uncompressed.size() % sizeof(int32_t) == 0);
	assert(this->header.n_additions * sizeof(int32_t) == this->data_uncompressed.size());
	const int32_t* const dat  = reinterpret_cast<const int32_t* const>(this->data_uncompressed.data());
	int32_t min_value = dat[0];
	int32_t max_value = dat[0];
	bool has_special = false;

	// Iterate over available data and search for either missingness
	// or sentinel node values. If a match is found trigger the
	// bool flag to signal the use of the recoding procedure.
	for (uint32_t j = 0; j < this->header.n_additions; ++j) {
		if (dat[j] == bcf_int32_missing || dat[j] == bcf_int32_vector_end) {
			has_special = true;
			continue;
		}
		if (dat[j] < min_value) min_value = dat[j];
		if (dat[j] > max_value) max_value = dat[j];
	}

	// If we have missing values then we have to use signed
	// primitives to accommodate this fact.
	uint8_t byte_width = 0;
	if (min_value < 0 || has_special == true) {
		byte_width = ceil((ceil(log2(abs(min_value) + 1 + 2)) + 1) / 8);  // One bit is used for sign, 2 values for missing and end-of-vector
		const uint8_t byte_width_max = ceil((ceil(log2(abs(max_value) + 1 + 2)) + 1) / 8);
		if (byte_width_max > byte_width) {
			byte_width = byte_width_max;
		}
	}
	else byte_width = ceil(ceil(log2(max_value + 1)) / 8);

	// Select the smallest primitive type (word width) that
	// can hold the target data range.
	if (byte_width == 0)     byte_width = 1;
	else if (byte_width >= 3 && byte_width <= 4) byte_width = 4;
	else if (byte_width > 4) byte_width = 8;

	// Setup buffers.
	this->data.reset();
	this->data.resize(this->data_uncompressed.size() + 65536);

	// Is non-negative
	// Also cannot have missing values
	if (min_value >= 0 && has_special == false) {
		this->header.data_header.controller.signedness = 0;

		if (byte_width == 1) {
			this->header.data_header.controller.type = YON_TYPE_8B;

			for (uint32_t j = 0; j < this->header.n_additions; ++j) {
				this->data += (uint8_t)dat[j];
			}

		} else if (byte_width == 2) {
			this->header.data_header.controller.type = YON_TYPE_16B;

			for (uint32_t j = 0; j < this->header.n_additions; ++j) {
				this->data += (uint16_t)dat[j];
			}

		} else if (byte_width == 4) {
			this->header.data_header.controller.type = YON_TYPE_32B;

			for (uint32_t j = 0; j < this->header.n_additions; ++j) {
				this->data += (uint32_t)dat[j];
			}

		} else if (byte_width == 8) {
			this->header.data_header.controller.type = YON_TYPE_64B;

			for (uint32_t j = 0; j < this->header.n_additions; ++j)
				this->data += (uint64_t)dat[j];

		} else {
			std::cerr << utility::timestamp("ERROR") << "Illegal primitive type!" << std::endl;
			exit(1);
		}
	}
	// Is negative or Has missing
	else {
		this->header.data_header.controller.signedness = true;

		if (byte_width == 1) {
			this->header.data_header.controller.type = YON_TYPE_8B;

			const int8_t missing = INT8_MIN;
			const int8_t eov     = INT8_MIN + 1;
			for (uint32_t j = 0; j < this->header.n_additions; ++j) {
				if (dat[j] == bcf_int32_missing)         this->data += missing;
				else if (dat[j] == bcf_int32_vector_end) this->data += eov;
				else this->data += (int8_t)dat[j];
			}

		} else if (byte_width == 2) {
			this->header.data_header.controller.type = YON_TYPE_16B;

			const int16_t missing = INT16_MIN;
			const int16_t eov     = INT16_MIN + 1;
			for (uint32_t j = 0; j < this->header.n_additions; ++j) {
				if (dat[j] == bcf_int32_missing)         this->data += missing;
				else if (dat[j] == bcf_int32_vector_end) this->data += eov;
				else this->data += (int16_t)dat[j];
			}

		} else if (byte_width == 4) {
			this->header.data_header.controller.type = YON_TYPE_32B;

			const int32_t missing = INT32_MIN;
			const int32_t eov     = INT32_MIN + 1;
			for (uint32_t j = 0; j < this->header.n_additions; ++j) {
				if (dat[j] == bcf_int32_missing)         this->data += missing;
				else if (dat[j] == bcf_int32_vector_end) this->data += eov;
				else this->data += (int32_t)dat[j];
			}

		} else {
			std::cerr << utility::timestamp("ERROR") << "Illegal primitive type!" << std::endl;
			exit(1);
		}
	}
	assert(this->data.size() % byte_width == 0);
	assert(this->header.n_additions * byte_width == this->data.size());

	memcpy(this->data_uncompressed.data(), this->data.data(), this->data.size());
	this->data_uncompressed.n_chars_ = this->data.size();
	this->header.data_header.uLength        = this->data_uncompressed.size();
	this->data.reset();
}

void yon1_dc_t::ReformatStride() {
	if (this->strides_uncompressed.size() == 0)
		return;

	if (this->header.data_header.HasMixedStride() == false)
		return;

	// Recode integer types
	if (!(this->header.stride_header.controller.type == YON_TYPE_32B &&
	   this->header.stride_header.controller.signedness == 0)) {
		std::cerr << utility::timestamp("ERROR") << "Illegal to have non-uint32_t values at this point: " << this->header.stride_header.controller.type << ":" << this->header.stride_header.controller.signedness << std::endl;
		exit(1);
	}

	// At this point all integers are uint32_t
	const uint32_t* const dat = reinterpret_cast<const uint32_t* const>(this->strides_uncompressed.data());
	assert(this->strides_uncompressed.size() % sizeof(uint32_t) == 0);
	assert(this->strides_uncompressed.size() / sizeof(uint32_t) == this->header.n_strides);

	uint32_t max = 1;
	for (uint32_t j = 0; j < this->header.n_strides; ++j) {
		if (dat[j] > max)
			max = dat[j];
	}

	uint8_t byte_width = ceil(log2(max + 1) / 8);
	if (byte_width == 0) byte_width = 1;
	if (byte_width >= 3 && byte_width <= 4) byte_width = 4;
	if (byte_width > 4)  byte_width = 8;

	// This cannot ever be uniform
	this->strides.reset();
	this->strides.resize(this->strides_uncompressed.size() + 65536);

	if (byte_width == 1) {
		this->header.stride_header.controller.type = YON_TYPE_8B;

		for (uint32_t j = 0; j < this->header.n_strides; ++j) {
			//assert((uint8_t)dat[j] == dat[j]);
			this->strides += (uint8_t)dat[j];
		}

	} else if (byte_width == 2) {
		this->header.stride_header.controller.type = YON_TYPE_16B;

		for (uint32_t j = 0; j < this->header.n_strides; ++j) {
			//assert((uint16_t)dat[j] == dat[j]);
			this->strides += (uint16_t)dat[j];
		}

	} else if (byte_width == 4) {
		this->header.stride_header.controller.type = YON_TYPE_32B;

		for (uint32_t j = 0; j < this->header.n_strides; ++j) {
			//assert((uint32_t)dat[j] == dat[j]);
			this->strides += (uint32_t)dat[j];
		}

	} else if (byte_width == 8) {
		this->header.stride_header.controller.type = YON_TYPE_64B;

		for (uint32_t j = 0; j < this->header.n_strides; ++j) {
			//assert((uint64_t)dat[j] == dat[j]);
			this->strides += (uint64_t)dat[j];
		}

	} else {
		std::cerr << utility::timestamp("ERROR") << "Illegal primitive type!" << std::endl;
		exit(1);
	}
	memcpy(this->strides_uncompressed.data(), this->strides.data(), this->strides.size());
	this->strides_uncompressed.n_chars_ = this->strides.size();
	this->header.stride_header.uLength         = this->strides_uncompressed.size();
	this->strides.reset();
}

uint32_t yon1_dc_t::GetObjectSize(void) const{
	// In case data is encrypted
	if (this->header.data_header.controller.encryption != YON_ENCRYPTION_NONE)
		return(this->data.size());

	uint32_t total_size = this->data.size();
	if (this->header.data_header.HasMixedStride())
		total_size += this->strides.size();

	return(total_size);
}

uint64_t yon1_dc_t::GetObjectSizeUncompressed(void) const{
	uint64_t total_size = this->data_uncompressed.size();
	if (this->header.data_header.HasMixedStride())
		total_size += this->strides_uncompressed.size();

	return(total_size);
}

void yon1_dc_t::UpdateContainer(bool reformat_data, bool reformat_stride) {
	// If the data container Has entries in it but has
	// no actual data then it is a BOOLEAN
	if (this->header.n_entries && this->data_uncompressed.size() == 0) {
		this->header.reset();
		this->header.data_header.controller.type        = YON_TYPE_BOOLEAN;
		this->header.data_header.controller.uniform     = true;
		this->header.data_header.controller.mixedStride = false;
		this->header.data_header.controller.encoder     = YON_ENCODE_NONE;
		this->header.data_header.controller.signedness  = 0;
		this->header.data_header.stride  = 1;
		this->header.data_header.uLength = 0;
		this->header.data_header.cLength = 0;
		this->header.n_strides           = 0;
		return;
	}

	if (this->data_uncompressed.size() == 0)
		return;

	// Check if stream is uniform in content
	if (this->header.data_header.controller.type != YON_TYPE_STRUCT) {
		this->CheckUniformity();
		// Reformat stream to use as small word size as possible
		if (reformat_data) this->ReformatInteger();
	}

	// Set uncompressed length
	this->header.data_header.uLength = this->data_uncompressed.size();

	// If we have mixed striding
	if (this->header.data_header.HasMixedStride()) {
		// Reformat stream to use as small word size as possible
		if (reformat_stride) this->ReformatStride();
		this->header.stride_header.uLength = this->strides_uncompressed.size();
	}
}

void yon1_dc_t::UpdateContainerFormat(bool reformat_data, bool reformat_stride, const uint32_t n_samples) {
	if (this->data_uncompressed.size() == 0)
		return;

	// Check if stream is uniform in content
	this->CheckUniformity();
	// Reformat stream to use as small word size as possible
	if (reformat_data) this->ReformatInteger();

	// Set uncompressed length
	this->header.data_header.uLength = this->data_uncompressed.size();

	// If we have mixed striding
	if (this->header.data_header.HasMixedStride()) {
		// Reformat stream to use as small word size as possible
		if (reformat_stride) this->ReformatStride();
		this->header.stride_header.uLength = this->strides_uncompressed.size();
	}
}

void yon1_dc_t::AddStride(const uint32_t value) {
	// If this is the first stride set
	if (this->header.n_strides == 0) {
		this->header.stride_header.controller.type = YON_TYPE_32B;
		this->header.stride_header.controller.signedness = false;
		this->header.data_header.stride = value;
	}

	// Check if there are different strides
	if (this->header.data_header.HasMixedStride() == false) {
		if (this->header.data_header.stride != value) {
			this->header.data_header.controller.mixedStride = true;
		}
	}

	// Add value
	this->strides_uncompressed += (uint32_t)value;
	++this->header.n_strides;
}

bool yon1_dc_t::Add(const uint8_t& value) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if (!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added uint8_t" << std::endl;
		exit(1);
		return false;
	}
	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool yon1_dc_t::Add(const uint16_t& value) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if (!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added uint16_t" << std::endl;
		exit(1);
		return false;
	}
	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool yon1_dc_t::Add(const uint32_t& value) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if (!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added uint32_t" << std::endl;
		exit(1);
		return false;
	}
	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool yon1_dc_t::Add(const int8_t& value) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}

	if (!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added int8_t" << std::endl;
		exit(1);
		return false;
	}

	if (value == bcf_int8_vector_end) {
		this->data_uncompressed += (int32_t)bcf_int32_vector_end;
		++this->header.n_additions;
		//std::cerr << "value is int8eov" << std::endl;
		return(true);
	}

	if (value == bcf_int8_missing) {
		this->data_uncompressed += (int32_t)bcf_int32_missing;
		++this->header.n_additions;
		//std::cerr << "value is int8miss" << std::endl;
		return(true);
	}

	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool yon1_dc_t::Add(const int16_t& value) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if (!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added int16_t" << std::endl;
		exit(1);
		return false;
	}

	if (value == bcf_int16_vector_end) {
		this->data_uncompressed += (int32_t)bcf_int32_vector_end;
		++this->header.n_additions;
		//std::cerr << "value is int16eov" << std::endl;
		return(true);
	}

	if (value == bcf_int16_missing) {
		this->data_uncompressed += (int32_t)bcf_int32_missing;
		++this->header.n_additions;
		//std::cerr << "value is int16miss" << std::endl;
		return(true);
	}

	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool yon1_dc_t::Add(const int32_t& value) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_32B);
		this->header.data_header.controller.signedness = false;
	}
	if (!this->CheckInteger()) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added int32_t" << std::endl;
		exit(1);
		return false;
	}

	if (value == bcf_int32_vector_end) {
		this->data_uncompressed += (int32_t)bcf_int32_vector_end;
		++this->header.n_additions;
		//std::cerr << "value is int32eov" << std::endl;
		return(true);
	}

	if (value == bcf_int32_missing) {
		this->data_uncompressed += (int32_t)bcf_int32_missing;
		++this->header.n_additions;
		//std::cerr << "value is int32miss" << std::endl;
		return(true);
	}

	this->data_uncompressed += (int32_t)value;
	++this->header.n_additions;
	return(true);
}

bool yon1_dc_t::Add(const uint64_t& value) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_64B);
		this->header.data_header.controller.signedness = false;
	}

	// Make checks
	if (!this->header.data_header.controller.CompareTypeSign(YON_TYPE_64B, false)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added uint64_t" << std::endl;
		return false;
	}

	this->data_uncompressed += (uint64_t)value;
	++this->header.n_additions;
	return(true);
}

bool yon1_dc_t::Add(const int64_t& value) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_64B);
		this->header.data_header.controller.signedness = true;
	}


	// Make checks
	if (!this->header.data_header.controller.CompareTypeSign(YON_TYPE_64B, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added int64_t" << std::endl;
		return false;
	}

	this->data_uncompressed += (int64_t)value;
	++this->header.n_additions;
	//++this->n_entries;
	return(true);
}

bool yon1_dc_t::Add(const float& value) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_FLOAT);
		this->header.data_header.controller.signedness = true;
	}

	// Make checks
	if (!this->header.data_header.controller.CompareTypeSign(YON_TYPE_FLOAT, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added FLOAT" << std::endl;
		return false;
	}

	this->data_uncompressed += (float)value;
	++this->header.n_additions;
	return(true);
}

bool yon1_dc_t::Add(const double& value) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_DOUBLE);
		this->header.data_header.controller.signedness = true;
	}

	// Make checks
	if (!this->header.data_header.controller.CompareTypeSign(YON_TYPE_DOUBLE, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added DOUBLE" << std::endl;
		return false;
	}

	this->data_uncompressed += (double)value;
	++this->header.n_additions;
	return(true);
}

bool yon1_dc_t::AddCharacter(const char& value) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_CHAR);
		this->header.data_header.controller.signedness = true;
	}

	// Make checks
	if (!this->header.data_header.controller.CompareTypeSign(YON_TYPE_CHAR, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added CHAR" << std::endl;
		return false;
	}

	this->data_uncompressed += (char)value;
	++this->header.n_additions;
	return(true);
}

bool yon1_dc_t::AddCharacter(const char* const string, const uint32_t l_string) {
	if (this->header.data_header.controller.encoder == YON_ENCODE_NONE && this->header.n_entries == 0) {
		this->header.data_header.SetType(YON_TYPE_CHAR);
		this->header.data_header.controller.signedness = true;
		//std::cerr << "triggering: string" << std::endl;
	}

	// Make checks
	if (!this->header.data_header.controller.CompareTypeSign(YON_TYPE_CHAR, true)) {
		std::cerr << "Primitive type -> local: " << (int)this->header.data_header.controller.type << " added CHAR" << std::endl;
		return false;
	}

	this->data_uncompressed.Add(string, l_string);
	this->header.n_additions += l_string;
	return(true);
}

std::ostream& operator<<(std::ostream& stream, const yon1_dc_t& entry) {
	stream << entry.data;
	if (entry.header.data_header.HasMixedStride())
		stream << entry.strides;

	return(stream);
}

std::istream& operator>>(std::istream& stream, yon1_dc_t& entry) {
	if (entry.header.data_header.controller.encryption == YON_ENCRYPTION_NONE) {
		entry.data.resize(entry.header.data_header.cLength);
		stream.read(entry.data.data(), entry.header.data_header.cLength);
		entry.data.n_chars_ = entry.header.data_header.cLength;

		if (entry.header.data_header.HasMixedStride()) {
			entry.strides.resize(entry.header.stride_header.cLength);
			stream.read(entry.strides.data(), entry.header.stride_header.cLength);
			entry.strides.n_chars_ = entry.header.stride_header.cLength;
		}
	} else { // Data is encrypted
		entry.data.resize(entry.header.data_header.eLength);
		stream.read(entry.data.data(), entry.header.data_header.eLength);
		entry.data.n_chars_ = entry.header.data_header.eLength;
	}
	return(stream);
}

}
