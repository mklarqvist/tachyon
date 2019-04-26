#include "index_record.h"

namespace tachyon {

yon1_idx_rec::yon1_idx_rec() :
	block_id(0), contig_id(-1), n_variants(0),
	byte_offset(0), byte_offset_end(0),
	min_position(0), max_position(0),
	min_bin(std::numeric_limits<int32_t>::max()), max_bin(0)
{}

bool yon1_idx_rec::operator!=(const self_type& other) const {
	if (this->block_id         != other.block_id)         return true;
	if (this->contig_id        != other.contig_id)        return true;
	if (this->n_variants       != other.n_variants)      return true;
	if (this->byte_offset      != other.byte_offset)     return true;
	if (this->byte_offset_end  != other.byte_offset_end) return true;
	if (this->min_position     != other.min_position)     return true;
	if (this->max_position     != other.max_position)     return true;
	if (this->min_bin          != other.min_bin)          return true;
	if (this->max_bin          != other.max_bin)          return true;
	return false;
}

bool yon1_idx_rec::operator<(const self_type& other) const {
	if (this->block_id  < other.block_id)  return true;
	if (this->block_id  > other.block_id)  return false;
	if (this->contig_id < other.contig_id) return true;
	if (this->contig_id > other.contig_id) return false;
	return true;
}

bool yon1_idx_rec::operator<=(const self_type& other) const {
	if (this->block_id  <= other.block_id)  return true;
	if (this->block_id  >  other.block_id)  return false;
	if (this->contig_id <= other.contig_id) return true;
	if (this->contig_id >  other.contig_id) return false;
	return true;
}

yon1_idx_rec::~yon1_idx_rec() {}

void yon1_idx_rec::reset(void) {
	this->block_id         = 0;
	this->contig_id        = -1;
	this->n_variants       = 0;
	this->byte_offset      = 0;
	this->byte_offset_end  = 0;
	this->min_position     = 0;
	this->max_position     = 0;
	this->min_bin          = std::numeric_limits<int32_t>::max();
	this->max_bin          = 0;
}

std::ostream& yon1_idx_rec::Print(std::ostream& stream) const {
	stream << block_id << '\t' << contig_id << '\t' << n_variants << "\tOffset: " << byte_offset << '-' << byte_offset_end << " Position: " << min_position << '-' << max_position << " Bins: " << min_bin << '-' << max_bin;
	return(stream);
}

std::ostream& operator<<(std::ostream& stream, const yon1_idx_rec& entry) {
	utility::SerializePrimitive(entry.block_id,     stream);
	utility::SerializePrimitive(entry.contig_id,    stream);
	utility::SerializePrimitive(entry.n_variants,   stream);
	utility::SerializePrimitive(entry.byte_offset,  stream);
	utility::SerializePrimitive(entry.byte_offset_end, stream);
	utility::SerializePrimitive(entry.min_position, stream);
	utility::SerializePrimitive(entry.max_position, stream);
	utility::SerializePrimitive(entry.min_bin,      stream);
	utility::SerializePrimitive(entry.max_bin,      stream);
	return(stream);
}

std::istream& operator>>(std::istream& stream, yon1_idx_rec& entry) {
	utility::DeserializePrimitive(entry.block_id,     stream);
	utility::DeserializePrimitive(entry.contig_id,    stream);
	utility::DeserializePrimitive(entry.n_variants,   stream);
	utility::DeserializePrimitive(entry.byte_offset,  stream);
	utility::DeserializePrimitive(entry.byte_offset_end, stream);
	utility::DeserializePrimitive(entry.min_position, stream);
	utility::DeserializePrimitive(entry.max_position, stream);
	utility::DeserializePrimitive(entry.min_bin,      stream);
	utility::DeserializePrimitive(entry.max_bin,      stream);
	return(stream);
}

}
