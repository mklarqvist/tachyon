#ifndef INDEX_INDEX_INDEX_ENTRY_H_
#define INDEX_INDEX_INDEX_ENTRY_H_

#include "index_record.h"

namespace tachyon{
namespace index{

struct VariantIndexMetaEntry{
public:
	typedef VariantIndexMetaEntry self_type;
	typedef yon1_idx_rec          entry_type;

public:
	VariantIndexMetaEntry() :
		contig_id(-1),
		n_blocks(0),
		n_variants(0),
		byte_offset_begin(0),
		byte_offest_end(0),
		min_position(0),
		max_position(0),
		start_block(0),
		end_block(0)
	{}

	~VariantIndexMetaEntry() {}

	void reset(void) {
		this->contig_id          = -1;
		this->n_blocks          = 0;
		this->n_variants        = 0;
		this->byte_offset_begin = 0;
		this->byte_offest_end   = 0;
		this->min_position       = 0;
		this->max_position       = 0;
		this->start_block       = 0;
		this->end_block         = 0;
	}

	void operator=(const entry_type& entry) {
		this->contig_id         = entry.contig_id;
		this->n_blocks          = 1;
		this->n_variants        = entry.n_variants;
		this->byte_offset_begin = entry.byte_offset;
		this->byte_offest_end   = entry.byte_offset_end;
		this->min_position      = entry.min_position;
		this->max_position      = entry.max_position;
		this->start_block       = entry.block_id;
		this->end_block         = entry.block_id;
	}

	inline bool operator==(const entry_type& entry) const { return(this->contig_id == entry.contig_id); }
	inline bool operator!=(const self_type& other) const { return(!(*this == other)); }

	bool operator==(const self_type& other) const {
		if (this->contig_id          != other.contig_id) return false;
		if (this->n_blocks          != other.n_blocks) return false;
		if (this->n_variants        != other.n_variants) return false;
		if (this->byte_offset_begin != other.byte_offset_begin) return false;
		if (this->byte_offest_end   != other.byte_offest_end) return false;
		if (this->min_position       != other.min_position) return false;
		if (this->max_position       != other.max_position) return false;
		if (this->start_block       != other.start_block) return false;
		if (this->end_block         != other.end_block) return false;
		return true;
	}

	void operator+=(const entry_type& entry) {
		++this->n_blocks;
		this->n_variants     += entry.n_variants;
		this->byte_offest_end = entry.byte_offset_end;
		this->max_position    = entry.max_position;
		this->end_block       = entry.block_id;
	}

	std::ostream& Print(std::ostream& stream) const {
		stream << contig_id << ":" << min_position << "-" << max_position << ", " << n_variants << "," << n_blocks << " @ " << byte_offset_begin << "-" << byte_offest_end << " blocks: " << start_block << "-" << end_block;
		return(stream);
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry) {
		stream.write(reinterpret_cast<const char*>(&entry.contig_id),         sizeof(int32_t));
		stream.write(reinterpret_cast<const char*>(&entry.n_blocks),         sizeof(uint32_t));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),       sizeof(uint32_t));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_begin),sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offest_end),  sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.min_position),      sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.max_position),      sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.start_block),      sizeof(uint32_t));
		stream.write(reinterpret_cast<const char*>(&entry.end_block),        sizeof(uint32_t));

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry) {
		stream.read(reinterpret_cast<char*>(&entry.contig_id),         sizeof(int32_t));
		stream.read(reinterpret_cast<char*>(&entry.n_blocks),         sizeof(uint32_t));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),       sizeof(uint32_t));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_begin),sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.byte_offest_end),  sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.min_position),      sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.max_position),      sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.start_block),      sizeof(uint32_t));
		stream.read(reinterpret_cast<char*>(&entry.end_block),        sizeof(uint32_t));

		return(stream);
	}

public:
	int32_t  contig_id;            // contig ID
	uint32_t n_blocks;           // number of index entries for this contig
	uint64_t n_variants;         // total number of variants
	uint64_t byte_offset_begin;  // start virtual offset
	uint64_t byte_offest_end;    // end virtual offset
	uint64_t min_position;        // smallest bp position observed
	uint64_t max_position;        // largest  bp position observed
	uint32_t start_block;        // start block
	uint32_t end_block;          // end block
};

}
}


#endif /* INDEX_INDEX_INDEX_ENTRY_H_ */
