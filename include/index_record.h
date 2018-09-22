/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TACHYON_INDEX_RECORD_H_
#define TACHYON_INDEX_RECORD_H_

#include <fstream>
#include <limits>

namespace tachyon {

struct yon1_idx_rec {
public:
	typedef yon1_idx_rec self_type;

public:
	yon1_idx_rec() :
		block_id(0),
		contig_id(-1),
		n_variants(0),
		byte_offset(0),
		byte_offset_end(0),
		min_position(0),
		max_position(0),
		min_bin(std::numeric_limits<int32_t>::max()),
		max_bin(0)
	{}

	bool operator!=(const self_type& other) const{
		if(this->block_id         != other.block_id)         return true;
		if(this->contig_id        != other.contig_id)        return true;
		if(this->n_variants       != other.n_variants)      return true;
		if(this->byte_offset      != other.byte_offset)     return true;
		if(this->byte_offset_end  != other.byte_offset_end) return true;
		if(this->min_position     != other.min_position)     return true;
		if(this->max_position     != other.max_position)     return true;
		if(this->min_bin          != other.min_bin)          return true;
		if(this->max_bin          != other.max_bin)          return true;
		return false;
	}

	inline bool operator==(const self_type& other) const{ return(!(*this != other)); }

	bool operator<(const self_type& other) const{
		if(this->block_id  < other.block_id)  return true;
		if(this->block_id  > other.block_id)  return false;
		if(this->contig_id < other.contig_id) return true;
		if(this->contig_id > other.contig_id) return false;
		return true;
	}

	bool operator<=(const self_type& other) const{
		if(this->block_id  <= other.block_id)  return true;
		if(this->block_id  >  other.block_id)  return false;
		if(this->contig_id <= other.contig_id) return true;
		if(this->contig_id >  other.contig_id) return false;
		return true;
	}

	~yon1_idx_rec(){}

	void reset(void){
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

	std::ostream& print(std::ostream& stream) const{
		stream << block_id << '\t' << contig_id << '\t' << n_variants << "\tOffset: " << byte_offset << '-' << byte_offset_end << " Position: " << min_position << '-' << max_position << " Bins: " << min_bin << '-' << max_bin;
		return(stream);
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.block_id),         sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.contig_id),        sizeof(int32_t));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),      sizeof(uint32_t));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset),     sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end), sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.min_position),     sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.max_position),     sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&entry.min_bin),          sizeof(int32_t));
		stream.write(reinterpret_cast<const char*>(&entry.max_bin),          sizeof(int32_t));

		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.block_id),         sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.contig_id),        sizeof(int32_t));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),      sizeof(uint32_t));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset),     sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end), sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.min_position),     sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.max_position),     sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&entry.min_bin),          sizeof(int32_t));
		stream.read(reinterpret_cast<char*>(&entry.max_bin),          sizeof(int32_t));

		return(stream);
	}

public:
	uint64_t block_id; // this data is duplicated
	int32_t  contig_id; // this data is duplicated
	uint32_t n_variants; // this data is duplicated
	uint64_t byte_offset; // tellg() position in stream for start of record in Tachyon file
	uint64_t byte_offset_end;// tellg() position in stream for start of record in Tachyon file
	uint64_t min_position;   // smallest bp position
	uint64_t max_position;  // largest  bp position
	int32_t  min_bin;
	int32_t  max_bin;
};

}

#endif
