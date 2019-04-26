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

#include "buffer.h"

namespace tachyon {

struct yon1_idx_rec {
public:
	typedef yon1_idx_rec self_type;

public:
	yon1_idx_rec();
	~yon1_idx_rec();

	// Basic comparators.
	bool operator<(const self_type& other) const;
	bool operator<=(const self_type& other) const;
	bool operator!=(const self_type& other) const;
	inline bool operator==(const self_type& other) const { return(!(*this != other)); }

	/**<
	 * Reset members of this struct to their original values
	 * without releasing memory.
	 */
	void reset(void);

	/**<
	 * Print out this record into the dst stream. Used for debugging purposed
	 * only.
	 * @param stream Dst stream.
	 * @return       Returns a reference to the dst stream.
	 */
	std::ostream& Print(std::ostream& stream) const;

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry);
	friend std::istream& operator>>(std::istream& stream, self_type& entry);

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
