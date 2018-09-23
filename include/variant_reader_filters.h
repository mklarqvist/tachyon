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
#ifndef TACHYON_VARIANT_READER_FILTERS_H_
#define TACHYON_VARIANT_READER_FILTERS_H_

#include "support/magic_constants.h"
#include "variant_record.h"

namespace tachyon {

struct VariantReaderFilters {
public:
	typedef VariantReaderFilters self_type;

public:
	VariantReaderFilters();
	VariantReaderFilters(const VariantReaderFilters& other) = delete;
	VariantReaderFilters& operator=(const self_type& other) = delete;
	~VariantReaderFilters();

	void AddWrapper(TACHYON_FILTER_FUNCTION filter_function);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const char& r_value, const TACHYON_COMPARATOR_TYPE& comparator);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const int8_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const int16_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const int32_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const int64_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const uint8_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const uint16_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const uint32_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const uint64_t& r_value, const TACHYON_COMPARATOR_TYPE& comparator);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const std::string& r_value, const TACHYON_COMPARATOR_TYPE& comparator);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const float& r_value, const TACHYON_COMPARATOR_TYPE& comparator);
	void Add(TACHYON_FILTER_FUNCTION filter_function, const double& r_value, const TACHYON_COMPARATOR_TYPE& comparator);

	// Capacity
	const size_t& size(void) const;
	const size_t& capacity(void) const;

	/**<
	 * Checks if any filter function require genotype data to be loaded and prepared
	 * @return Returns TRUE if genotype data is required or FALSE otherwise
	 */
	bool HasRequireGenotypes(void) const;
	void SetRequireGenotypes(const bool set = true);

	/**<
	 * Iteratively apply filters in the filter pointer vector
	 * @param objects  Target objects container structure
	 * @param position Target position (relative loci) in the container
	 * @return         Returns TRUE if passes filtering or FALSE otherwise
	 */
	bool Filter(yon1_vnt_t& objects, const uint32_t position) const;

private:
	// Pimpl idiom
	class VariantReaderFiltersImpl;
	std::unique_ptr<VariantReaderFiltersImpl> mImpl;
};

}

#endif
