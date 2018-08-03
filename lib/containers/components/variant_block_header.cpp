#include "variant_block_header.h"

#include <cmath>
#include <cassert>


namespace tachyon{
namespace containers{

VariantBlockHeader::VariantBlockHeader() :
	l_offset_footer(0),
	block_hash(0),
	contigID(-1),
	minPosition(0),
	maxPosition(0),
	n_variants(0)
{}

VariantBlockHeader::~VariantBlockHeader(){}

}
}
