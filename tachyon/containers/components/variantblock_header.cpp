#include <cmath>
#include <cassert>
#include "variantblock_header.h"


namespace tachyon{
namespace containers{

VariantBlockHeader::VariantBlockHeader() :
	l_offset_footer(0),
	contigID(-1),
	minPosition(0),
	maxPosition(0),
	n_variants(0)
{}

VariantBlockHeader::~VariantBlockHeader(){}

}
}
