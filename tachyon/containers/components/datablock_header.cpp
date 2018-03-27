#include "../components/datablock_header.h"

#include <cmath>
#include <cassert>


namespace tachyon{
namespace containers{

DataBlockHeader::DataBlockHeader() :
	l_offset_footer(0),
	contigID(-1),
	minPosition(0),
	maxPosition(0),
	n_variants(0)
{}

DataBlockHeader::~DataBlockHeader(){}

}
}
