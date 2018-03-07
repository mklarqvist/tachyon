#include <cstring>

#include "../../support/MagicConstants.h"
#include "read_header_magic.h"

namespace tachyon{
namespace core{

ReadHeaderMagic::ReadHeaderMagic() :
	major_version(tachyon::constants::TACHYON_VERSION_MAJOR),
	minor_version(tachyon::constants::TACHYON_VERSION_MINOR),
	release_version(tachyon::constants::TACHYON_VERSION_RELEASE),
	controller(0),
	l_literals(0),
	l_header(0),
	l_header_uncompressed(0)
{
	memcpy(&this->magic_string[0],
		   &tachyon::constants::FILE_HEADER[0],
			tachyon::constants::FILE_HEADER_LENGTH);
}

ReadHeaderMagic::ReadHeaderMagic(const self_type& other) :
	major_version(other.major_version),
	minor_version(other.minor_version),
	release_version(other.release_version),
	controller(other.controller),
	l_literals(other.l_literals),
	l_header(other.l_header),
	l_header_uncompressed(other.l_header_uncompressed)
{
	memcpy(&this->magic_string[0],
		   &other.magic_string[0],
			tachyon::constants::FILE_HEADER_LENGTH);
}

}
}
