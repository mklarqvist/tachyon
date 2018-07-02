#include "header_magic.h"

namespace tachyon{
namespace core{

HeaderMagic::HeaderMagic() :
	major_version(tachyon::constants::TACHYON_VERSION_MAJOR),
	minor_version(tachyon::constants::TACHYON_VERSION_MINOR),
	release_version(tachyon::constants::TACHYON_VERSION_RELEASE),
	controller(0),
	n_samples(0),
	n_contigs(0),
	n_info_values(0),
	n_format_values(0),
	n_filter_values(0),
	l_literals(0),
	l_header(0),
	l_header_uncompressed(0)
{
	//memcpy(&this->magic_string[0],
	//	   &tachyon::constants::FILE_HEADER[0],
	//		tachyon::constants::FILE_HEADER_LENGTH);
}

HeaderMagic::HeaderMagic(const self_type& other) :
	major_version(other.major_version),
	minor_version(other.minor_version),
	release_version(other.release_version),
	controller(other.controller),
	n_samples(other.n_samples),
	n_contigs(other.n_contigs),
	n_info_values(other.n_info_values),
	n_format_values(other.n_format_values),
	n_filter_values(other.n_filter_values),
	l_literals(other.l_literals),
	l_header(other.l_header),
	l_header_uncompressed(other.l_header_uncompressed)
{
	//memcpy(&this->magic_string[0],
	//	   &other.magic_string[0],
	//		tachyon::constants::FILE_HEADER_LENGTH);
}

}
}
