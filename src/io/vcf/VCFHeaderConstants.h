#ifndef VCF_VCFHEADERCONSTANTS_H_
#define VCF_VCFHEADERCONSTANTS_H_

#include <string>
#include "../../support/TypeDefinitions.h"

namespace tachyon{
namespace vcf{
namespace constants{

const std::string HEADER_COLUMN = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
const char VCF_DELIMITER = '\t';

const std::string HEADER_VCF_FORMAT  = "##fileformat=VCFv";
const std::string HEADER_VCF_VERSION = "##fileformat=VCFv4";
const std::string HEADER_FILEFORMAT  = "##fileformat=";
const std::string HEADER_CONTIG      = "##contig=";
const std::string HEADER_ALT         = "##alt=";
const std::string HEADER_INFO        = "##info=";
const std::string HEADER_FILTER      = "##filter=";
const std::string HEADER_FORMAT      = "##format=";

}
}
}


#endif /* VCFHEADERCONSTANTS_H_ */
