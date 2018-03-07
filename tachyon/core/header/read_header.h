#ifndef CORE_HEADER_READ_HEADER_H_
#define CORE_HEADER_READ_HEADER_H_

#include <cstring>

#include "../../support/type_definitions.h"
#include "../../support/MagicConstants.h"
#include "read_header_magic.h"

namespace tachyon{
namespace core{

/**<
 * This class describes data mandatory data in the
 * Tachyon header for Read data
 */
class ReadHeader{
private:
	typedef ReadHeader       self_type;
	typedef ReadHeaderMagic  magic_type;

public:
	explicit ReadHeader(void){}
	~ReadHeader(){}

	// write
	std::ofstream& write(std::ofstream& stream){
		stream << this->header_magic;
		stream.write(&this->literals[0], this->literals.size());
		return(stream);
	}

private:

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.header_magic;
		entry.literals.resize(entry.header_magic.l_literals);
		stream.read(&entry.literals[0], entry.header_magic.l_literals);
		return(stream);
	}

public:
	magic_type  header_magic;
	std::string literals;
};

}
}



#endif /* CORE_HEADER_READ_HEADER_H_ */
