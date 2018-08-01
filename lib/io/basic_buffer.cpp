#include "basic_buffer.h"

namespace tachyon{
namespace io{

void SerializeString(const std::string& string, io::BasicBuffer& buffer){
	uint32_t size_helper = string.size();
	buffer += size_helper;
	buffer += string;
}

void DeserializeString(std::string& string, io::BasicBuffer& buffer){
	uint32_t size_helper;
	buffer >> size_helper;
	string.resize(size_helper);
	buffer.read(&string[0], size_helper);
}

}
}
