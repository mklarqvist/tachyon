#include "io/basic_buffer.h"
#include "keychain_key.h"

namespace tachyon{

io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const KeychainKey& key){
	return(key.AddToBuffer(buffer));
}

std::ostream& operator<<(std::ostream& stream, const KeychainKey& key){
	return(key.WriteToStream(stream));
}

std::istream& operator>>(std::istream& stream, KeychainKey& key){
	return(key.ReadFromStream(stream));
}

}
