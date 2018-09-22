#ifndef BASIC_HELPERS_H_
#define BASIC_HELPERS_H_

#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <iostream>
#include <cstring>
#include <regex>

#include "buffer.h"

namespace tachyon{
namespace utility{

int isBigEndian(void);

std::vector<std::string> &split(std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(std::string &s, char delim, const bool keepEmpty = false);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim, const bool keepEmpty = false);
std::vector<std::string> splitLastOf(const std::string& s, const char delim, const bool includeLast = false);

std::string remove_whitespace(std::string& string);
std::string remove_excess_whitespace(const std::string& string);

std::string timestamp(const std::string type);
std::string timestamp(const std::string type, const std::string type2);
std::string datetime();
std::string NumberThousandsSeparator(std::string number);

std::string BasePath(const std::string& input);
std::string BaseName(const std::string& input);
std::string ExtensionName(const std::string& input);
std::vector<std::string> FilePathBaseExtension(const std::string& input);

template <class T>
std::string ToPrettyString(const T& data){
	return utility::NumberThousandsSeparator(std::to_string(data));
}

template <class T>
std::string ToPrettyString(const std::vector<T>& data){
	std::string ret;
	for(uint32_t i = 0; i < data.size() - 1; ++i){
		ret += utility::NumberThousandsSeparator(std::to_string(data[i]));
		ret += ", ";
	}
	ret += std::to_string(data[data.size()-1]);
	return ret;
}

inline std::string SecondsToTimestring(const double& value){
	uint32_t internalVal = value;
	std::string retVal;
	const uint32_t hours = internalVal / 3600;
	if(hours > 0) retVal += std::to_string(hours) + "h";
	internalVal %= 3600;
	const uint32_t min = internalVal / 60;
	if(min > 0) retVal += std::to_string(min) + "m";
	internalVal %= 60;
	const uint32_t sec = internalVal;
	retVal += std::to_string(sec) + "s";

	return(retVal);
}

int32_t char2int(const char& input);
bool HexToBytes(const std::string& hex, uint8_t* target);

template <class T>
std::string ToPrettyDiskString(const T value){
	if(value > 1E9){
		return(std::to_string((double)value/1e9) + " GB");
	} else if(value > 1E6){
		return(std::to_string((double)value/1e6) + " MB");
	} else if(value > 1E3){
		return(std::to_string((double)value/1e3) + " KB");
	} else
		return(std::to_string(value) + " B");
}

void SerializeString(const std::string& string, std::ostream& stream);
void DeserializeString(std::string& string, std::istream& stream);

template <class T>
void SerializePrimitive(const T& value, std::ostream& stream){
	stream.write((const char*)&value, sizeof(T));
}

template <class T>
void DeserializePrimitive(T& value, std::istream& stream){
	stream.read((char*)&value, sizeof(T));
}

}
}

#endif /* HELPERS_H_ */
