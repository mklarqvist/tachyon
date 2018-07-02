#ifndef BASIC_HELPERS_H_
#define BASIC_HELPERS_H_

#include <string>
#include <vector>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <atomic>
#include <iostream>
#include <cstring>

#include "type_definitions.h"

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
inline void flushErrPointer(std::ostream& stream){ stream << "\33[2K\r\033[0m"; }
std::string datetime();
std::string NumberThousandsSeparator(std::string number);

std::string basePath(const std::string& input);
std::string baseName(const std::string& input);
std::string extensionName(const std::string& input);
std::vector<std::string> filePathBaseExtension(const std::string& input);

template <class T>
T roundUp(T numToRound, int multiple){
	if(multiple == 0) return numToRound;

	int remainder = numToRound % multiple;
	if (remainder == 0) return numToRound;
	return numToRound + multiple - remainder;
}

template <class T>
std::string ToPrettyString(const T& data){
	return utility::NumberThousandsSeparator(std::to_string(data));
}

template <class T>
std::string ToPrettyString(const std::vector<T>& data){
	std::string ret;
	for(U32 i = 0; i < data.size() - 1; ++i){
		ret += utility::NumberThousandsSeparator(std::to_string(data[i]));
		ret += ", ";
	}
	ret += std::to_string(data[data.size()-1]);
	return ret;
}

inline std::string secondsToTimestring(const double& value){
	U32 internalVal = value;
	std::string retVal;
	const U32 hours = internalVal / 3600;
	if(hours > 0) retVal += std::to_string(hours) + "h";
	internalVal %= 3600;
	const U32 min = internalVal / 60;
	if(min > 0) retVal += std::to_string(min) + "m";
	internalVal %= 60;
	const U32 sec = internalVal;
	retVal += std::to_string(sec) + "s";

	return(retVal);
}

S32 char2int(const char& input);
bool HexToBytes(const std::string& hex, uint8_t* target);

template <class T>
std::string toPrettyDiskString(const T value){
	if(value > 1E9){
		return(std::to_string((double)value/1e9) + " GB");
	} else if(value > 1E6){
		return(std::to_string((double)value/1e6) + " MB");
	} else if(value > 1E3){
		return(std::to_string((double)value/1e3) + " KB");
	} else
		return(std::to_string(value) + " B");
}

}
}

#endif /* HELPERS_H_ */
