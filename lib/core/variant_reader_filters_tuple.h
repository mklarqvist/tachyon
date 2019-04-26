#ifndef CORE_VARIANT_READER_FILTERS_TUPLE_H_
#define CORE_VARIANT_READER_FILTERS_TUPLE_H_

#include <regex>

#include "utility.h"
#include "tachyon.h"

namespace tachyon{

struct VariantReaderFiltersTupleInterface {
public:
	VariantReaderFiltersTupleInterface() : filter(false) {}
	VariantReaderFiltersTupleInterface(const bool filter) : filter(filter) {}
	VariantReaderFiltersTupleInterface(const VariantReaderFiltersTupleInterface& other) : filter(other.filter) {}
	virtual ~VariantReaderFiltersTupleInterface() {}

	// Not possible to have a templated pure virtual functionality
	virtual bool applyFilter(const bool& l_value)  const =0;
	virtual bool applyFilter(const uint8_t& l_value)  const =0;
	virtual bool applyFilter(const uint16_t& l_value)   const =0;
	virtual bool applyFilter(const uint32_t& l_value)   const =0;
	virtual bool applyFilter(const uint64_t& l_value)   const =0;
	virtual bool applyFilter(const int8_t& l_value) const =0;
	virtual bool applyFilter(const int16_t& l_value)   const =0;
	virtual bool applyFilter(const int32_t& l_value)   const =0;
	virtual bool applyFilter(const int64_t& l_value)   const =0;
	virtual bool applyFilter(const float& l_value) const =0;
	virtual bool applyFilter(const double& l_value)const =0;
	virtual bool applyFilter(const std::string& l_value) const =0;

public:
	bool filter;
};

template <class ValueClass>
struct VariantReaderFiltersTuple : public VariantReaderFiltersTupleInterface{
public:
	typedef VariantReaderFiltersTuple<ValueClass> self_type;
	typedef bool (self_type::*filter_function)(const ValueClass& target, const ValueClass& limit) const;

	// Bring forward interface function into scope to overload
	using VariantReaderFiltersTupleInterface::applyFilter;

public:
	VariantReaderFiltersTuple() :
		r_value(0),
		comparator(&self_type::__filterGreaterEqual)
	{}

	VariantReaderFiltersTuple(const ValueClass& r_value) :
		VariantReaderFiltersTupleInterface(true),
		r_value(r_value),
		comparator(&self_type::__filterGreaterEqual)
	{}

	VariantReaderFiltersTuple(const ValueClass& r_value, const TACHYON_COMPARATOR_TYPE& comparator) :
		VariantReaderFiltersTupleInterface(true),
		r_value(r_value),
		comparator(nullptr)
	{
		switch(comparator) {
		case(YON_CMP_GREATER):       this->comparator = &self_type::__filterGreater;      break;
		case(YON_CMP_GREATER_EQUAL): this->comparator = &self_type::__filterGreaterEqual; break;
		case(YON_CMP_LESS):          this->comparator = &self_type::__filterLesser;       break;
		case(YON_CMP_LESS_EQUAL):    this->comparator = &self_type::__filterLesserEqual;  break;
		case(YON_CMP_EQUAL):         this->comparator = &self_type::__filterEqual;        break;
		case(YON_CMP_NOT_EQUAL):     this->comparator = &self_type::__filterNotEqual;     break;
		case(YON_CMP_REGEX):
		default:
			std::cerr << utility::timestamp("ERROR","FILTER") << "Numerical filtering operations do not support regular expression operations..." << std::endl;
			this->comparator = &self_type::__filterEqual;
		}
	}

	VariantReaderFiltersTuple(const self_type& other) :
		VariantReaderFiltersTupleInterface(other),
		r_value(other.r_value),
		comparator(other.comparator)
	{}

	void operator()(const ValueClass& r_value) {
		this->filter = true;
		this->r_value = r_value;
	}

	void operator()(const ValueClass& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
		this->filter  = true;
		this->r_value = r_value;

		switch(comparator) {
		case(YON_CMP_GREATER):       this->comparator = &self_type::__filterGreater;      break;
		case(YON_CMP_GREATER_EQUAL): this->comparator = &self_type::__filterGreaterEqual; break;
		case(YON_CMP_LESS):          this->comparator = &self_type::__filterLesser;       break;
		case(YON_CMP_LESS_EQUAL):    this->comparator = &self_type::__filterLesserEqual;  break;
		case(YON_CMP_EQUAL):         this->comparator = &self_type::__filterEqual;        break;
		case(YON_CMP_NOT_EQUAL):     this->comparator = &self_type::__filterNotEqual;     break;
		case(YON_CMP_REGEX):
		default:
			std::cerr << utility::timestamp("ERROR","FILTER") << "Numerical filtering operations do not support regular expression operations..." << std::endl;
			this->comparator = &self_type::__filterEqual;
			return;
		}
	}

	~VariantReaderFiltersTuple() = default;

	inline bool applyFilter(const bool& l_value)  const { return((this->*comparator)(l_value, r_value)); }
	inline bool applyFilter(const uint8_t& l_value)  const { return((this->*comparator)(l_value, r_value)); }
	inline bool applyFilter(const uint16_t& l_value)   const { return((this->*comparator)(l_value, r_value)); }
	inline bool applyFilter(const uint32_t& l_value)   const { return((this->*comparator)(l_value, r_value)); }
	inline bool applyFilter(const uint64_t& l_value)   const { return((this->*comparator)(l_value, r_value)); }
	inline bool applyFilter(const int8_t& l_value) const { return((this->*comparator)(l_value, r_value)); }
	inline bool applyFilter(const int16_t& l_value)   const { return((this->*comparator)(l_value, r_value)); }
	inline bool applyFilter(const int32_t& l_value)   const { return((this->*comparator)(l_value, r_value)); }
	inline bool applyFilter(const int64_t& l_value)   const { return((this->*comparator)(l_value, r_value)); }
	inline bool applyFilter(const float& l_value) const { return((this->*comparator)(l_value, r_value)); }
	inline bool applyFilter(const double& l_value)const { return((this->*comparator)(l_value, r_value)); }
	inline bool applyFilter(const std::string& l_value) const { return(false); }

	// Comparator functions
	inline bool __filterLesser(const ValueClass& target, const ValueClass& limit) const {return(target < limit);}
	inline bool __filterLesserEqual(const ValueClass& target, const ValueClass& limit) const {return(target <= limit);}
	inline bool __filterGreater(const ValueClass& target, const ValueClass& limit) const {return(target > limit);}
	inline bool __filterGreaterEqual(const ValueClass& target, const ValueClass& limit) const {return(target >= limit);}
	inline bool __filterEqual(const ValueClass& target, const ValueClass& limit) const {return(target == limit);}
	inline bool __filterNotEqual(const ValueClass& target, const ValueClass& limit) const {return(target != limit);}

public:
	ValueClass      r_value;
	filter_function comparator;
};

template <>
struct VariantReaderFiltersTuple<std::string> : public VariantReaderFiltersTupleInterface{
public:
	typedef VariantReaderFiltersTuple<std::string> self_type;
	typedef bool (self_type::*filter_function)(const std::string& target, const std::string& limit) const;

public:
	VariantReaderFiltersTuple() :
		comparator(&self_type::__filterGreaterEqual)
	{}

	VariantReaderFiltersTuple(const std::string& r_value) :
		VariantReaderFiltersTupleInterface(true),
		r_value(r_value),
		r_value_regex(std::regex(r_value)),
		comparator(&self_type::__filterGreaterEqual)
	{}

	VariantReaderFiltersTuple(const std::string& r_value, const TACHYON_COMPARATOR_TYPE& comparator) :
		VariantReaderFiltersTupleInterface(true),
		r_value(r_value),
		r_value_regex(std::regex(r_value)),
		comparator(nullptr)
	{
		switch(comparator) {
		case(YON_CMP_GREATER):       this->comparator = &self_type::__filterGreater;      break;
		case(YON_CMP_GREATER_EQUAL): this->comparator = &self_type::__filterGreaterEqual; break;
		case(YON_CMP_LESS):          this->comparator = &self_type::__filterLesser;       break;
		case(YON_CMP_LESS_EQUAL):    this->comparator = &self_type::__filterLesserEqual;  break;
		case(YON_CMP_EQUAL):         this->comparator = &self_type::__filterEqual;        break;
		case(YON_CMP_NOT_EQUAL):     this->comparator = &self_type::__filterNotEqual;     break;
		case(YON_CMP_REGEX):         this->comparator = &self_type::__filterRegex;        break;
		}
	}

	VariantReaderFiltersTuple(const self_type& other) :
		VariantReaderFiltersTupleInterface(other),
		r_value(other.r_value),
		r_value_regex(other.r_value_regex),
		comparator(other.comparator)
	{}

	void operator()(const std::string& r_value) {
		this->filter = true;
		this->r_value = r_value;
		this->r_value_regex = std::regex(r_value);
	}

	void operator()(const std::string& r_value, const TACHYON_COMPARATOR_TYPE& comparator) {
		this->filter  = true;
		this->r_value = r_value;
		this->r_value_regex = std::regex(r_value);

		switch(comparator) {
		case(YON_CMP_GREATER):       this->comparator = &self_type::__filterGreater;      break;
		case(YON_CMP_GREATER_EQUAL): this->comparator = &self_type::__filterGreaterEqual; break;
		case(YON_CMP_LESS):          this->comparator = &self_type::__filterLesser;       break;
		case(YON_CMP_LESS_EQUAL):    this->comparator = &self_type::__filterLesserEqual;  break;
		case(YON_CMP_EQUAL):         this->comparator = &self_type::__filterEqual;        break;
		case(YON_CMP_NOT_EQUAL):     this->comparator = &self_type::__filterNotEqual;     break;
		case(YON_CMP_REGEX):         this->comparator = &self_type::__filterRegex;        break;
		}
	}

	~VariantReaderFiltersTuple() = default;

	inline bool applyFilter(const bool& l_value)  const { return false; }
	inline bool applyFilter(const uint8_t& l_value)  const { return false; }
	inline bool applyFilter(const uint16_t& l_value)   const { return false; }
	inline bool applyFilter(const uint32_t& l_value)   const { return false; }
	inline bool applyFilter(const uint64_t& l_value)   const { return false; }
	inline bool applyFilter(const int8_t& l_value) const { return false; }
	inline bool applyFilter(const int16_t& l_value)   const { return false; }
	inline bool applyFilter(const int32_t& l_value)   const { return false; }
	inline bool applyFilter(const int64_t& l_value)   const { return false; }
	inline bool applyFilter(const float& l_value) const { return false; }
	inline bool applyFilter(const double& l_value)const { return false; }
	inline bool applyFilter(const std::string& l_value) const { return((this->*comparator)(l_value, r_value)); }

	// Comparator functions
	inline bool __filterLesser(const std::string& target, const std::string& limit) const {return(target.size() < limit.size());}
	inline bool __filterLesserEqual(const std::string& target, const std::string& limit) const {return(target.size() <= limit.size());}
	inline bool __filterGreater(const std::string& target, const std::string& limit) const {return(target.size() > limit.size());}
	inline bool __filterGreaterEqual(const std::string& target, const std::string& limit) const {return(target.size() >= limit.size());}
	inline bool __filterEqual(const std::string& target, const std::string& limit) const {return(target == limit);}
	inline bool __filterNotEqual(const std::string& target, const std::string& limit) const {return(target != limit);}

	inline bool __filterRegex(const std::string& target, const std::string& limit) const {
		return(std::regex_search(target, this->r_value_regex, std::regex_constants::match_any));
	}

public:
	std::string     r_value;
	std::regex      r_value_regex;
	filter_function comparator;
};


}



#endif /* CORE_VARIANT_READER_FILTERS_TUPLE_H_ */
