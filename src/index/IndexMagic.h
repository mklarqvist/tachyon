#ifndef INDEX_INDEXMAGIC_H_
#define INDEX_INDEXMAGIC_H_

namespace Tachyon{
namespace IO{

template <U16 length>
struct MAGICBase{
	typedef MAGICBase self_type;

public:
	MAGICBase(){}	// for reading
	MAGICBase(const char* target){ memcpy(&this->MAGIC[0], target, length); } // for writing
	MAGICBase(const self_type& other){ memcpy(&this->MAGIC[0], &other.MAGIC[0], length); }
	virtual ~MAGICBase(){}

	friend std::istream& operator>>(std::istream& stream, self_type& base){
		stream.read(base.MAGIC, length);
		return(stream);
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& base){
		stream.write(base.MAGIC, length);
		return stream;
	}

	virtual inline bool validate(const char* match) const{ return(strncmp(&this->MAGIC[0], match, length) == 0); }

public:
	char MAGIC[length];
};

template <U16 length>
struct TomahawkHeader : public MAGICBase<length>{
	typedef TomahawkHeader self_type;
	typedef MAGICBase<length> parent_type;

	TomahawkHeader() : version(0), samples(0){} // for reading
	TomahawkHeader(const char* target, const U64 samples) :
		version(Tomahawk::Constants::PROGRAM_VERSION),
		samples(samples)
	{
		memcpy(&this->MAGIC[0], target, length);
	} // for writing

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& header){
		stream.write(header.MAGIC, length);
		stream.write(reinterpret_cast<const char*>(&Tomahawk::Constants::PROGRAM_VERSION), sizeof(float));
		stream.write(reinterpret_cast<const char*>(&header.samples), sizeof(U64));
		return stream;
	}

	friend std::istream& operator>>(std::istream& stream, self_type& header){
		stream.read(header.MAGIC, length);
		stream.read(reinterpret_cast<char *>(&header.version), sizeof(float));
		stream.read(reinterpret_cast<char *>(&header.samples), sizeof(U64));
		return(stream);
	}

public:
	float version;
	U64 samples;
};

}
}

#endif /* TOTEMPOLEMAGIC_H_ */
