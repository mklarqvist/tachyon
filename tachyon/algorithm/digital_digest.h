#ifndef ALGORITHM_DigitalDigest_H_
#define ALGORITHM_DigitalDigest_H_

#include <cstring>
#include <openssl/sha.h>

namespace tachyon{
namespace core{

struct DigitalDigest{
private:
	typedef DigitalDigest             self_type;
	typedef containers::DataContainer container_type;
	typedef io::BasicBuffer           buffer_type;

public:
	DigitalDigest(void){ this->initialize(); }
	DigitalDigest(const self_type& other) :
		data_context(other.data_context),
		stride_context(other.stride_context)
	{
		memcpy(&this->data_digest[0],   other.data_digest,   64);
		memcpy(&this->stride_digest[0], other.stride_digest, 64);
	}

	~DigitalDigest(void){}

	/**<
	 * Initializes the SHA512 context. Calling this function
	 * is mandatory!
	 * @return
	 */
	inline bool initialize(){
		if(!SHA512_Init(&this->data_context))
			return false;

		if(!SHA512_Init(&this->stride_context))
			return false;

		return true;
	}

	/**<
	 *
	 * @param container
	 * @return
	 */
	inline bool update(const buffer_type& data_buffer, const buffer_type& stride_buffer, const bool has_strides = true){
		if(!SHA512_Update(&this->data_context, (BYTE*)data_buffer.data(), data_buffer.size()))
			return false;

		if(has_strides){
			if(!SHA512_Update(&this->stride_context, (BYTE*)stride_buffer.data(), stride_buffer.size()))
				return false;
		}
		return true;
	}

	/**<
	 *
	 * @return
	 */
	inline bool finalize(){
		if(!SHA512_Final(&this->data_digest[0], &this->data_context))
			return false;

		if(!SHA512_Final(&this->stride_digest[0], &this->stride_context))
			return false;

		return true;
	}

	/**<
	 *
	 */
	inline void clear(void){
		this->finalize();
		memset(&this->data_digest[0],   0, 64);
		memset(&this->stride_digest[0], 0, 64);
		this->initialize();
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << std::hex;
		for(U32 j = 0; j < 64; ++j)
			stream << (int)entry.data_digest[j];
		stream << std::dec << '\t' << std::hex;

		for(U32 j = 0; j < 64; ++j)
			stream << (int)entry.stride_digest[j];

		stream << std::hex;
		return(stream);
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char* const>(&entry.data_digest),   64);
		stream.write(reinterpret_cast<const char* const>(&entry.stride_digest), 64);
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.data_digest),   64);
		stream.read(reinterpret_cast<char*>(&entry.stride_digest), 64);
		return(stream);
	}

public:
	SHA512_CTX data_context;
	SHA512_CTX stride_context;
	BYTE       data_digest[64];
	BYTE       stride_digest[64];
};

struct DigitalDigestPair{
private:
	typedef DigitalDigestPair     self_type;
	typedef DigitalDigest         digest_type;

public:
	DigitalDigestPair(){}

	DigitalDigestPair(const self_type& other) :
		uncompressed(other.uncompressed),
		compressed(other.compressed)
	{

	}

	~DigitalDigestPair(){}

	inline bool finalize(void){
		if(!this->uncompressed.finalize())
			return false;

		if(!this->compressed.finalize())
			return false;

		return true;
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.compressed << '\t' << entry.uncompressed;
		return(stream);
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream << entry.compressed;
		stream << entry.uncompressed;
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.compressed;
		stream >> entry.uncompressed;
		return(stream);
	}


public:
	digest_type uncompressed;
	digest_type compressed;
};

}
}

#endif /* ALGORITHM_DigitalDigest_H_ */
