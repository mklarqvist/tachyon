#ifndef ALGORITHM_DIGEST_DIGEST_H_
#define ALGORITHM_DIGEST_DIGEST_H_

#include <fstream>

#include "data_container.h"

namespace tachyon{
namespace algorithm{

struct DigitalDigest{
private:
	typedef DigitalDigest self_type;
	typedef yon1_dc_t     container_type;
	typedef yon_buffer_t  buffer_type;

public:
	DigitalDigest(void) :
		hasInitialized(true),
		hasFinished(false)
	{
		this->initialize();
	}

	DigitalDigest(const self_type& other) :
		hasInitialized(other.hasInitialized),
		hasFinished(other.hasFinished),
		data_context(other.data_context),
		stride_context(other.stride_context)
	{
		memcpy(&this->data_digest[0],   other.data_digest,   64);
		memcpy(&this->stride_digest[0], other.stride_digest, 64);
	}

	DigitalDigest& operator=(const self_type& other) {
		this->hasInitialized = other.hasInitialized;
		this->hasFinished    = other.hasFinished;
		this->data_context   = other.data_context;
		this->stride_context = other.stride_context;
		memcpy(&this->data_digest[0],   &other.data_digest[0],   64);
		memcpy(&this->stride_digest[0], &other.stride_digest[0], 64);
		return(*this);
	}

	~DigitalDigest(void) {}

	/**<
	 * Initializes the SHA512 context. Calling this function
	 * is mandatory!
	 * @return
	 */
	inline bool initialize() {
		this->hasInitialized = true;

		if (!SHA512_Init(&this->data_context))   return false;
		if (!SHA512_Init(&this->stride_context)) return false;

		return true;
	}

	/**<
	 *
	 * @param data_buffer
	 * @param stride_buffer
	 * @param has_strides
	 * @return
	 */
	inline bool update(const buffer_type& data_buffer, const buffer_type& stride_buffer, const bool has_strides = true) {
		if (!this->hasInitialized) this->initialize();

		if (!SHA512_Update(&this->data_context, (const uint8_t*)data_buffer.data(), data_buffer.size()))
			return false;

		if (has_strides) {
			if (!SHA512_Update(&this->stride_context, (const uint8_t*)stride_buffer.data(), stride_buffer.size()))
				return false;
		}
		return true;
	}

	/**<
	 *
	 * @return
	 */
	inline bool finalize() {
		if (!this->hasInitialized) return false;
		if (this->hasFinished) return true;

		if (!SHA512_Final(&this->data_digest[0], &this->data_context))
			return false;

		if (!SHA512_Final(&this->stride_digest[0], &this->stride_context))
			return false;

		return true;
	}

	/**<
	 *
	 */
	void clear(void) {
		this->hasFinished = false;
		this->finalize();
		memset(&this->data_digest[0],   0, 64);
		memset(&this->stride_digest[0], 0, 64);
		this->initialize();
		this->hasInitialized = true;
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry) {
		stream.write(reinterpret_cast<const char* const>(&entry.data_digest),   64);
		stream.write(reinterpret_cast<const char* const>(&entry.stride_digest), 64);
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry) {
		stream.read(reinterpret_cast<char*>(&entry.data_digest),   64);
		stream.read(reinterpret_cast<char*>(&entry.stride_digest), 64);
		return(stream);
	}

public:
	bool       hasInitialized;
	bool       hasFinished;
	SHA512_CTX data_context;
	SHA512_CTX stride_context;
	uint8_t       data_digest[64];
	uint8_t       stride_digest[64];
};

struct DigitalDigestPair{
private:
	typedef DigitalDigestPair         self_type;
	typedef DigitalDigest             digest_type;
	typedef yon1_dc_t container_type;

public:
	DigitalDigestPair() {}

	DigitalDigestPair(const self_type& other) :
		uncompressed(other.uncompressed),
		compressed(other.compressed)
	{

	}

	DigitalDigestPair& operator=(const self_type& other) {
		this->uncompressed = other.uncompressed;
		this->compressed = other.compressed;
		return(*this);
	}

	~DigitalDigestPair() {}

	inline bool finalize(void) {
		if (!this->uncompressed.finalize())
			return false;

		if (!this->compressed.finalize())
			return false;

		return true;
	}

	void operator+=(const container_type& container) {
		this->compressed.update(container.data, container.strides, container.header.HasMixedStride());
		this->uncompressed.update(container.data_uncompressed, container.strides_uncompressed, container.header.HasMixedStride());
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry) {
		stream << entry.compressed;
		stream << entry.uncompressed;
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry) {
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



#endif /* ALGORITHM_DIGEST_DIGEST_H_ */
