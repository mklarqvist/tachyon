#ifndef ALGORITHM_DIGEST_VARIANT_DIGEST_MANAGER_H_
#define ALGORITHM_DIGEST_VARIANT_DIGEST_MANAGER_H_

#include "openssl/md5.h"
#include "digest_manager.h"
#include "variant_container.h"

namespace tachyon{
namespace algorithm{

class VariantDigestManager : public DigestManager{
private:
	typedef VariantDigestManager  self_type;
	typedef DigestManager         parent_type;
	typedef yon1_vb_t             variant_block_type;

public:
	VariantDigestManager();
	VariantDigestManager(const size_type base_capacity);
	VariantDigestManager(const size_type base_capacity, const size_type capacity_info, const size_type capacity_format);
	VariantDigestManager(const self_type& other);
	~VariantDigestManager();

	void Finalize(void);

	inline const_reference atINFO(const uint32_t position) const { return(this->__entries_info[position]); }
	inline const_reference atFORMAT(const uint32_t position) const { return(this->__entries_format[position]); }
	inline reference atINFO(const uint32_t position) { return(this->__entries_info[position]); }
	inline reference atFORMAT(const uint32_t position) { return(this->__entries_format[position]); }

	void operator+=(const variant_block_type& block);

	template <class T>
	static void GenerateMd5(const T& value, uint8_t* dst) {
	   // uint8_t hash[MD5_DIGEST_LENGTH];
		MD5_CTX md5;
		MD5_Init(&md5);
		MD5_Update(&md5, value, sizeof(T));
		MD5_Final(dst, &md5);
	}

	static void GenerateMd5(const char* data, const uint32_t l_data, uint8_t* dst) {
	   // uint8_t hash[MD5_DIGEST_LENGTH];
		MD5_CTX md5;
		MD5_Init(&md5);
		MD5_Update(&md5, data, l_data);
		MD5_Final(dst, &md5);
	}

	friend std::ostream& operator<<(std::ostream& out, const self_type& container) {
		const parent_type* const parent = reinterpret_cast<const parent_type* const>(&container);
		out << *parent;

		out.write((const char* const)reinterpret_cast<const size_type* const>(&container.n_entries_info_), sizeof(size_type));
		out.write((const char* const)reinterpret_cast<const size_type* const>(&container.n_entries_format_), sizeof(size_type));
		for (size_type i = 0; i < container.n_entries_info_; ++i)   out << container.atINFO(i);
		for (size_type i = 0; i < container.n_entries_format_; ++i) out << container.atFORMAT(i);

		return(out);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& container) {
		parent_type* parent = reinterpret_cast<parent_type*>(&container);
		stream >> *parent;

		stream.read((char*)reinterpret_cast<size_type*>(&container.n_entries_info_), sizeof(size_type));
		stream.read((char*)reinterpret_cast<size_type*>(&container.n_capacity_format), sizeof(size_type));

		delete [] container.__entries_info;
		delete [] container.__entries_format;
		container.__entries_info = new value_type[container.n_entries_info_];
		container.__entries_format = new value_type[container.n_entries_format_];
		for (size_type i = 0; i < container.n_entries_info_; ++i)   stream >> container.atINFO(i);
		for (size_type i = 0; i < container.n_entries_format_; ++i) stream >> container.atFORMAT(i);

		return(stream);
	}

private:
	size_type n_entries_info_;
	size_type n_entries_format_;
	size_type n_capacity_info_;
	size_type n_capacity_format;
	pointer   __entries_info;
	pointer   __entries_format;
};

}
}


#endif /* ALGORITHM_DIGEST_VARIANT_DIGEST_MANAGER_H_ */
