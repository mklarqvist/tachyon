/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TACHYON_ENCRYPTION_H_
#define TACHYON_ENCRYPTION_H_

#include <memory>

// encryption
#include "tachyon.h"
#include "buffer.h"
#include "variant_container.h"

namespace tachyon {

/**<
 * This structure represent the basic key-unit in an encryption keychain. All
 * possible encryption types extend this interface with their own specific
 * requirements.
 */
struct KeychainKey {
public:
	KeychainKey();
	virtual ~KeychainKey() = default;

	friend yon_buffer_t& operator+=(yon_buffer_t& buffer, const KeychainKey& key);
	friend std::ostream& operator<<(std::ostream& stream, const KeychainKey& key);
	friend std::istream& operator>>(std::istream& stream, KeychainKey& key);

	inline TACHYON_ENCRYPTION GetType(void) const { return(TACHYON_ENCRYPTION(this->encryption_type)); }
	inline uint64_t& GetIdx(void) { return(this->field_id); }
	inline const uint64_t& GetIdx(void) const { return(this->field_id); }

	// Virtual functions.
	virtual std::ostream& print(std::ostream& stream) =0;
	virtual KeychainKey* Clone() const =0;
	virtual yon_buffer_t& AddToBuffer(yon_buffer_t& buffer) const =0;
	virtual std::ostream& WriteToStream(std::ostream& stream) const =0;
	virtual std::istream& ReadFromStream(std::istream& stream) =0;

public:
	uint8_t   encryption_type;
	uint64_t  field_id;
};

template <uint8_t KeyLength = 32, uint8_t IVLength = 16>
struct KeychainKeyBase : public KeychainKey {
public:
	typedef KeychainKeyBase self_type;
	typedef KeychainKey     parent_type;
	typedef yon_buffer_t buffer_type;

public:
	KeychainKeyBase() = default;
	virtual ~KeychainKeyBase() = default;

	KeychainKeyBase(const KeychainKeyBase& other) :
		parent_type(other)
	{
		memcpy(&this->key[0], &other.key[0], KeyLength);
		memcpy(&this->iv[0],  &other.iv[0],  IVLength);
	}

	KeychainKeyBase& operator=(const KeychainKeyBase& other) {
		this->field_id = other.field_id;
		this->encryption_type = other.encryption_type;
		memcpy(&this->key[0], &other.key[0], KeyLength);
		memcpy(&this->iv[0],  &other.iv[0],  IVLength);
		return(*this);
	}

	virtual inline KeychainKey* Clone() const { return new self_type(*this); }

	virtual yon_buffer_t& AddToBuffer(yon_buffer_t& buffer) const {
		buffer += this->encryption_type;
		buffer += this->field_id;
		buffer.Add((const char*)&this->key[0], KeyLength);
		buffer.Add((const char*)&this->iv[0],  IVLength);
		return(buffer);
	}

	virtual std::ostream& WriteToStream(std::ostream& stream) const {
		stream.write(reinterpret_cast<const char*>(&this->encryption_type), sizeof(uint8_t));
		stream.write(reinterpret_cast<const char*>(&this->field_id), sizeof(uint64_t));
		stream.write((const char*)&this->key[0], KeyLength);
		stream.write((const char*)&this->iv[0],  IVLength);
		return(stream);
	}

	virtual std::istream& ReadFromStream(std::istream& stream) {
		//stream.read(reinterpret_cast<char*>(&key.encryption_type), sizeof(uint8_t));
		stream.read(reinterpret_cast<char*>(&this->field_id), sizeof(uint64_t));
		stream.read((char*)&this->key[0], KeyLength);
		stream.read((char*)&this->iv[0],  IVLength);
		return(stream);
	}

	virtual std::ostream& print(std::ostream& stream) {
		for(uint32_t i = 0; i < KeyLength; ++i) stream << std::hex << (int)this->key[i];
		stream << '\t';
		for(uint32_t i = 0; i < IVLength; ++i) stream << std::hex << (int)this->iv[i];
		return(stream);
	}

	inline bool InitiateRandom(void) {
		RAND_bytes(&this->key[0], KeyLength);
		RAND_bytes(&this->iv[0],  IVLength);
		return true;
	}

public:
	uint8_t key[KeyLength];  // 256 bit key
	uint8_t iv[IVLength];   // 128 bit initiation vector
};

template <uint8_t KeyLength = 32, uint8_t IVLength = 16, uint8_t TagLength = 16>
struct KeychainKeyGCM : public KeychainKeyBase<KeyLength, IVLength>{
public:
	typedef KeychainKeyGCM self_type;
	typedef KeychainKeyBase<KeyLength, IVLength> parent_type;

public:
	KeychainKeyGCM() {
		this->encryption_type = YON_ENCRYPTION_AES_256_GCM;
	}
	~KeychainKeyGCM() {}

	KeychainKeyGCM(const KeychainKeyGCM& other) :
		parent_type(other)
	{
		this->encryption_type = YON_ENCRYPTION_AES_256_GCM;
		memcpy(&this->tag[0], &other.tag[0], TagLength);
	}

	KeychainKeyGCM& operator=(const KeychainKeyGCM& other) {
		this->field_id = other.field_id;
		this->encryption_type = other.encryption_type;
		memcpy(&this->key[0], &other.key[0], KeyLength);
		memcpy(&this->iv[0],  &other.iv[0],  IVLength);
		memcpy(&this->tag[0], &other.tag[0], TagLength);
		return(*this);
	}

	inline KeychainKey* Clone() const { return new self_type(*this); }

	std::ostream& print(std::ostream& stream) {
		for(uint32_t i = 0; i < KeyLength; ++i) stream << std::hex << (int)this->key[i];
		stream << '\t';
		for(uint32_t i = 0; i < IVLength; ++i) stream << std::hex << (int)this->iv[i];
		stream << '\t';
		for(uint32_t i = 0; i < TagLength; ++i) stream << std::hex << (int)this->tag[i];
		stream << std::dec;
		return(stream);
	}

	yon_buffer_t& AddToBuffer(yon_buffer_t& buffer) const {
		buffer += this->encryption_type;
		buffer += this->field_id;
		buffer.Add((const char*)&this->key[0], KeyLength);
		buffer.Add((const char*)&this->iv[0],  IVLength);
		buffer.Add((const char*)&this->tag[0], TagLength);
		return(buffer);
	}

	std::ostream& WriteToStream(std::ostream& stream) const {
		stream.write(reinterpret_cast<const char*>(&this->encryption_type), sizeof(uint8_t));
		stream.write(reinterpret_cast<const char*>(&this->field_id), sizeof(uint64_t));
		stream.write((const char*)&this->key[0], KeyLength);
		stream.write((const char*)&this->iv[0],  IVLength);
		stream.write((const char*)&this->tag[0], TagLength);
		return(stream);
	}

	std::istream& ReadFromStream(std::istream& stream) {
		//stream.read(reinterpret_cast<char*>(&key.encryption_type), sizeof(uint8_t));
		stream.read(reinterpret_cast<char*>(&this->field_id), sizeof(uint64_t));
		stream.read((char*)&this->key[0], KeyLength);
		stream.read((char*)&this->iv[0],  IVLength);
		stream.read((char*)&this->tag[0], TagLength);
		return(stream);
	}

public:
	uint8_t tag[TagLength];  // validation tag for gcm
};

/**<
 * Primary encryption/decryption keychain for a tachyon
 * archive. The keychain does NOT contain any information
 * to ascertain the relationship between the keychain and
 * the correct archive.
 */
class Keychain {
public:
    typedef Keychain           self_type;
	typedef std::size_t        size_type;
    typedef KeychainKey*       value_type;
    typedef value_type&        reference;
    typedef const value_type&  const_reference;
    typedef value_type*        pointer;
    typedef const value_type*  const_pointer;
    typedef yon1_vb_t          variant_block_type;
    typedef yon1_dc_t          container_type;
    typedef std::unordered_map<uint64_t, uint32_t> htable_type;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
    Keychain();
    Keychain(const uint32_t start_capacity);
    Keychain(const self_type& other);
    Keychain(Keychain&& other) = delete;
    Keychain& operator=(Keychain&& other) = delete;
    Keychain& operator=(const self_type& other);
    ~Keychain();
    Keychain& operator+=(const self_type& other);

   	// Element access
	inline reference at(const size_type& position) { return(this->entries_[position]); }
	inline const_reference at(const size_type& position) const { return(this->entries_[position]); }
	inline reference operator[](const size_type& position) { return(this->entries_[position]); }
	inline const_reference operator[](const size_type& position) const { return(this->entries_[position]); }
	inline pointer data(void) { return(this->entries_); }
	inline const_pointer data(void) const { return(this->entries_); }
	inline reference front(void) { return(this->entries_[0]); }
	inline const_reference front(void) const { return(this->entries_[0]); }
	inline reference back(void) { return(this->entries_[this->n_entries_ - 1]); }
	inline const_reference back(void) const { return(this->entries_[this->n_entries_ - 1]); }

	// Capacity
	inline bool empty(void) const { return(this->n_entries_ == 0); }
	inline const size_type& size(void) const { return(this->n_entries_); }
	inline const size_type& capacity(void) const { return(this->n_capacity_); }

	// Iterator
	inline iterator begin() { return iterator(&this->entries_[0]); }
	inline iterator end() { return iterator(&this->entries_[this->n_entries_]); }
	inline const_iterator begin() const { return const_iterator(&this->entries_[0]); }
	inline const_iterator end() const { return const_iterator(&this->entries_[this->n_entries_]); }
	inline const_iterator cbegin() const { return const_iterator(&this->entries_[0]); }
	inline const_iterator cend() const { return const_iterator(&this->entries_[this->n_entries_]); }

	void operator+=(KeychainKey& keychain);

	void resize(const size_type new_capacity);
	void resize(void);

	/**<
	 * Generate a random 64-bit hash identifier from random bytes. If this
	 * identifier has not been used previous then return it. Otherwise
	 * keep trying until an unused one is found.
	 * @param store Flag toggling whether the generated identifier should be stored or not.
	 * @return      Returns the 64-bit hash identifier.
	 */
	uint64_t GetRandomHashIdentifier(const bool store = true);

	/**<
	 * Accessor for checking if a target value is found in the hash table
	 * of hashed identifier. If the target is found it updates the match
	 * reference in-place and returns TRUE. Otherwise return FALSE and
	 * perform no update.
	 * @param value Src value of interest.
	 * @param match Dst value reference of match if found.
	 * @return      Returns TRUE if found or FALSE otherwise.
	 */
	bool GetHashIdentifier(const uint64_t& value, uint32_t& match);

private:
	/**<
	 * Construct a hash-table mapping from a 64-bit field identifier
	 * to the local offset of that entry in the keychain. Used to map
	 * from block->field->keychain entry.
	 * @return Returns TRUE if successful or FALSE otherwise.
	 */
	bool BuildHashTable(void);

	friend std::ostream& operator<<(std::ostream& stream, const self_type& keychain);
	friend std::istream& operator>>(std::istream& stream, self_type& keychain);

private:
    size_type   n_entries_;
    size_type   n_capacity_;
    pointer     entries_;
    htable_type htable_;

    // Pimpl idiom to mask spinlock dependency.
	class KeychainImpl;
	std::unique_ptr<KeychainImpl> mImpl;
};

class EncryptionDecorator {
public:
	typedef EncryptionDecorator self_type;
	typedef yon1_vb_t     variant_block_type;
	typedef yon1_dc_t     stream_container;
	typedef yon_buffer_t  buffer_type;
	typedef Keychain      keychain_type;

public:
	EncryptionDecorator();
	~EncryptionDecorator();

	/**<
	 * Encrypt a provided variant block using the desired encryption method
	 * and store it in the dst keychain.
	 * @param block           Src variant block to encrypt available fields in.
	 * @param keychain        Dst keychain to store the encryption keys.
	 * @param encryption_type Target encryption algorithm to use.
	 * @return                Returns TRUE upon success or FALSE otherwise.
	 */
	bool Encrypt(variant_block_type& block, keychain_type& keychain, TACHYON_ENCRYPTION encryption_type);

	/**<
	 * Decrypt a provided variant block using the provided keychain.
	 * @param block    Src variant block to decrypt data for.
	 * @param keychain Src keychain holding the decryption keys.
	 * @return         Returns TRUE upon success or FALSE otherwise.
	 */
	bool Decrypt(variant_block_type& block, keychain_type& keychain);

private:
	// Pimpl idiom
	class EncryptionDecoratorImpl;
	std::unique_ptr<EncryptionDecoratorImpl> mImpl;
};

}

#endif /* ALGORITHM_ENCRYPTION_ENCRYPTION_H_ */
