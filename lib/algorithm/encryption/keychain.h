#ifndef ALGORITHM_ENCRYPTION_KEYCHAIN_H_
#define ALGORITHM_ENCRYPTION_KEYCHAIN_H_

#include "containers/variant_block.h"
#include "io/basic_buffer.h"
#include "keychain_key.h"
#include "containers/components/generic_iterator.h"
#include "algorithm/spinlock.h"

namespace tachyon{

/**<
 * Primary encryption/decryption keychain for a tachyon
 * archive. The keychain does NOT contain any information
 * to ascertain the relationship between the keychain and
 * the correct archive.
 */
class Keychain {
public:
    typedef Keychain                  self_type;
	typedef std::size_t               size_type;
    typedef KeychainKey*              value_type;
    typedef value_type&               reference;
    typedef const value_type&         const_reference;
    typedef value_type*               pointer;
    typedef const value_type*         const_pointer;
    typedef containers::VariantBlock  variant_block_type;
    typedef containers::DataContainer container_type;
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
	inline reference at(const size_type& position){ return(this->entries_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->entries_[position]); }
	inline reference operator[](const size_type& position){ return(this->entries_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->entries_[position]); }
	inline pointer data(void){ return(this->entries_); }
	inline const_pointer data(void) const{ return(this->entries_); }
	inline reference front(void){ return(this->entries_[0]); }
	inline const_reference front(void) const{ return(this->entries_[0]); }
	inline reference back(void){ return(this->entries_[this->n_entries_ - 1]); }
	inline const_reference back(void) const{ return(this->entries_[this->n_entries_ - 1]); }

	// Capacity
	inline bool empty(void) const{ return(this->n_entries_ == 0); }
	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	// Iterator
	inline iterator begin(){ return iterator(&this->entries_[0]); }
	inline iterator end(){ return iterator(&this->entries_[this->n_entries_]); }
	inline const_iterator begin() const{ return const_iterator(&this->entries_[0]); }
	inline const_iterator end() const{ return const_iterator(&this->entries_[this->n_entries_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->entries_[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->entries_[this->n_entries_]); }

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
    algorithm::SpinLock spinlock_;
};

}

#endif /* ALGORITHM_ENCRYPTION_KEYCHAIN_H_ */
