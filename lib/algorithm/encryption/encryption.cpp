#include <openssl/evp.h>
#include <openssl/aes.h>
#include <openssl/rand.h>

#include "third_party/xxhash/xxhash.h"

#include "algorithm/spinlock.h"
#include "encryption.h"

namespace tachyon {

KeychainKey::KeychainKey(): field_id(0), encryption_type(YON_ENCRYPTION_NONE){}

// keys
io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const KeychainKey& key){
	return(key.AddToBuffer(buffer));
}

std::ostream& operator<<(std::ostream& stream, const KeychainKey& key){
	return(key.WriteToStream(stream));
}

std::istream& operator>>(std::istream& stream, KeychainKey& key){
	return(key.ReadFromStream(stream));
}

class Keychain::KeychainImpl {
public:
	KeychainImpl() = default;
	~KeychainImpl() = default;

public:
	algorithm::SpinLock spinlock_;
};


// Keychain
Keychain::Keychain() :
	n_entries_(0),
	n_capacity_(100000),
	entries_(new value_type[this->n_capacity_]),
	mImpl(new Keychain::KeychainImpl)
{
}

Keychain::Keychain(const uint32_t start_capacity) :
	n_entries_(0),
	n_capacity_(start_capacity),
	entries_(new value_type[this->n_capacity_]),
	mImpl(new Keychain::KeychainImpl)
{}

Keychain::Keychain(const self_type& other) :
	n_entries_(other.n_entries_),
	n_capacity_(other.n_capacity_),
	entries_(new value_type[this->n_capacity_]),
	mImpl(new Keychain::KeychainImpl)
{
	for(uint32_t i = 0; i < this->size(); ++i)
		this->entries_[i] = other.entries_[i]->Clone();
	this->BuildHashTable();
}

Keychain::~Keychain(){
	for(int i = 0; i < this->size(); ++i)
		delete this->entries_[i];
	delete [] this->entries_;
}

Keychain& Keychain::operator=(const self_type& other){
	for(int i = 0; i < this->size(); ++i)
		delete this->entries_[i];
	delete [] this->entries_;

	this->n_entries_ = other.n_entries_;
	this->n_capacity_ = other.n_capacity_;
	this->entries_ = new value_type[this->n_capacity_];

	// Clone data.
	for(uint32_t i = 0; i < this->size(); ++i)
		this->entries_[i] = other.entries_[i]->Clone();

	// Construct hash table for hashed field-id lookups.
	this->BuildHashTable();
	return(*this);
}

Keychain& Keychain::operator+=(const self_type& other){
	if(this->size() + other.size() >= this->capacity()){
		this->resize(this->size() + other.size() + 1024);
	}

	uint32_t offset = 0;
	for(int i = this->size(); i < this->size() + other.size(); ++i, ++offset){
		this->entries_[i] = other.entries_[offset]->Clone();
	}
	this->n_entries_ += other.size();

	// Reconstruct the hash-table for hashed field-id lookups.
	this->BuildHashTable();

	return(*this);
}

void Keychain::operator+=(KeychainKey& keychain){
	this->mImpl->spinlock_.lock();
	if(this->size() + 1 == this->capacity())
		this->resize();

	this->entries_[this->n_entries_++] = keychain.Clone();
	this->mImpl->spinlock_.unlock();
}

void Keychain::resize(const size_type new_capacity){
	if(new_capacity < this->capacity()){
		this->n_entries_ = new_capacity;
		return;
	}

	pointer old = this->entries_;
	this->entries_ = new value_type[new_capacity];
	for(uint32_t i = 0; i < this->size(); ++i) this->entries_[i] = old[i];
	delete [] old;
	this->n_capacity_ =  new_capacity;
}

void Keychain::resize(void){ this->resize(this->capacity()*2); }

uint64_t Keychain::GetRandomHashIdentifier(const bool store){
	uint8_t RANDOM_BYTES[32];
	uint64_t value = 0;

	// Spinlock for parallel execution without locking with
	// a mutex lock.
	this->mImpl->spinlock_.lock();
	while(true){
		RAND_bytes(&RANDOM_BYTES[0], 32);
		value = XXH64(&RANDOM_BYTES[0], 32, 1337);

		if(value == 0) continue;

		htable_type::const_iterator it = this->htable_.find(value);
		if(it == this->htable_.end()){
			if(store)
				this->htable_[value] = this->size();
			break;
		}
	}
	this->mImpl->spinlock_.unlock();

	return value;
}

bool Keychain::GetHashIdentifier(const uint64_t& value, uint32_t& match){
	htable_type::const_iterator it = this->htable_.find(value);
	if(it != this->htable_.end()){
		match = it->second;
		return(true);
	}

	return false;
}

bool Keychain::BuildHashTable(void){
	this->htable_.clear();

	for(uint32_t i = 0; i < this->size(); ++i)
		this->htable_[this->at(i)->field_id] = i;

	return true;
}

std::ostream& operator<<(std::ostream& stream, const Keychain& keychain){
	stream.write(TACHYON_MAGIC_HEADER.data(), TACHYON_MAGIC_HEADER_LENGTH);
	stream.write(reinterpret_cast<const char*>(&keychain.n_entries_),  sizeof(size_t));
	stream.write(reinterpret_cast<const char*>(&keychain.n_capacity_), sizeof(size_t));

	// Iteratively write keychain keys.
	for(uint32_t i = 0; i < keychain.size(); ++i)
		stream << *keychain[i];

	return(stream);
}

std::istream& operator>>(std::istream& stream, Keychain& keychain){
	char header[TACHYON_MAGIC_HEADER_LENGTH];
	stream.read(&header[0], TACHYON_MAGIC_HEADER_LENGTH);
	if(strncmp(&header[0], &TACHYON_MAGIC_HEADER[0], TACHYON_MAGIC_HEADER_LENGTH) != 0){
		std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Illegal keychain..." << std::endl;
		return(stream);
	}

	stream.read(reinterpret_cast<char*>(&keychain.n_entries_),  sizeof(size_t));
	stream.read(reinterpret_cast<char*>(&keychain.n_capacity_), sizeof(size_t));

	for(int i = 0; i < keychain.size(); ++i)
		delete keychain.entries_[i];
	delete [] keychain.entries_;

	keychain.entries_ = new Keychain::value_type[keychain.n_capacity_];
	for(uint32_t i = 0; i < keychain.size(); ++i){
		uint8_t type = 0;
		stream.read(reinterpret_cast<char*>(&type), sizeof(uint8_t));
		if(type == YON_ENCRYPTION_AES_256_GCM)
			keychain[i] = new KeychainKeyGCM<>;
		else {
			std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Unknown encryption type!" << std::endl;
			exit(1);
		}
		stream >> *keychain[i];
	}

	keychain.BuildHashTable();

	return(stream);
}

// decorator impl
class EncryptionDecorator::EncryptionDecoratorImpl {
public:
	EncryptionDecoratorImpl() = default;
	~EncryptionDecoratorImpl() = default;

    // Pimpl
	bool EncryptAES256(variant_block_type& block, keychain_type& keychain);
	bool EncryptAES256(stream_container& container, keychain_type& keychain);
	bool DecryptAES256(stream_container& container, keychain_type& keychain);

public:
	buffer_type buffer;
};

//decorator
bool EncryptionDecorator::Encrypt(variant_block_type& block, keychain_type& keychain, TACHYON_ENCRYPTION encryption_type){
	if(encryption_type == YON_ENCRYPTION_NONE){
		return true;
	} else if(encryption_type == YON_ENCRYPTION_AES_256_GCM){
		return(this->mImpl->EncryptAES256(block, keychain));
	} else {
		std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Not implemented yet!" << std::endl;
	}
	return(false);
}

bool EncryptionDecorator::Decrypt(variant_block_type& block, keychain_type& keychain){
	uint32_t temp_match = 0;

	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i){
		if(block.base_containers[i].IsEncrypted() == false) continue;
		if(block.base_containers[i].header.data_header.controller.encryption == YON_ENCRYPTION_AES_256_GCM){
			if(keychain.GetHashIdentifier(block.base_containers[i].header.identifier, temp_match)){
				this->mImpl->DecryptAES256(block.base_containers[i], keychain);
			}

		} else {
			std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Not implemented yet!" << std::endl;
		}
	}

	for(uint32_t i = 0; i < block.footer.n_info_streams; ++i){
		if(block.info_containers[i].IsEncrypted() == false) continue;
		if(block.info_containers[i].header.data_header.controller.encryption == YON_ENCRYPTION_AES_256_GCM){
			if(keychain.GetHashIdentifier(block.info_containers[i].header.identifier, temp_match)){
				this->mImpl->DecryptAES256(block.info_containers[i], keychain);
			}
		} else {
			std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Not implemented yet!" << std::endl;
		}
	}

	for(uint32_t i = 0; i < block.footer.n_format_streams; ++i){
		if(block.format_containers[i].IsEncrypted() == false) continue;
		if(block.format_containers[i].header.data_header.controller.encryption == YON_ENCRYPTION_AES_256_GCM){
			if(keychain.GetHashIdentifier(block.format_containers[i].header.identifier, temp_match)){
				this->mImpl->DecryptAES256(block.format_containers[i], keychain);
			}
		} else {
			std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Not implemented yet!" << std::endl;
		}
	}

	return true;
}

bool EncryptionDecorator::EncryptionDecoratorImpl::EncryptAES256(variant_block_type& block, keychain_type& keychain){
	block.header.block_hash = keychain.GetRandomHashIdentifier(true);

	// Iterate over available basic containers.
	for(uint32_t i = 1; i < YON_BLK_N_STATIC; ++i){
		if(!this->EncryptAES256(block.base_containers[i], keychain)){
			std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl;
			return false;
		}
	}

	// Iterate over Info containers.
	for(uint32_t i = 0; i < block.footer.n_info_streams; ++i){
		if(!this->EncryptAES256(block.info_containers[i], keychain)){
			std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl;
			return false;
		}
	}

	// Iterate over Format containers.
	for(uint32_t i = 0; i < block.footer.n_format_streams; ++i){
		if(!this->EncryptAES256(block.format_containers[i], keychain)){
			std::cerr << utility::timestamp("ERROR","ENCRYPTION") << "Failed to encrypt!" << std::endl;
			return false;
		}
	}

	return(true);
}

bool EncryptionDecorator::EncryptionDecoratorImpl::EncryptAES256(stream_container& container, keychain_type& keychain){
	KeychainKeyGCM<> entry;
	entry.InitiateRandom();
	entry.encryption_type = YON_ENCRYPTION_AES_256_GCM;

	EVP_CIPHER_CTX *ctx = NULL;
	int32_t len = 0;

	if(!(ctx = EVP_CIPHER_CTX_new())){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption context..." << std::endl;
		return false;
	}

	if(1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL)){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption..." << std::endl;
		return false;
	}

	// 16 uint8_ts = 128 bits
	if(1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, 16, NULL)){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption tag..." << std::endl;
		return false;
	}

	if(1 != EVP_EncryptInit_ex(ctx, NULL, NULL, entry.key, entry.iv)){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption key and IV..." << std::endl;
		return false;
	}

	// Mask header in encrypted message
	this->buffer.reset();
	io::BasicBuffer temp(65536);
	temp << container.header;
	this->buffer.resize(container.data.size() + container.strides.size() + temp.size() + 65536);

	if(1 != EVP_EncryptUpdate(ctx, (uint8_t*)this->buffer.data(), &len, (uint8_t*)temp.data(), temp.size())){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
		return false;
	}
	this->buffer.n_chars_ = len;

	if(1 != EVP_EncryptUpdate(ctx, (uint8_t*)&this->buffer[this->buffer.size()], &len, (uint8_t*)container.data.data(), container.data.size())){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
		return false;
	}
	this->buffer.n_chars_ += len;

	if(1 != EVP_EncryptUpdate(ctx,
		                      (uint8_t*)&this->buffer[this->buffer.size()],
	                          &len,
		                      (uint8_t*)container.strides.data(),
		                      container.strides.size()))
	{
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the encryption model..." << std::endl;
		return false;
	}
	this->buffer.n_chars_ += len;

	if(1 != EVP_EncryptFinal_ex(ctx,
	                            (uint8_t*)this->buffer.data() + len,
		                        &len))
	{
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to finalise the encryption..." << std::endl;
		return false;
	}
	this->buffer.n_chars_ += len;

	container.data.resize(this->buffer.size());
	if(1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_GET_TAG, 16, &entry.tag[0])){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to retrieve the GCM tag..." << std::endl;
		return false;
	}

	EVP_CIPHER_CTX_free(ctx);

	// Trigger encryption flag
	container.header.reset(); // reset data
	container.header.data_header.controller.encryption = YON_ENCRYPTION_AES_256_GCM;
	memcpy(container.data.data(), this->buffer.data(), this->buffer.size());
	container.data.n_chars_ = this->buffer.size();
	container.header.data_header.eLength = this->buffer.size();

	const uint64_t hashID = keychain.GetRandomHashIdentifier();
	container.header.identifier = hashID;
	entry.field_id = hashID;
	keychain += entry; // add key to keychain

	return(true);
}

bool EncryptionDecorator::EncryptionDecoratorImpl::DecryptAES256(stream_container& container, keychain_type& keychain){
	if(container.data.size() == 0)
		return true;

	if(container.header.data_header.controller.encryption == YON_ENCRYPTION_NONE)
		return true;

	if(container.header.data_header.controller.encryption != YON_ENCRYPTION_AES_256_GCM){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Wrong decryption algorithm used..." << std::endl;
		return false;
	}

	uint32_t match = 0;
	if(keychain.GetHashIdentifier(container.header.identifier, match) == false){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Did not find ID in keychain..." << std::endl;
		return false;
	}
	KeychainKeyGCM<>& entry = *reinterpret_cast<KeychainKeyGCM<>*>(keychain[match]);

	EVP_CIPHER_CTX *ctx = NULL;
	int32_t len = 0, plaintext_len = 0, ret;

	/* Create and initialise the context */
	if(!(ctx = EVP_CIPHER_CTX_new())){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the encryption context..." << std::endl;
		return(false);
	}

	/* Initialise the decryption operation. */
	if(!EVP_DecryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL)){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the decryption..." << std::endl;
		return(false);
	}

	if(!EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, 16, NULL)){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the IV context..." << std::endl;
		return(false);
	}

	/* Initialise key and IV */
	if(!EVP_DecryptInit_ex(ctx, NULL, NULL, entry.key, entry.iv)){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to initialise the decryption..." << std::endl;
		return(false);
	}

	this->buffer.reset();
	if(container.data.size()){
		this->buffer.resize(container.data.size() + 65536);
		if(!EVP_DecryptUpdate(ctx, (uint8_t*)this->buffer.data(), &len, (uint8_t*)container.data.data(), container.data.size())){
			std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to update the decryption..." << std::endl;
			return(false);
		}
		plaintext_len = len;
	}

	/* Set expected tag value. Works in OpenSSL 1.0.1d and later */
	if(!EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_TAG, 16, entry.tag)){
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to set expected TAG value..." << std::endl;
		return(false);
	}

	/* Finalise the decryption. A positive return value indicates success,
	 * anything else is a failure - the plaintext is not trustworthy.
	 */
	ret = EVP_DecryptFinal_ex(ctx, (uint8_t*)this->buffer.data() + len, &len);
	EVP_CIPHER_CTX_free(ctx);

	if(ret > 0){
		plaintext_len += len;
		this->buffer >> container.header; // unmask encrypted header

		container.data.reset();
		container.strides.reset();
		container.data.resize(container.header.data_header.cLength + 65536);
		container.strides.resize(container.header.stride_header.cLength + 65536);
		memcpy(container.data.data(), &this->buffer[this->buffer.iterator_position_], container.header.data_header.cLength);
		memcpy(container.strides.data(), &this->buffer[this->buffer.iterator_position_ + container.header.data_header.cLength], container.header.stride_header.cLength);
		container.header.data_header.controller.encryption = YON_ENCRYPTION_NONE;
		container.data.n_chars_    = container.header.data_header.cLength;
		container.strides.n_chars_ = container.header.stride_header.cLength;
		container.header.data_header.eLength = 0;
		return(true);
	} else {
		std::cerr << utility::timestamp("ERROR", "ENCRYPTION") << "Failed to validate decryption..." << std::endl;
		return(false);
	}
}

// Ctor for decorator
EncryptionDecorator::EncryptionDecorator() : mImpl(new EncryptionDecorator::EncryptionDecoratorImpl){}
EncryptionDecorator::~EncryptionDecorator(){}

}
