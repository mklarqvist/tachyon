#ifndef ALGORITHM_ENCRYPTION_KEYCHAIN_KEY_H_
#define ALGORITHM_ENCRYPTION_KEYCHAIN_KEY_H_

namespace tachyon{
namespace encryption{

template <BYTE KeyLength = 32, BYTE IVLength = 16>
struct KeychainKey{
public:
	typedef KeychainKey     self_type;
	typedef io::BasicBuffer buffer_type;

public:
	KeychainKey() : field_id(0), encryption_type(YON_ENCRYPTION_NONE){}
	virtual ~KeychainKey(){}

	KeychainKey(const KeychainKey& other) :
		field_id(other.field_id),
		encryption_type(other.encryption_type)
	{
		memcpy(&this->key[0], &other.key[0], KeyLength);
		memcpy(&this->iv[0],  &other.iv[0],  IVLength);
	}

	KeychainKey& operator=(const KeychainKey& other){
		this->field_id = other.field_id;
		this->encryption_type = other.encryption_type;
		memcpy(&this->key[0], &other.key[0], KeyLength);
		memcpy(&this->iv[0],  &other.iv[0],  IVLength);
		return(*this);
	}

	virtual void print(void){
		for(U32 i = 0; i < KeyLength; ++i) std::cerr << std::hex << (int)this->key[i];
		std::cerr << '\t';
		for(U32 i = 0; i < IVLength; ++i) std::cerr << std::hex << (int)this->iv[i];
	}

	inline bool initiateRandom(void){
		RAND_bytes(&this->key[0], KeyLength);
		RAND_bytes(&this->iv[0],  IVLength);
		return true;
	}

private:
	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const self_type& key){
		buffer += key.field_id;
		buffer += key.encryption_type;
		buffer.Add((const char*)&key.key[0], KeyLength);
		buffer.Add((const char*)&key.iv[0],  IVLength);
		return(buffer);
	}

public:
	U64  field_id;
	BYTE encryption_type;
	BYTE key[KeyLength];  // 256 bit key
	BYTE iv[IVLength];   // 128 bit initiation vector
};

template <BYTE KeyLength = 32, BYTE IVLength = 16, BYTE TagLength = 16>
struct KeychainKeyGCM : public KeychainKey<KeyLength, IVLength>{
public:
	typedef KeychainKeyGCM                   self_type;
	typedef KeychainKey<KeyLength, IVLength> parent_type;

public:
	KeychainKeyGCM(){}
	~KeychainKeyGCM(){}

	KeychainKeyGCM(const KeychainKeyGCM& other) :
		parent_type(other)
	{
		memcpy(&this->tag[0], &other.tag[0], TagLength);
	}

	KeychainKeyGCM& operator=(const KeychainKeyGCM& other){
		this->field_id = other.field_id;
		this->encryption_type = other.encryption_type;
		memcpy(&this->key[0], &other.key[0], KeyLength);
		memcpy(&this->iv[0],  &other.iv[0],  IVLength);
		memcpy(&this->tag[0], &other.tag[0], TagLength);
		return(*this);
	}

	void print(void){
		for(U32 i = 0; i < KeyLength; ++i) std::cerr << std::hex << (int)this->key[i];
		std::cerr << '\t';
		for(U32 i = 0; i < IVLength; ++i) std::cerr << std::hex << (int)this->iv[i];
		std::cerr << '\t';
		for(U32 i = 0; i < TagLength; ++i) std::cerr << std::hex << (int)this->tag[i];
		std::cerr << std::dec;
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& key){
		stream.write(reinterpret_cast<const char*>(&key.field_id), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&key.encryption_type), sizeof(BYTE));
		stream.write((char*)&key.key[0], KeyLength);
		stream.write((char*)&key.iv[0],  IVLength);
		stream.write((char*)&key.tag[0], TagLength);
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& key){
		stream.read(reinterpret_cast<char*>(&key.field_id), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&key.encryption_type), sizeof(BYTE));
		stream.read((char*)&key.key[0], KeyLength);
		stream.read((char*)&key.iv[0],  IVLength);
		stream.read((char*)&key.tag[0], TagLength);
		return(stream);
	}

	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const self_type& key){
		buffer += key.field_id;
		buffer += key.encryption_type;
		buffer.Add((const char*)&key.key[0], KeyLength);
		buffer.Add((const char*)&key.iv[0],  IVLength);
		buffer.Add((const char*)&key.tag[0], TagLength);
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& key){
		buffer >> key.field_id;
		buffer >> key.encryption_type;
		buffer.read((char*)&key.key[0], KeyLength);
		buffer.read((char*)&key.iv[0],  IVLength);
		buffer.read((char*)&key.tag[0], TagLength);
		return(buffer);
	}

public:
	BYTE tag[TagLength];  // validation tag for gcm
};

}
}



#endif /* ALGORITHM_ENCRYPTION_KEYCHAIN_KEY_H_ */
