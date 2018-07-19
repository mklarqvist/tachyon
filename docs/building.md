# Installing Tachyon
## Dependencies
You will need to have installed the following dependencies:
* [zstd][zstd]: A compression library developed at Facebook
* [openssl][openssl]: An open-source library for encryption/decryption

## Building from source
If the required external dependencies listed above are installed then building is trivial. Note the added `--recursive` flag to the clone request. This flag is required to additionally pull down the latest third-party dependencies.
```bash
git clone --recursive https://github.com/mklarqvist/tachyon
cd tachyon
make
```
Tachyon comes bundled with several API-examples in the `lib_example` directory. They are built by default but should you want to rebuild them execute the command:
```bash
make examples
```

## Building without admin privilidges
If you have no super-user (`sudo`) powers required to install software on your machine:

## Linux/MacOSX
```bash
git clone --recursive https://github.com/mklarqvist/tachyon
cd tachyon
# If you do NOT have ZSTD available
git clone https://github.com/facebook/zstd
cd zstd
make
cd ..
# If you do NOT have OpenSSL installed
git clone git://git.openssl.org/openssl.git
cd openssl
./config
make
cd ..
# Build Tachyon
make
```
### MacOSX
Installation using [Homebrew](https://brew.sh/):
```bash
brew update
# If you do NOT have OpenSSL installed
brew install openssl
# If you do NOT have ZSTD installed
brew install zstd
# Install Tachyon
git clone --recursive https://github.com/mklarqvist/tachyon
cd tachyon
make
```