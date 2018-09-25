# Installing Tachyon
## Dependencies
You will need to have installed the following dependencies:
* [zstd][zstd]: A compression library developed at Facebook
* [openssl][openssl]: An open-source library for encryption/decryption
* [htslib][htslib]: C library for high-throughput sequencing data formats 

## Building from source
If the required external dependencies listed above are installed then building is trivial. Note the added `--recursive` flag to the clone request. This flag is required to additionally pull down the latest third-party dependencies.
```bash
git clone --recursive https://github.com/mklarqvist/tachyon
cd tachyon
make
```

## Installation
For Ubuntu, Debian, and Mac systems, installation is easy: just run
```bash
git clone --recursive https://github.com/mklarqvist/tachyon
cd tachyon
sudo ./install.sh
```
Note the added `--recursive` flag to the clone request. This flag is required to additionally pull down the latest third-party dependencies. The install.sh file depends extensively on apt-get, so it is unlikely to run without extensive modifications on non-Debian-based systems.
If you do not have super-user privileges required to install new packages on your system then run
```bash
./install.sh local
```
In this situation, all required dependencies are downloaded and built in the current directory. This approach will require additional effort if you intend to move the compiled libraries to a new directory.

[openssl]:  https://www.openssl.org/
[zstd]:     https://github.com/facebook/zstd
[tomahawk]: https://github.com/mklarqvist/tomahawk
[htslib]:   https://github.com/samtools/htslib