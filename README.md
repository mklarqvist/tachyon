[![Build Status](https://travis-ci.org/mklarqvist/tachyon.svg?branch=master)](https://travis-ci.org/mklarqvist/tachyon)
[![Release](https://img.shields.io/badge/Release-beta_0.1-blue.svg)](https://github.com/mklarqvist/Tachyon/releases)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

<div align="center">
<img src="https://github.com/mklarqvist/tachyon/blob/master/yon_logo.png"><br><br>
</div>

# Exploring population-scale sequence variant data
Tachyon, or `YON` for short, is an open source C++ software library for reading, writing, and manipulating sequence variant data in a lossless and bit-exact representation. It is completely compatible with BCF/VCF. It was developed with a focus on enabling fast experimentation and storage of population-scaled datasets.

## Highlights of Tachyon
* **Self-indexing**: Tachyon always builds the best possible quad-tree, linear, and meta-index given the input data (irrespective of sorting). There are no external indices as data are stored in the file itself.
* **Integrity checking**: The `YON` specification enforces validity checks for each data field and across all fields through checksum validation. This approach guarantees file integrity when compressing/decompressing and encrypting/decrypting. These checksums are stored internally.
* **Encryption**: Natively supports block-wise, field-wise, and entry-wise encryption with all commonly used encryption models and paradigms through [openssl][openssl].
* **Compression**: Tachyon files are generally many fold (in many cases many 10- to 100-folds) smaller than the current standard file-format.
* **Field-specific layout**: In principle, Tachyon is implemented as a standard column-oriented management system with several layers of domain-specific heuristics providing fast and flexible data queries. This memory layout enables extremely rapid field-specific queries.  
* **High-level API**: User-friendly C++/C API for quering, manipulating, and exploring sequence data with minimal programming experience
* **Comaptibility**: We strive to provide API calls to return YON data streams to any of the current standard file-formats (`VCF` and `BCF`). This allows for immediate use of Tachyon without disrupting the existing ecosystem of tools.

---

## Project status
Tachyon is under active development and the specification and/or the API interfaces may change at any time!   
**Commits may break functionality! THERE IS NO STABILITY PROMISE WHATSOEVER!**  

Current limitations imposed during development:
* Importing is restricted to `BCF`
* Output is restricted to `VCF`, `JSON`, and custom field slicing

## Installation
### Dependencies
You will need to have installed the following dependencies:
* [zstd][zstd]: A compression library developed at Facebook
* [openssl][openssl]: An open-source library for encryption/decryption

### Building from source
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

### Building without admin privilidges
If you have no super-user (`sudo`) powers required to install software on your machine:

### Linux/MacOSX
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

## Documentation

* Overview.
* [Getting started](docs/getting_started.md)
* [Summary of example programs](docs/example_programs.md).
* Python API Reference.


### Contributing

Interested in contributing? See CONTRIBUTING.

### Support

The Genomics team in Google Brain actively supports Nucleus and are always interested in improving its quality. If you run into an issue, please report the problem on our Issue tracker. Be sure to add enough detail to your report that we can reproduce the problem and fix it. We encourage including links to snippets of BAM/VCF/etc files that provoke the bug, if possible. Depending on the severity of the issue we may patch Nucleus immediately with the fix or roll it into the next release.

### Version
This is Tachyon 0.1.0.Tachyon follows [semantic versioning](https://semver.org/).

### History
Tachyon grew out of the [Tomahawk][tomahawk] project for calculating genome-wide linkage-disequilibrium.

### Author
Marcus D. R. Klarqvist (<mk819@cam.ac.uk>)  
Department of Genetics, University of Cambridge  
Wellcome Trust Sanger Institute

### Acknowledgements
[James Bonfield](https://github.com/jkbonfield), Wellcome Trust Sanger Institute  
[Petr Daněček](https://github.com/pd3), Wellcome Trust Sanger Institute  
[Richard Durbin](https://github.com/richarddurbin), Wellcome Trust Sanger Institute, and Department of Genetics, University of Cambridge  

### License
Tachyon is licensed under [MIT](LICENSE)

[openssl]:  https://www.openssl.org/
[zstd]:     https://github.com/facebook/zstd
[tomahawk]: https://github.com/mklarqvist/tomahawk
[msprime]:  https://github.com/jeromekelleher/msprime