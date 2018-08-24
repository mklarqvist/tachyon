[![Build Status](https://travis-ci.org/mklarqvist/tachyon.svg?branch=master)](https://travis-ci.org/mklarqvist/tachyon)
[![Release](https://img.shields.io/badge/Release-beta_0.3.0-blue.svg)](https://github.com/mklarqvist/Tachyon/releases)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

<div align="center">
<img src="https://github.com/mklarqvist/tachyon/blob/master/yon_logo.png"><br><br>
</div>

Tachyon, or `YON` for short, is an open source C++ software library for reading, writing, and manipulating sequence variant data in a lossless and bit-exact representation. It is completely compatible with BCF/VCF. It was developed with a focus on enabling fast experimentation and storage of population-scaled datasets.

## Highlights of Tachyon
* **Self-indexing**: Tachyon always builds the best possible quad-tree, linear, and meta-index given the input data (irrespective of sorting). There are no external indices as data are stored in the file itself.
* **Integrity checking**: The `YON` specification enforces validity checks for each data field and across all fields through checksum validation. This approach guarantees file integrity when compressing/decompressing and encrypting/decrypting. These checksums are stored internally.
* **Encryption**: Natively supports block-wise, field-wise, and entry-wise encryption with all commonly used encryption models and paradigms through [openssl][openssl].
* **Compression**: Tachyon files are generally many fold (in many cases many 10- to 100-folds) smaller than the current standard file-format.
* **Field-specific layout**: In principle, Tachyon is implemented as a standard column-oriented management system with several layers of domain-specific heuristics providing fast and flexible data queries. This memory layout enables extremely rapid field-specific queries.  
* **Performance**: The file-format is designed as independent blocks of data into independent byte streams. This approach is inherently amenable to paralellization through scatter-gather approaches on multiple cores or multiple machines.
* **Comaptibility**: We strive to provide API calls to return YON data streams to any of the current standard file-formats (`VCF`, `VCF.GZ`, and `BCF`). This allows for immediate use of Tachyon without disrupting the existing ecosystem of tools.

---  

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

## Documentation

* Overview.
* [Building and installing](docs/building.md)
* [Getting started](docs/getting_started.md)
* [Summary of example programs](docs/example_programs.md).
* [Performance benchmarks](docs/benchmarks.md)

### Contributing

Interested in contributing? Fork and submit a pull request and it will be reviewed.

### Support
We are actively developing Tachyon and are always interested in improving its quality. If you run into an issue, please report the problem on our Issue tracker. Be sure to add enough detail to your report that we can reproduce the problem and address it. We have not reached version 1.0 and as such the specification and/or the API interfaces may change.

### Version
This is Tachyon 0.3.0. Tachyon follows [semantic versioning](https://semver.org/).

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