[![Release](https://img.shields.io/badge/Release-beta_0.1-blue.svg)](https://github.com/mklarqvist/Tachyon/releases)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)


# Tachyon
<div align="center">
<img src="https://github.com/mklarqvist/Tachyon/blob/master/yon_logo.png"><br><br>
</div>

Tachyon is an open source software library for storing and querying variant data. Tachyon efficiently stores data fields by column and implicitly represent genotypic data by exploiting intrinsic genetic properties. Most genotype-specific algorithms were originally developed for [Tomahawk][tomahawk] for the purpose of calculating linkage-disequilibrium in large-scale cohorts.

Marcus D. R. Klarqvist (<mk819@cam.ac.uk>)  
Department of Genetics, University of Cambridge  
Wellcome Trust Sanger Institute

### Notice
Tachyon is under active development and the specification and/or the API interfaces may change at any time! Commits may break functionality!  
THERE IS NO API STABILITY PROMISE WHATSOEVER!

### Installation instructions
Building requires [zstd][zstd] and [openssl][openssl]
```bash
git clone --recursive https://github.com/mklarqvist/Tachyon
cd Tachyon/build
make
```

[openssl]: https://www.openssl.org/
[zstd]: https://github.com/facebook/zstd
[tomahawk]: https://github.com/mklarqvist/tomahawk

### License
[MIT](LICENSE)
