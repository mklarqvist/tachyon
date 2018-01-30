[![Release](https://img.shields.io/badge/Release-beta_0.2-blue.svg)](https://github.com/mklarqvist/Tachyon/releases)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)


# Tachyon
<div align="center">
<img src="https://github.com/mklarqvist/Tachyon/blob/master/yon_logo.png"><br><br>
</div>

Tachyon is an open source software library for storing and querying variant data. Tachyon efficiently stores data fields by column and implicitly represent genotypic data by exploiting intrinsic genetic properties. Most genotype-specific algorithms were originally developed for [Tomahawk][tomahawk] for the purpose of calculating linkage-disequilibrium in large-scale cohorts.

The Tachyon specification has complete backward-compatibility with `VCF`/`BCF` and is many 10s- to 100s-fold smaller compared to `BCF`. The library supports additional functionality like field-specific and granular encryption (symmetric, assymetric, and eventually homomorphic).

### Author
Marcus D. R. Klarqvist (<mk819@cam.ac.uk>)  
Department of Genetics, University of Cambridge  
Wellcome Trust Sanger Institute

### Notice
Tachyon is under active development and the specification and/or the API interfaces may change at any time! Commits may break functionality!  
THERE IS NO API STABILITY PROMISE WHATSOEVER!  
Documentation is currently sparse.

### Installation instructions
Building requires [zstd][zstd] and [openssl][openssl]
```bash
git clone --recursive https://github.com/mklarqvist/Tachyon
cd Tachyon/build
make
```

### C++ API Examples
```c++
/**<
 * Tachyon: in this example we will load data from
 * the FORMAT field GL into the iterable template class
 * `FormatContainer`. This is a complete example!
 */
#include <tachyon/tachyon.h>

std::string my_input_file = "somefile.yon"; // Change me to an actual file that exists on your filesystem
tachyon::TachyonReader reader;
reader.getSettings().loadFormat("GL");
reader.open(my_input_file);

/**<
 *  The `FormatContainer` class stores the data for each variant 
 *  for each individual as container[variant][sample][data]
 */
while(reader.get_next_block()){ // As long as there is YON blocks available
    containers::FormatContainer<float>* gp_container = this->get_balanced_format_container<float>("GL", meta);
    if(gp_container != nullptr){
        for(U32 variant = 0; variant < gp_container->size(); ++variant){
            for(U32 sample = 0; sample < gp_container->at(variant).size(); ++sample){
                // Write the data to `cout` in `VCF` formatting
                util::to_vcf_string(std::cout, gp_container->at(variant).at(sample)) << ' ';
            }
            std::cout << '\n';
        }
        std::cout << '\n';
    }
    delete gp_container;
}
```

```c++
/**<
 * Tachyon: in this example we will load data from
 * the INFO field SVLEN into the iterable template class
 * `InfoContainer`. This is a complete example!
 */
#include <tachyon/tachyon.h>

std::string my_input_file = "somefile.yon"; // Change me to an actual file that exists on your filesystem
tachyon::TachyonReader reader;
reader.getSettings().loadFormat("GL");
reader.open(my_input_file);

/**<
 * The `InfoContainer` class stores the data for each variant as
 * container[variant][data]. Both `InfoContainer` and `FormatContainer`
 * supports variant-balancing of the classes. Balancing refers to filling
 * variant sites in the file with empty objects if no tdaargetta is present
 * at that site. 
 */
while(reader.get_next_block()){
    // Variant-balanced
	containers::InfoContainer<U32>* info_balanced   = this->get_balanced_info_container<U32>("SVLEN", meta);
	// Not variant-balanced
	containers::InfoContainer<U32>* info_unbalanced = this->get_info_container<U32>("SVLEN");
    // Print the sizes of the two containers
    // The size of the balanced container is always the number of variants
    // In contrast, the unbalanced one returns a container only for the 
    // variants with data available. 
    std::cout << info_balanced->size() << "/" << info_unbalanced->size() << '\n';
    delete info_balanced;
    delete info_unbalanced;
}
```

[openssl]:  https://www.openssl.org/
[zstd]:     https://github.com/facebook/zstd
[tomahawk]: https://github.com/mklarqvist/tomahawk

### License
[MIT](LICENSE)
