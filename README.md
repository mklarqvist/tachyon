[![Release](https://img.shields.io/badge/Release-beta_0.2-blue.svg)](https://github.com/mklarqvist/Tachyon/releases)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)


# Tachyon
<div align="center">
<img src="https://github.com/mklarqvist/tachyon/blob/master/yon_logo.png"><br><br>
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
git clone --recursive https://github.com/mklarqvist/tachyon
cd tachyon/build
make
```

### ABI examples
Import a `bcf` file to `yon` with a block-size of `-c` number of variants and/or `-C` number of base-pairs. If both `-c` and `-C` are set then tachyon breaks whenever either condition is satisfied.
```bash
tachyon import -i <file.bcf> -o <outfile.yon> -c <variants checkpoint> -C <base pair checkpoint>
```

Viewing a `yon` file
```bash
tachyon view -i <input.yon>
```

### C++ API Examples
```c++
/**<
 * Tachyon: in this example we will load data from
 * the FORMAT field GL into the iterable template class
 * `FormatContainer`. This is a complete example!
 */
#include <tachyon/tachyon_reader.h>

std::string my_input_file = "somefile.yon"; // Change me to an actual file that exists on your filesystem
tachyon::TachyonReader reader;
reader.getSettings().loadFormat("GL");
reader.open(my_input_file);

/**<
 *  The `FormatContainer` class stores the data for each variant 
 *  for each individual as container[variant][sample][data]
 */
while(reader.get_next_block()){ // As long as there are YON blocks available
    // Meta container
    containers::MetaContainer meta(this->block);
    // FORMAT container with float return type primitive
    containers::FormatContainer<float>* gp_container = this->get_balanced_format_container<float>("GL", meta);
    if(gp_container != nullptr){
        for(U32 variant = 0; variant < gp_container->size(); ++variant){
            for(U32 sample = 0; sample < gp_container->at(variant).size(); ++sample){
                // Write the data to `cout` in `VCF` formatting
                utility::to_vcf_string(std::cout, gp_container->at(variant).at(sample)) << ' ';
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
#include <tachyon/tachyon_reader.h>

std::string my_input_file = "somefile.yon"; // Change me to an actual file that exists on your filesystem
tachyon::TachyonReader reader;
reader.getSettings().loadInfo("SVLEN");
reader.open(my_input_file);

/**<
 * The `InfoContainer` class stores the data for each variant as
 * container[variant][data]. Both `InfoContainer` and `FormatContainer`
 * supports variant-balancing of the classes. Balancing refers to filling
 * variant sites in the file with empty objects if no target data is present
 * at that site. 
 */
while(reader.get_next_block()){ // As long as there are YON blocks available
    // Meta container
    containers::MetaContainer meta(this->block);
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

More advanced example using genotype summary statistics
```c++
/**<
* In this example we will use the genotype summary statistics
* to calculate allele-specific bias using a Fisher's 2x2 exact
* test
*/
#include <tachyon/tachyon_reader.h>

std::string my_input_file = "somefile.yon"; // Change me to an actual file that exists on your filesystem
tachyon::TachyonReader reader;
reader.getSettings().loadFormat("GL");
reader.open(my_input_file);

while(reader.get_next_block()){ // As long as there are YON blocks available
    containers::GenotypeContainer gt(this->block);
    math::Fisher fisher(); // Math class for Fisher's exact test and Chi-squared
    containers::GenotypeSum gt_summary; // Genotype summary statistics
    for(U32 i = 0; i < gt.size(); ++i){ // Foreach variant
        // If there's > 5 alleles continue
        if(gt[i].getMeta().getNumberAlleles() >= 5) continue;
        // Calculate summary statistics
        gt[i].getSummary(gt_summary);

        // Calculate total number of alt-alleles (allele 1, where 0 is ref)
        const U64 total = gt_summary.getAlleleA(1) + gt_summary.getAlleleB(1);
        const double p = fisher.fisherTest(gt_summary.getAlleleA(1), total, gt_summary.getAlleleB(1), total); // P-value for allele-bias
        if(p < 1e-3){ // If P < 0.001 report it
            gt[i].getMeta().toVCFString(stream, this->header, this->block.index_entry.contigID, this->block.index_entry.minPosition);
            std::cout << '\t' << gt_summary << '\t' << p << '\t' << ((gt_summary.getAlleleA(1) == 0 || gt_summary.getAlleleB(1) == 0) ? 1 : 0) << '\n';
        }
        gt_summary.clear(); // Recycle summary object
    }
}
```

```c++
/**<
* The higher-order primitive GTObject allows for single-genotype
* manipulation whenever summary statistics is insufficient. There
* are three levels of return types:
* 1) The literal internal representation (preffered use). This
*    function returns the implementation representation as a generic
*    GTObject. 
* 2) Unpacked (one genotype -> one sample) GTObject but in permuted
*    order. These genotypes are in the internal sorted order and does
*    not match with the tachyon sample header
* 3) Unpacked (one genotype -> one sample) GTObject in original order.
*    These genotype objects are returned in the same order as described
*    in the tachyon sample header. 
*/
#include <tachyon/tachyon_reader.h>

std::string my_input_file = "somefile.yon"; // Change me to an actual file that exists on your filesystem
tachyon::TachyonReader reader;
reader.getSettings().loadFormat("GL");
reader.open(my_input_file);

while(reader.get_next_block()){ // As long as there are YON blocks available
    containers::GenotypeContainer gt(this->block);
    for(U32 i = 0; i < gt.size(); ++i){
        // All of these functions are in relative terms very expensive!
        // Avoid using them unless you have to!
        // Vector of literal genotype representations (lower level)
        std::vector<core::GTObject> objects     = gt[i].getLiteralObjects();
        // Vector of genotype objects (high level permuted)
        std::vector<core::GTObject> objects_all = gt[i].getObjects(this->header.n_samples);
        // Vector of genotype objects (high level unpermuted - original)
        std::vector<core::GTObject> objects_true = gt[i].getObjects(this->header.n_samples, this->block.ppa_manager);

        // Print the difference
        std::cerr << objects.size() << '\t' << objects_all.size() << '\t' << objects_true.size() << std::endl;
        // Dump data
        gt[i].getMeta().toVCFString(stream, this->header, this->block.index_entry.contigID, this->block.index_entry.minPosition);
        for(U32 j = 0; j < objects_true.size(); ++j){
            std::cout << (int)objects_true[j].alleles[0].first << (objects_true[j].alleles[1].second ? "|" : "/") << (int)objects_true[j].alleles[1].first << ' ';
        }
        std::cout << '\n';
    }
}
```

[openssl]:  https://www.openssl.org/
[zstd]:     https://github.com/facebook/zstd
[tomahawk]: https://github.com/mklarqvist/tomahawk

### License
[MIT](LICENSE)
