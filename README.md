[![Release](https://img.shields.io/badge/Release-beta_0.1-blue.svg)](https://github.com/mklarqvist/Tachyon/releases)
[![Build Status](https://travis-ci.org/mklarqvist/tachyon.svg?branch=master)](https://travis-ci.org/mklarqvist/tachyon)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)


# Tachyon
<div align="center">
<img src="https://github.com/mklarqvist/tachyon/blob/master/yon_logo.png"><br><br>
</div>

Tachyon is an open source software library for storing and querying sequencing/sequence variant data. Tachyon efficiently stores data fields by column and implicitly represent genotypic data by exploiting intrinsic genetic properties. Most genotype-specific algorithms were originally developed for [Tomahawk][tomahawk] for the purpose of calculating linkage-disequilibrium and identity-by-state in large-scale cohorts.

### Author
Marcus D. R. Klarqvist (<mk819@cam.ac.uk>)  
Department of Genetics, University of Cambridge  
Wellcome Trust Sanger Institute

## Notice!
Tachyon is under active development and the specification and/or the API interfaces may change at any time!   
Commits may break functionality!  
**THERE IS NO STABILITY PROMISE WHATSOEVER!**  

## Introduction

### Why a new framework?
There are a large number of field-specific file formats that have reached near-universal standard such as FASTA/FASTQ, SAM/BAM/uBAM/CRAM, VCF/BCF/gVCF. They all have unique file-specifications and different native toolsets interacting with them in addition to the multitude of tools developed for answering specific scientific questions. Tachyon was developed as a format-agnostic framework encapsulating all these previous standard formats in a unified database-like system constructed from generalized data-agnostic containers. The primary incentive of this project is to empower the research community with the tools required to query and interact with sequencing data.

---

### Guiding principles
* **Uniformity**. Tachyon is designed from data-agnostic STL-like containers and as such is decoupled from the higher-order file type-specific implementation details. This forces all underlying specifcations to share the same API calls.
* **User friendliness**. Tachyon is an API designed to be used by human beings and simultaneously a backbone specification consumed by machines. It puts user experience front and center. Tachyon follows best practices for reducing cognitive load: it offers consistent & simple APIs, it minimizes the number of user actions required for common use cases, and it provides clear and actionable feedback upon user error.
* **Modularity**. Tachyon is composed entirely of STL-like data containers that in turn are abstracted into higher-order type-specific containers. These containers are all standalone and can be plugged together with little restriction in any context.
* **Easy extensibility**. Describing novel specifcations and/or adding new functionality is simple as all basic containers are data-agnostic. Allowing for easy extensibility and to create new tools and modules allows for near-unlimited expressiveness making Tachyon suitable for advanced research.

---

### Highlights of Tachyon
* **Self-indexing**: Tachyon always builds the best possible index and super-index (indexing of the index for even faster queries) given the input data (irrespective of sorting). There are no external indices as data are stored in the file itself.
* **Integrity checking**: The `YON` specification enforces validity checks for each data field and across all fields through checksum validation. Guaranteed file integrity when compressing/decompressing and encrypting/decrypting.
* **Encryption**: Natively supports block-wise, field-wise, and entry-wise encryption with all commonly used encryption models and paradigms
* **Compression**: Tachyon files are generally many fold (in many cases many 10-100-fold) smaller than standard formats
* **Field-specific layout**: In principle, Tachyon is implemented as a standard column-oriented database management system with several layers of domain-specific heuristics providing  
* **High-level API**: User-friendly C++/C API for quering, manipulating, and exploring sequence data with minimal programming experience
* **Comaptibility**: We strive to provide API calls to return YON data streams to any of the file-formats we are supplanting. This allows for immediate use of Tachyon without disrupting the existing ecosystem of tools

---

## Getting started
### Dependencies
You will need the following dependencies:
* [zstd][zstd]: A compression library developed at Facebook
* [openssl][openssl]: An open-source library for encryption/decryption

### Building from source
Assuming the dependencies are installed then building is trivial:
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

## C++ API Examples
### Standard containers
```c++
/**<
 * Tachyon: https://github.com/mklarqvist/tachyon 
 * In this example we will load data from
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
while(reader.nextBlock()){ // As long as there are YON blocks available
    // Meta container
    containers::MetaContainer meta(reader.block);
    // FORMAT container with float return type primitive
    containers::FormatContainer<float>* gp_container = reader.get_balanced_format_container<float>("GL", meta);
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
 * Tachyon: https://github.com/mklarqvist/tachyon 
 * In this example we will load data from
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
while(reader.nextBlock()){ // As long as there are YON blocks available
    // Meta container
    containers::MetaContainer meta(reader.block);
    // Variant-balanced
    containers::InfoContainer<U32>* info_balanced   = reader.get_balanced_info_container<U32>("SVLEN", meta);
    // Not variant-balanced 
    containers::InfoContainer<U32>* info_unbalanced = reader.get_info_container<U32>("SVLEN");
    // Print the sizes of the two containers
    // The size of the balanced container is always the number of variants
    // In contrast, the unbalanced one returns a container only for the 
    // variants with data available. 
    std::cout << info_balanced->size() << "/" << info_unbalanced->size() << '\n';
    delete info_balanced;
    delete info_unbalanced;
}
```

```c++
/**<
 * Tachyon: https://github.com/mklarqvist/tachyon 
 * In this example we will load data from
 * the INFO field MEINFO into the iterable template class
 * `InfoContainer`. This is a complete example!
 */
#include <tachyon/tachyon_reader.h>

std::string my_input_file = "somefile.yon"; // Change me to an actual file that exists on your filesystem
tachyon::TachyonReader reader;
reader.getSettings().loadInfo("MEINFO");
reader.open(my_input_file);

/**<
 * The `InfoContainer` class has to be templated as std::string if the underlying data is of type `char`
 */
while(reader.nextBlock()){ // As long as there are YON blocks available
    // Meta container
    containers::MetaContainer meta(reader.block);
    containers::InfoContainer<std::string>* meinfo_container = reader.get_balanced_info_container<std::string>("MEINFO", meta);
    std::cout << meinfo_container->size() << std::endl;
    delete meinfo_container;
}
```

```c++
/**<
 * Tachyon: https://github.com/mklarqvist/tachyon 
 * In this example we will load data from
 * the VEP-generated CSQ string and tokenize it using
 * built-in utility functions
 */
#include <tachyon/tachyon_reader.h>

std::string my_input_file = "somefile.yon"; // Change me to an actual file that exists on your filesystem
tachyon::TachyonReader reader;
reader.getSettings().loadInfo("CSQ");
reader.open(my_input_file);

while(reader.nextBlock()){ // As long as there are YON blocks available
    containers::MetaContainer meta(reader.block);
    containers::InfoContainer<std::string>* csq_data = reader.get_balanced_info_container<U32>("std::string", meta);
    if(it2 != nullptr){
        for(size_t i = 0; i < csq_data->size(); ++i){
            if(csq_data->at(i).size() == 0) continue;
            else {
                meta.at(i).toVCFString(std::cout, this->header);
                std::cout.put('\t');
                // Tokenize CSQ string into a vector of strings
                std::vector<std::string> ret = utility::split((*csq_data)[i],'|');
                // Dump element 6 out of ret.size() total elements
                std::cout << ret[6] << "\n";
            }
        }
    }
    delete csq_data;
}
```

---

### Genotype containers / objects
More advanced example using genotype summary statistics
```c++
/**<
* Tachyon: https://github.com/mklarqvist/tachyon 
* In this example we will use the genotype summary statistics
* to calculate strand-specific bias of an alelle using a Fisher's 
* 2x2 exact test. This is a complete example!
*/
#include <tachyon/tachyon_reader.h>

std::string my_input_file = "somefile.yon"; // Change me to an actual file that exists on your filesystem
tachyon::TachyonReader reader;
reader.getSettings().loadGenotypes(true);
reader.open(my_input_file);

while(reader.nextBlock()){ // As long as there are YON blocks available
    containers::GenotypeContainer gt(reader.block);
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
            gt[i].getMeta().toVCFString(std::cout, reader.header, reader.block.index_entry.contigID, reader.block.index_entry.minPosition);
            std::cout << '\t' << gt_summary << '\t' << p << '\t' << ((gt_summary.getAlleleA(1) == 0 || gt_summary.getAlleleB(1) == 0) ? 1 : 0) << '\n';
        }
        gt_summary.clear(); // Recycle summary object
    }
}
```

```c++
/**<
* Tachyon: https://github.com/mklarqvist/tachyon 
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
*
* This is a complete example!
*/
#include <tachyon/tachyon_reader.h>

std::string my_input_file = "somefile.yon"; // Change me to an actual file that exists on your filesystem
tachyon::TachyonReader reader;
reader.getSettings().loadGenotypes(true);
reader.open(my_input_file);

while(reader.nextBlock()){ // As long as there are YON blocks available
    containers::GenotypeContainer gt(reader.block);
    for(U32 i = 0; i < gt.size(); ++i){
        // All of these functions are in relative terms very expensive!
        // Avoid using them unless you absolutely have to!
        // Vector of literal genotype representations (lower level)
        std::vector<core::GTObject> objects     = gt[i].getLiteralObjects();
        // Vector of genotype objects (high level permuted)
        std::vector<core::GTObject> objects_all = gt[i].getObjects(reader.header.getSampleNumber());
        // Vector of genotype objects (high level unpermuted - original)
        std::vector<core::GTObject> objects_true = gt[i].getObjects(reader.header.getSampleNumber(), reader.block.ppa_manager);

        // Print the difference
        std::cerr << objects.size() << '\t' << objects_all.size() << '\t' << objects_true.size() << std::endl;
        // Dump data
        gt[i].getMeta().toVCFString(std::cout, reader.header, reader.block.index_entry.contigID, reader.block.index_entry.minPosition);
        utility::to_vcf_string(std::cout, objects_true) << '\n';
    }
}
```

---

### Math objects
```c++
/**<
* Tachyon: https://github.com/mklarqvist/tachyon 
* In this example we will calculate the identity-by-descent
* between samples using the math::SquareMatrix object.
*
* This is a complete example!
*/
#include <tachyon/tachyon_reader.h>

std::string my_input_file = "somefile.yon"; // Change me to an actual file that exists on your filesystem
tachyon::TachyonReader reader;
reader.getSettings().loadGenotypes(true);
reader.open(my_input_file);

// Allocate two matrices of size N*N
// The second object is a temporary matrix used when
// genotype data is permuted.
tachyon::math::SquareMatrix<double> square(reader.header.getSampleNumber());
tachyon::math::SquareMatrix<double> square_temporary(reader.header.getSampleNumber());
U64 n_alleles = 0;
while(reader.nextBlock()){ // As long as there are YON blocks available
    containers::GenotypeContainer gt(reader.block);
    // Compare genotypes pairwise for each M in the current YON block
    for(U32 i = 0; i < gt.size(); ++i)
        gt[i].comparePairwise(square_temporary);

    // Add the data from the temporary matrix to the main matrix
    // in the unpermuted genotype order
    square.addUpperTriagonal(square_temporary, reader.block.ppa_manager);
    square_temporary.clear(); // Recycle memory

    // 2 * (Upper triagonal + diagonal) * number of variants
    // This is equivalent to (choose(N, 2) + N) * M_block
    const U64 updates = 2*((reader.header.getSampleNumber()*reader.header.getSampleNumber() - reader.header.getSampleNumber())/2 + reader.header.getSampleNumber()) * gt.size();
    n_aleles += 2*reader.header.getSampleNumber()*gt.size(); // 2*N*M_block
}
square /= n_alleles; // Divide matrix by the number of observed alleles
std::cout << square << std::endl; // Print output
```
Lets plot this output matrix in `R`
```R
# In this example we have pre-computed the identity-by-state
# matrix for the 2,504 samples from the 1000 genomes project
# using genotypes from chromosme 20.
# Calculating this matrix takes roughly ~1 hour on a single
# core on a laptop
#
# Upper triagonal similiary
diff<-read.delim("1kgp3_chr20_ibs_matrix.txt",h=F)
# Square matrix
diff2<-matrix(0,ncol(diff)-1,ncol(diff)-1)

# Utility function to convert upper triagonal matrix into
# square matrix with empty diagonal removed
helper<-function(dataset,position){ t(dataset[position,-position])+dataset[-position,position] }
for(i in 1:(ncol(diff)-1)) diff2[,i]<-helper(diff,i)

# Load sample meta data (including labels)
# Available online at http://www.internationalgenome.org/
groupings<-read.delim("integrated_call_samples_v3.20130502.ALL.panel")

# Generate some colours
library(RColorBrewer)
#colors = rainbow(length(unique(groupings$super_pop)))
colors = brewer.pal(length(unique(groupings$super_pop)), "Accent")
names(colors) = unique(groupings$super_pop)

# Load and run t-SNE with various perplexities
library(Rtsne)
tsneP10 <- Rtsne(diff2, dims = 5, perplexity=10,verbose=TRUE, max_iter = 500)
tsneP20 <- Rtsne(diff2, dims = 5, perplexity=20,verbose=TRUE, max_iter = 500)
tsneP30 <- Rtsne(diff2, dims = 5, perplexity=30,verbose=TRUE, max_iter = 500)
tsneP40 <- Rtsne(diff2, dims = 5, perplexity=40,verbose=TRUE, max_iter = 500)
tsneP50 <- Rtsne(diff2, dims = 5, perplexity=50,verbose=TRUE, max_iter = 500)

# Plot some data
plot(tsneP10$Y[,1],tsneP10$Y[,2],pch=20,cex=.8,col=colors[groupings$super_pop])
legend("topright",legend = names(colors),fill=colors,cex=.6)
plot(tsneP20$Y[,1],tsneP20$Y[,2],pch=20,cex=.8,col=colors[groupings$super_pop])
legend("topright",legend = names(colors),fill=colors,cex=.6)
plot(tsneP30$Y[,1],tsneP30$Y[,2],pch=20,cex=.8,col=colors[groupings$super_pop])
legend("topright",legend = names(colors),fill=colors,cex=.6)
plot(tsneP40$Y[,1],tsneP40$Y[,2],pch=20,cex=.8,col=colors[groupings$super_pop])
legend("topright",legend = names(colors),fill=colors,cex=.6)
plot(tsneP50$Y[,1],tsneP50$Y[,2],pch=20,cex=.8,col=colors[groupings$super_pop])
legend("topright",legend = names(colors),fill=colors,cex=.6)
```
Generated output  
![screenshot](examples/1kgp3_chr20_ibs.png)

[openssl]:  https://www.openssl.org/
[zstd]:     https://github.com/facebook/zstd
[tomahawk]: https://github.com/mklarqvist/tomahawk

### Acknowledgements
[James Bonfield](https://github.com/jkbonfield), Wellcome Trust Sanger Institute  
[Petr Daněček](https://github.com/pd3), Wellcome Trust Sanger Institute  
[Richard Durbin](https://github.com/richarddurbin), Wellcome Trust Sanger Institute, and Department of Genetics, University of Cambridge  

### License
[MIT](LICENSE)
