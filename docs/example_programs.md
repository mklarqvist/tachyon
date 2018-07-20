# Example programs
After running `make examples`, the executables versions of the example programs are located in `lib_example`. You can
run these examples using the bundled example dataset located at `examples/example_dataset.yon`. For example, to run `meta_container`
you would have to run a command like
```bash
lib_example/meta_container examples/example_dataset.yon
``` 

Here is a summary of example programs included with Tachyon:
* `calculate_depth_profile <input.yon>`  
   If the input tachyon file has the FORMAT field DP set, then the output will be a matrix of average, standard deviation, minimum and maximum, and total number of non-zero depth for each individual. This example program demonstrates the power of the `SummaryStatistics` objects.
* `format_container_balanced <input.yon>`  
   If the input tachyon file has the FORMAT field GQ set, then the output will print a VCF-string for each variant site. This example uses the balanced FORMAT container that contains
   empty entries to match the number of variants in a block.
* `format_container_raw <input.yon>`  
   If the input tachyon file has the FORMAT field PGT set, then the output will print a VCF-string for each site that has data. This example uses the raw FORMAT container that, unlike its balanced counterpart, do not store empty records. This means that the container has no knowledge of what records go with what sites. This is generally useful when site-information is not directly required. This example also demonstrates how to use the `FormatContainer` with strings.
* `genotype_container <input.yon>`  
   If the input tachyon file has genotypes available, then this example will demonstrate the various internal representations of genotype containers. 1) Return the literal encoded objects: this is useful for most low-level operations but require considerable technical insight; 2) Vector of genotypes for a site but in, potentially, permuted order: this is useful when ordering is not important. Retrieving permuted genotypic vectors are much more efficient that retrieving unpermuted ones; 3) Vector of genotypes for a site in original (unpermuted) order. This example will first print out the number of elements in each of these three containers and then print the content of each.
* `genotype_likelihoods <input.yon>`  
   If the input tachyon file has the FORMAT field PL set, this example will print out the genotype likelihoods for each genotype. This example demonstrates the use of a floating value `FormatContainer`.
* `info_container_balance_comparison <input.yon>`  
   This example program will print out the difference between two `InfoContainer` storing the same data where one is balanced and one is not. The input file has to have the INFO field InbreedingCoeff set.
* `info_container_balanced <input.yon>`  
   If the input tachyon file has the INFO field DP set, then this program will print out the sum depth at each site. Because this is a balanced container the output may contain empty values
* `meta_container <input.yon>`  
   This example demonstrates the use of the `MetaContainer` that stores site-specific and required internal data. This example will print out VCF-strings for the site-specific information for each site.