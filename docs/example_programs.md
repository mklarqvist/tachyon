# Example programs
After running `make examples`, the executables versions of the example programs are located in `lib_example`. You can
run these examples using the bundled example dataset located at `examples/example_dataset.yon`. For example, to run `meta_container`
you would have to run a command like
```bash
lib_example/meta_container examples/example_dataset.yon
``` 

Here is a summary of example programs included with Tachyon:
* `calculate_depth_profile <input.yon>`  
   If input.vcf is a VCF file in which the variant calls have allele depths (in their 'AD' FORMAT fields), then output.vcf will be the same but with the allele depths summed into the 'AD' INFO field of the variants. This program is a good example of writing out a file with a modified header, as well as demonstrating the variant_utils and variantcall_utils routines for setting and getting INFO and FORMAT information.
* `format_container_balanced <input.yon>`  
* `genotype_container <input.yon>`  
* `genotype_likelihoods <input.yon>`  
* `info_container_balance_comparison <input.yon>`  
* `info_container_balanced <input.yon>`  
* `meta_container <input.yon>`  