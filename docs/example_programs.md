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
* `format_container_raw <input.yon>`  
* `genotype_container <input.yon>`  
* `genotype_likelihoods <input.yon>`  
* `info_container_balance_comparison <input.yon>`  
* `info_container_balanced <input.yon>`  
* `meta_container <input.yon>`  