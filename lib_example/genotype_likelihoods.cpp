/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/

#include "variant_reader.h"

int main(int argc, char** argv){
	if(argc < 2){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Have to provide an input name..." << std::endl;
		return(1);
	}

	std::string my_input_file(argv[1]);
	tachyon::VariantReader reader;

	reader.getBlockSettings().loadFORMAT("PL");

	if(!reader.open(my_input_file)){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to open file: " << my_input_file << "..." << std::endl;
		return(1);
	}

	/**<
	 *  The `FormatContainer` class stores the data for each variant
	 *  for each individual as container[variant][sample][data]. In this example
	 *  we will use the FORMAT `PL` that generally describes the genotype likelihoods
	 *  as determined by the variant caller
	 */
	while(reader.nextBlock()){ // As long as there are YON blocks available
		// Meta container
		tachyon::containers::MetaContainer meta(reader.getCurrentBlock().getBlock());

	    // FORMAT container with double return type primitive
	    tachyon::containers::FormatContainer<double>* pl_container = reader.getCurrentBlock().get_balanced_format_container<double>("PL", meta);
	    if(pl_container == nullptr) continue;

	    // Iterate over PL container
		for(U32 variant = 0; variant < pl_container->size(); ++variant){
			tachyon::utility::to_vcf_string(std::cout, '\t', meta[variant], reader.getGlobalHeader());
			std::cout << '\t';
			for(U32 sample = 0; sample < pl_container->at(variant).size(); ++sample){
				// Write the data to `cout` in `VCF` formatting
				std::cout << ' ';
				for(U32 entry = 0; entry < pl_container->at(variant).at(sample).size(); ++entry){
					std::cout << "," << pl_container->at(variant).at(sample).at(entry);
				}
				//tachyon::utility::to_vcf_string(std::cout, pl_container->at(variant).at(sample)) << ',';
			}
			std::cout << '\n';
		}
		std::cout << '\n';
	    delete pl_container;
	}

	return(0);
}
