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
	std::string my_input_file = "/home/mk21/Repos/tachyon/examples/example_dataset.yon"; // Change me to an actual file that exists on your filesystem
	tachyon::VariantReader reader;

	if(!reader.open(my_input_file)){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to open file: " << my_input_file << "..." << std::endl;
		return(1);
	}

	reader.getBlockSettings().loadFORMAT("GQ");

	/**<
	 *  The `FormatContainer` class stores the data for each variant
	 *  for each individual as container[variant][sample][data]
	 */
	while(reader.nextBlock()){ // As long as there are YON blocks available
		// Meta container
		tachyon::containers::MetaContainer meta(reader.block);

	    // FORMAT container with U32 return type primitive
	    tachyon::containers::FormatContainer<U32>* dp_container = reader.get_balanced_format_container<U32>("GQ", meta);

	    if(dp_container != nullptr){
	        for(U32 variant = 0; variant < dp_container->size(); ++variant){
	            for(U32 sample = 0; sample < dp_container->at(variant).size(); ++sample){
	                // Write the data to `cout` in `VCF` formatting
	                tachyon::utility::to_vcf_string(std::cout, dp_container->at(variant).at(sample)) << ' ';
	            }
	            std::cout << '\n';
	        }
	        std::cout << '\n';
	    }
	    delete dp_container;

	}

	return(0);
}
