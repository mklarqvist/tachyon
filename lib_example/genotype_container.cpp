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

	if(!reader.open(my_input_file)){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to open file: " << my_input_file << "..." << std::endl;
		return(1);
	}

	reader.getBlockSettings().loadGenotypes(true).loadAllMeta(true).loadPermutationArray(true);

	/**<
	 *  The `GenotypeContaier` class stores the genotype information in a site-centric
	 *  fashion
	 */
	while(reader.nextBlock()){ // As long as there are YON blocks available
		// Meta container
		tachyon::containers::MetaContainer meta(reader.variant_container.getBlock());

	    // Genotype container
		tachyon::containers::GenotypeContainer gt(reader.variant_container.getBlock(), meta);

		for(U32 i = 0; i < gt.size(); ++i){
			// All of these functions are in relative terms very expensive!
			// Avoid using them unless you absolutely have to!
			// Vector of literal genotype representations (lower level)
			std::vector<tachyon::core::GTObject> objects     = gt[i].getLiteralObjects();
			// Vector of genotype objects (high level permuted)
			std::vector<tachyon::core::GTObject> objects_all = gt[i].getObjects(reader.global_header.getSampleNumber());
			// Vector of genotype objects (high level unpermuted - original)
			std::vector<tachyon::core::GTObject> objects_true = gt[i].getObjects(reader.global_header.getSampleNumber(), reader.variant_container.getBlock().ppa_manager);

			// Print the difference
			std::cout << objects.size() << '\t' << objects_all.size() << '\t' << objects_true.size() << std::endl;

			// Dump data
			tachyon::utility::to_vcf_string(std::cout, '\t', gt[i].getMeta(), reader.global_header);
			std::cout << objects[0].n_objects << ":" << objects[0];
			for(U32 sample = 1; sample < objects.size(); ++sample){
				std::cout << '\t' << objects[sample].n_objects << ":" << objects[sample];
			}
			std::cout << '\n';

			tachyon::utility::to_vcf_string(std::cout, '\t', gt[i].getMeta(), reader.global_header);
			std::cout << objects_all[0];
			for(U32 sample = 1; sample < objects_all.size(); ++sample){
				std::cout << '\t' << objects_all[sample];
			}
			std::cout << '\n';

			tachyon::utility::to_vcf_string(std::cout, '\t', gt[i].getMeta(), reader.global_header);
			std::cout << objects_true[0];
			for(U32 sample = 1; sample < objects_true.size(); ++sample){
				std::cout << '\t' << objects_true[sample];
			}
			std::cout << "\n\n";
		}
		std::cout.flush();
	}

	return(0);
}
