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

	reader.getBlockSettings().loadINFO("InbreedingCoeff");

	if(!reader.open(my_input_file)){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to open file: " << my_input_file << "..." << std::endl;
		return(1);
	}

	/**<
	 * The `InfoContainer` class stores the data for each variant as
	 * container[variant][data]. Both `InfoContainer` and `FormatContainer`
	 * supports variant-balancing of the classes. Balancing refers to filling
	 * variant sites in the file with empty objects if no target data is present
	 * at that site.
	 */
	while(reader.nextBlock()){ // As long as there are YON blocks available
		// Meta container
		tachyon::containers::MetaContainer meta(reader.variant_container.getBlock());

	    // FORMAT container with U32 return type primitive
	    tachyon::containers::InfoContainer<U32>* dp_balanced   = reader.variant_container.get_balanced_info_container<U32>("InbreedingCoeff", meta);
	    tachyon::containers::InfoContainer<U32>* dp_unbalanced = reader.variant_container.get_info_container<U32>("InbreedingCoeff");

	    assert(dp_balanced != nullptr);
	    assert(dp_unbalanced != nullptr);

	    std::cout << "balanced: " << dp_balanced->size() << "\tunbalanced: " << dp_unbalanced->size() << std::endl;

	    delete dp_balanced;
	    delete dp_unbalanced;
	}

	return(0);
}
