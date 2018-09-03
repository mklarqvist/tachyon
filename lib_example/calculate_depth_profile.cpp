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

	reader.GetBlockSettings().LoadFormat("DP");

	if(!reader.open(my_input_file)){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to open file: " << my_input_file << "..." << std::endl;
		return(1);
	}

	std::vector<tachyon::math::SummaryStatistics> depth_data(reader.GetGlobalHeader().GetNumberSamples());

	/**<
	 *  In this example we will write a simple program to calculate
	 *  the depth profile for each individual in a file
	 */
	while(reader.NextBlock()){ // As long as there are YON blocks available
		// Meta container
		tachyon::containers::MetaContainer meta(reader.GetCurrentContainer());

	    // FORMAT container with U32 return type primitive
	    tachyon::containers::FormatContainer<uint32_t>* dp_container = reader.GetCurrentContainer().get_balanced_format_container<uint32_t>("DP", meta);

	    if(dp_container != nullptr){
	        for(int variant = 0; variant < dp_container->size(); ++variant){
	            for(int sample = 0; sample < dp_container->at(variant).size(); ++sample){
	               if(dp_container->at(variant).at(sample).size()){
	            		depth_data[sample].addNonzero(dp_container->at(variant).at(sample)[0]);
	            	}
	            }
	        }
	    }
	    delete dp_container;
	}

	std::cout << "Sample\tMean\tSD\tMin\tMax\tN\n";
	for(int i = 0; i < reader.GetGlobalHeader().GetNumberSamples(); ++i){
		depth_data[i].calculate();
		std::cout << reader.GetGlobalHeader().samples_[i] << "\t" << depth_data[i].mean << "\t" << depth_data[i].getStandardDeviation() << "\t" << depth_data[i].min << "\t" << depth_data[i].max << "\t" << depth_data[i].getCount() << "\n";
	}
	std::cout.flush();

	return(0);
}
