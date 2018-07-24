/*
Copyright (C) 2017-2018 Genome Research Ltd.
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

#ifndef STATS_H_
#define STATS_H_

#include <iostream>
#include <getopt.h>

#include "utility.h"
#include "variant_reader.h"

void stats_usage(void){

}

int stats(int argc, char** argv){
	if(argc <= 2){
		programHelp();
		return(1);
	}

	int c;
	if(argc < 2){
		stats_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",    required_argument, 0, 'i' },
		{"output",   optional_argument, 0, 'o' },
		{"keychain", optional_argument, 0, 'k' },
		{"silent",   no_argument,       0, 's' },
		{0,0,0,0}
	};

	std::string input;
	std::string output;
	std::string keychain_file;
	SILENT = 0;

	while ((c = getopt_long(argc, argv, "i:o:k:s?", long_options, &option_index)) != -1){
		switch (c){
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			input = std::string(optarg);
			break;
		case 'o':
			output = std::string(optarg);
			break;
		case 'k':
			keychain_file = std::string(optarg);
			break;
		case 's':
			SILENT = 1;
			break;
		default:
			std::cerr << tachyon::utility::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
			return(1);
		}
	}

	if(input.length() == 0){
		std::cerr << tachyon::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	// Print messages
	if(!SILENT){
		programMessage();
		std::cerr << tachyon::utility::timestamp("LOG") << "Calling view..." << std::endl;
	}

	tachyon::VariantReader reader;

	if(!reader.open(input)){
		std::cerr << "failed to open" << std::endl;
		return 1;
	}

	tachyon::algorithm::Timer timer, timer2;
	timer.Start(); timer2.Start();

	reader.getBlockSettings().contig(true, true);
	reader.getBlockSettings().positions(true, true);
	reader.getBlockSettings().controller(true, true);
	reader.getBlockSettings().set_membership(true, true);
	reader.getBlockSettings().loadGenotypes(true);
	reader.getBlockSettings().ppa(true, true);
	reader.getBlockSettings().alleles(true, true);

	U32 block_counter = 0;
	std::vector<tachyon::core::TsTvObject> global_titv(reader.getGlobalHeader().getSampleNumber());
	while(reader.nextBlock()){
		reader.getTiTVRatios(global_titv);
		//reader.getGenotypeSummary(std::cout);
		std::cerr << block_counter++ << "/" << reader.getIndex().size() << " in " << timer.ElapsedString() << " " << timer2.ElapsedString() << " " << timer2.Elapsed().count()/(block_counter+1)*reader.getIndex().size() << std::endl;
		timer.Start();
	}

	std::cout << "Sample\tTransversions\tTransitions\tTiTV\tAA\tAT\tAG\tAC\tTA\tTT\tTG\tTC\tGA\tGT\tGG\tGC\tCA\tCT\tCG\tCC\ttotalVariants\tn_insertions\n";
	for(U32 i = 0; i < global_titv.size(); ++i){
		std::cout << reader.getGlobalHeader().samples[i].name << '\t' << global_titv[i] << '\n';
	}

	return 0;
}

#endif /* STATS_H_ */
