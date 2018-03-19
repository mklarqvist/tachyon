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
#include <iostream>
#include <getopt.h>

#include "utility.h"
#include "variant_reader.h"

void view_usage(void){

}

int view(int argc, char** argv){
	if(argc <= 2){
		programMessage();
		programHelpDetailed();
		return(1);
	}

	int c;
	if(argc < 2){
		import_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",		required_argument, 0,  'i' },
		{"output",		optional_argument, 0,  'o' },
		{"noHeader",		no_argument, 0,  'H' },
		{"onlyHeader",		no_argument, 0,  'h' },
		{"dropFormat",	no_argument, 0,  'G' },
		{"dropInfo",	no_argument, 0,  'I' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	std::string input;
	std::string output;
	SILENT = 0;
	bool dropFormat = false;
	bool dropInfo   = false;
	bool headerOnly = false;
	bool showHeader = true;

	while ((c = getopt_long(argc, argv, "i:o:GIshH?", long_options, &option_index)) != -1){
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
		case 'G':
			dropFormat = true;
			break;
		case 'I':
			dropInfo = true;
			break;
		case 's':
			SILENT = 1;
			break;
		case 'h':
			headerOnly = true;
			break;
		case 'H':
			showHeader = false;
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
	//reader.getSettings().loadGenotypes(true);
	//reader.getSettings().loadINFO(true);
	//reader.getSettings().loadAll(true);

	if(!reader.open(input)){
		std::cerr << "failed to open" << std::endl;
		return 1;
	}

	if(headerOnly){
		std::cerr << reader.header.literals.size() << std::endl;
		std::cout << reader.header.literals << std::endl;
		return(0);
	}

	if(showHeader){
		reader.header.writeHeaderString(std::cout, false);
	}

	U64 n_variants = 0;
	tachyon::algorithm::Timer timer;
	timer.Start();

	reader.getSettings().loadAll(true);
	//reader.getSettings().loadMeta(true);

	//tachyon::math::SquareMatrix<double> square(reader.header.n_samples);
	//tachyon::math::SquareMatrix<double> square_temporary(reader.header.n_samples);
	U32 n_blocks = 0;
	//U64 square_division = 0;
	//std::vector<tachyon::core::TsTvObject> global_titv(reader.header.getSampleNumber());
	while(reader.nextBlock()){
		//n_variants += reader.countVariants();
		//n_variants += reader.getTiTVRatios(std::cout, global_titv);
		n_variants += reader.iterate_all_info();
		//n_variants += reader.iterate_genotypes(std::cout);
		//n_variants += reader.iterate_all_info();
		//square_division += reader.calculateIBS(square, square_temporary);
		//std::cerr << n_blocks << '\t' << 597 << std::endl;
		//n_variants += reader.getTiTVRatios(std::cout, global_titv);
		++n_blocks;
		//if(n_blocks == 50) break;
	}

	/*
	std::cout << "Sample\tTransversions\tTransitions\tTiTV\tAA\tAT\tAG\tAC\tTA\tTG\tTC\tTT\tGA\tGT\tGG\tGC\tCA\tCT\tCG\tCC\totalVariantsn";
	for(U32 i = 0; i < global_titv.size(); ++i){
		std::cout << reader.header.samples[i].name << '\t' << global_titv[i] << '\n';
	}
	*/


	//square /= square_division;
	std::cerr << "Blocks: " << n_blocks << std::endl;
	//std::cout << square << std::endl;
	std::cerr << "Variants: " << tachyon::utility::ToPrettyString(n_variants) << " genotypes: " << tachyon::utility::ToPrettyString(n_variants*reader.header.getSampleNumber()) << '\t' << timer.ElapsedString() << '\t' << tachyon::utility::ToPrettyString((U64)((double)n_variants*reader.header.getSampleNumber()/timer.Elapsed().count())) << std::endl;

	return 0;
}
