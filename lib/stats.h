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

#include "program_utils.h"
#include "variant_reader.h"

void stats_usage(void) {
	programMessage(true);
	std::cerr <<
	"About:  Calculate comprehensive per-sample statistics. Data is written to\n"
	"        stdout as a JSON object.\n"
	"Usage:  " << tachyon::TACHYON_PROGRAM_NAME << " stats [options] -i <in.yon>\n\n"
	"Options:\n"
	"  -i FILE   input YON file (required)\n" << std::endl;
}

int stats(int argc, char** argv) {
	if (argc < 2) {
		programHelp();
		return(1);
	}

	int c;
	if (argc == 2) {
		stats_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",    required_argument, 0, 'i' },
		{"silent",   no_argument,       0, 's' },
		{0,0,0,0}
	};

	tachyon::VariantReaderSettings settings;
	SILENT = 0;

	while ((c = getopt_long(argc, argv, "i:s?", long_options, &option_index)) != -1) {
		switch (c) {
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			settings.input = std::string(optarg);
			break;
		case 'o':
			settings.output = std::string(optarg);
			break;
		case 'k':
			settings.keychain_file = std::string(optarg);
			break;
		case 's':
			SILENT = 1;
			break;
		default:
			std::cerr << tachyon::utility::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
			return(1);
		}
	}

	if (settings.input.length() == 0) {
		std::cerr << tachyon::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	// Print messages
	if (!SILENT) {
		programMessage();
		std::cerr << tachyon::utility::timestamp("LOG") << "Calling stats..." << std::endl;
	}

	tachyon::VariantReader reader;
	reader.GetSettings() = settings;

	if (!reader.open(settings.input)) {
		std::cerr << "failed to open" << std::endl;
		return 1;
	}

	reader.GetBlockSettings().LoadMinimumVcf(true).LoadGenotypes(true);
	reader.Stats();

	return 0;
}

#endif /* STATS_H_ */
