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
*/
#include <iostream>
#include <getopt.h>

#include "program_utils.h"
#include "variant_importer.h"

void import_usage(void) {
	programMessage();
	std::cerr <<
	"Brief:  Convert Vcf/Bcf records into a Yon archive.\n"
	"Usage:  " << tachyon::TACHYON_PROGRAM_NAME << " import [options] -i <input> -o <output.yon>\n\n"
	"Options:\n"
	"  -i FILE  input Vcf/Bcf/Vcf.gz file (required, \"-\" for piping from stdin)\n"
	"  -o FILE  output file prefix (required, \"-\" for piping to stdout)\n"
	"  -c INT   Import checkpoint size in number of variants (default: 1000)\n"
	"  -C FLOAT Import checkpoint size in bases (default: 5 Mb)\n"
	"  -L INT   Compression level 1-20 (default: 6)\n"
	"  -t INT   Number of compression threads (default: all available)\n"
	"  -p/-P    Permute/Do not permute diploid genotypes\n"
	"  -e       Encrypt data with AES-256\n"
	"  -t       Number of consumer threads for import (default: max available)\n"
	"  -T       Number of (extra) htslib threads for decompression (default: max available)\n"
	"  -s       Hide all program messages [null]\n";
}

int import(int argc, char** argv) {
	if (argc <= 2) {
		import_usage();
		return(1);
	}

	int c;
	if (argc < 2) {
		import_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",               required_argument, 0, 'i' },
		{"output",              optional_argument, 0, 'o' },
		{"checkpoint-variants", optional_argument, 0, 'c' },
		{"checkpoint-bases",    optional_argument, 0, 'C' },
		{"compression-level",   optional_argument, 0, 'L' },
		{"permute",             no_argument,       0, 'p' },
		{"encrypt",             no_argument,       0, 'e' },
		{"no-permute",          no_argument,       0, 'P' },
		{"threads",             optional_argument, 0, 't' },
		{"hts-threads",         optional_argument, 0, 'T' },
		{"silent",              no_argument,       0, 's' },
		{0,0,0,0}
	};

	SILENT = 0;

	tachyon::VariantImporterSettings settings;

	while ((c = getopt_long(argc, argv, "i:o:c:C:L:t:T:sepP?", long_options, &option_index)) != -1) {
		switch (c) {
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			settings.input_file = std::string(optarg);
			break;
		case 'o':
			settings.output_prefix = std::string(optarg);
			break;
		case 'e':
			settings.encrypt_data = true;
			break;
		case 'c':
			settings.checkpoint_n_snps = atoi(optarg);
			if (settings.checkpoint_n_snps <= 0) {
				std::cerr << tachyon::utility::timestamp("ERROR") << "Cannot set checkpoint to <= 0..." << std::endl;
				return(1);
			}
			break;
		case 'C':
			settings.checkpoint_bases = atof(optarg);
			if (settings.checkpoint_bases <= 0) {
				std::cerr << tachyon::utility::timestamp("ERROR") << "Cannot set checkpoint to <= 0..." << std::endl;
				return(1);
			}
			break;
		case 'L':
			settings.compression_level = atoi(optarg);
			if (settings.compression_level <= 0) {
				std::cerr << tachyon::utility::timestamp("ERROR") << "Cannot set compression level to <= 0..." << std::endl;
				return(1);
			}
			break;
		case 'p': settings.permute_genotypes = true;  break;
		case 'P': settings.permute_genotypes = false; break;
		case 't': settings.n_threads = atoi(optarg); break;
		case 'T': settings.htslib_extra_threads = atoi(optarg); break;
		case 's':
			SILENT = 1;
			break;

		default:
			import_usage();
			std::cerr << tachyon::utility::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
			return(1);
		}
	}

	if (settings.input_file.length() == 0) {
		import_usage();
		std::cerr << tachyon::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	// Print messages
	if (!SILENT) {
		programMessage();
		std::cerr << tachyon::utility::timestamp("LOG") << "Calling import..." << std::endl;
	}

	tachyon::VariantImporter importer(settings);
	if (!importer.Build())
		return 1;

	return 0;
}
