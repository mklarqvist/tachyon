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

#include <regex>

#include "tachyon.h"
#include "program_utils.h"
#include "variant_reader.h"

void view_usage(void){
	programMessage(true);
	std::cerr <<
	"About:  Convert YON->VCF/BCF/YON; provides subsetting and slicing functionality\n"
	"Usage:  " << tachyon::TACHYON_PROGRAM_NAME << " view [options] -i <in.yon>\n\n"
	"Options:\n"
	"  -i FILE   input YON file (required)\n"
	"  -o FILE   output file (- for stdout)[-]\n"
	"  -O <y|b|u|z|v> y: tachyon archive, b: compressed BCF, u: uncompressed BCF, \n"
	"                 z: compressed VCF,  v: uncompressed VCF [v]\n"
	"  -c INT    import checkpoint size in number of variants (default: 1000)\n"
	"  -C FLOAT  import checkpoint size in bases (default: 5 Mb)\n"
	"  -L INT    compression level 1-20 (default: 6)\n"
	"  -t INT    number of compression threads (default: all available)\n"
	"  -p/-P     permute/do not permute diploid genotypes\n"
	"  -k FILE   keychain file with encryption keys (required if the file is encrypted)\n"
	"  -f STRING interpreted filter string for slicing output (see manual)\n"
	"  -r STRING interval string (e.g. chr20:10e6-11e6 or chr11:451021)\n"
	"  -b STRING groupings file for Occ-functionality\n"
	//"  -R STRING path to file with interval strings\n"
	"  -G        drop all FORMAT fields from output\n"
	"  -X        annotate FORMAT:GT data and add these statistics to the INFO column\n"
	"            If -b is set then compute these statistics for each grouping."
	"  -h/H      header only / no header\n\n"

	//"Subset options:\n"
	//"  -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" prefix)\n"
	//"  -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n\n"

	"Filter options:\n"
	"  -a/A, --ref-match/--alt-match <REGEX>       regular expression pattern for the reference allele -a or for any alternative alleles -A\n"
	"  -n,   --name-match <REGEX>                  regular expression pattern for the locus name\n"
    "  -s/S, --min-ac/--max-ac <int>               minimum/maximum count for non-reference least frequent\n"
    "                                                 (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]\n"
    //"  -g,   --genotype [^]<hom|het|miss>          require one or more hom/het/missing genotype or, if prefixed with \"^\", exclude sites with hom/het/missing genotypes\n"
    "  -z/Z, --known/--novel                       select known/novel sites only (ID is not/is '.')\n"
	"  -q/Q, --min-quality/--max-quality           minimum/maximum quality value\n"
    "  -m/M, --min-alleles/--max-alleles <int>     minimum/maximum number of alleles listed in REF and ALT\n"
    "  -y/Y, --phased/--exclude-phased             select/exclude sites where all samples are phased\n"
	"  -j/J  --mixed-phasing/--no-mixed-phasing    select sites with both phased and unphased samples\n"
	"  -w/W  --mixed-ploidy/--no-mixed-ploidy      select sites with mixed ploidy\n"
	"  -d/D, --min-af/--max-af <float>             minimum/maximum frequency for non-reference least frequent\n"
    "                                                 (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]\n"
    "  -u/U, --uncalled/--exclude-uncalled         select/exclude sites without a called genotype\n"
	"  -e/E, --remove-unseen/--keep-unseen         select/exclude sites with unseen alternative allele(s)\n\n";
    //"  -v/V, --types/--exclude-types <list>        select/exclude comma-separated list of variant types: snps,indels,mnps,ref,bnd,other [null]\n\n";
}

int view(int argc, char** argv){
	if(argc < 2){
		programHelp();
		return(1);
	}

	int c;
	if(argc == 2){
		view_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",             required_argument, 0,  'i' },
		{"output",            optional_argument, 0,  'o' },
		{"occ",               optional_argument, 0,  'b' },
		{"keychain",          optional_argument, 0,  'k' },
		{"filter",            optional_argument, 0,  'f' },
		{"output-type",       optional_argument, 0,  'O' },
		{"checkpoint-variants", optional_argument, 0, 'c' },
		{"checkpoint-bases",    optional_argument, 0, 'C' },
		{"compression-level",   optional_argument, 0, 'L' },
		{"permute",             no_argument,       0, 'p' },
		{"no-permute",          no_argument,       0, 'P' },
		{"threads",             optional_argument, 0, 't' },

		{"annotate-genotype", no_argument,       0,  'X' },
		{"region",            optional_argument, 0,  'r' },
		{"noHeader",          no_argument,       0,  'H' },
		{"onlyHeader",        no_argument,       0,  'h' },
		{"dropFormat",        no_argument,       0,  'G' },
		{"af-min",            optional_argument, 0,  'd' },
		{"af-max",            optional_argument, 0,  'D' },
		{"ac-min",            optional_argument, 0,  's' },
		{"ac-max",            optional_argument, 0,  'S' },
		{"alleles-min",       optional_argument, 0,  'm' },
		{"alleles-max",       optional_argument, 0,  'M' },
		{"known",             no_argument,       0,  'z' },
		{"novel",             no_argument,       0,  'Z' },
		{"phased",            no_argument,       0,  'y' },
		{"exclude-phased"    ,no_argument,       0,  'Y' },
		{"mixed-phase",       no_argument,       0,  'j' },
		{"no-mixed-phase",    no_argument,       0,  'J' },
		{"uncalled",          no_argument,       0,  'u' },
		{"exclude-uncalled",  no_argument,       0,  'U' },
		{"ref-match",         optional_argument, 0,  'a' },
		{"alt-match",         optional_argument, 0,  'A' },
		{"name-match",        optional_argument, 0,  'n' },
		{"mixed-ploidy",      no_argument,       0,  'w' },
		{"no-mixed-ploidy",   no_argument,       0,  'W' },
		{"remove-unseen",     no_argument,       0,  'e' },
		{"keep-unseen",       no_argument,       0,  'E' },
		{"min-quality",       optional_argument, 0,  'q' },
		{"max-quality",       optional_argument, 0,  'Q' },
		{0,0,0,0}
	};

	tachyon::VariantReaderSettings settings;
	tachyon::yon_vb_settings block_settings;
	std::vector<std::string> interpret_commands;
	std::vector<std::string> interval_strings;

	SILENT = 0;
	std::string temp;
	tachyon::VariantReader reader;
	tachyon::VariantReaderFilters& filters = reader.GetFilterSettings();

	while ((c = getopt_long(argc, argv, "i:o:k:f:O:r:GshHX?l:L:m:M:pPuUc:C:jJzZa:A:n:wWeEq:Q:b:", long_options, &option_index)) != -1){
		switch (c){
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			settings.input = std::string(optarg);
			break;
		case 'o':
			settings.output = std::string(optarg);
			break;
		case 'c':
			settings.checkpoint_n_snps = atoi(optarg);
			if(settings.checkpoint_n_snps <= 0){
				std::cerr << tachyon::utility::timestamp("ERROR") << "Cannot set checkpoint to <= 0..." << std::endl;
				return(1);
			}
			break;
		case 'C':
			settings.checkpoint_bases = atof(optarg);
			if(settings.checkpoint_bases <= 0){
				std::cerr << tachyon::utility::timestamp("ERROR") << "Cannot set checkpoint to <= 0..." << std::endl;
				return(1);
			}
			break;
		case 'L':
			settings.compression_level = atoi(optarg);
			if(settings.compression_level <= 0){
				std::cerr << tachyon::utility::timestamp("ERROR") << "Cannot set compression level to <= 0..." << std::endl;
				return(1);
			}
			break;
		case 'p': settings.permute_genotypes = true;  break;
		case 'P': settings.permute_genotypes = false; break;
		case 'b':
			settings.group_file = std::string(optarg);
			break;
		case 'k':
			settings.keychain_file = std::string(optarg);
			break;
		case 'd':
			filters.Add(tachyon::YON_FILTER_ALLELE_FREQUENCY, atof(optarg), tachyon::YON_CMP_GREATER);
			break;
		case 'D':
			filters.Add(tachyon::YON_FILTER_ALLELE_FREQUENCY, atof(optarg), tachyon::YON_CMP_LESS_EQUAL);
			break;
		case 'q':
			filters.Add(tachyon::YON_FILTER_QUALITY, atof(optarg), tachyon::YON_CMP_GREATER);
			break;
		case 'Q':
			filters.Add(tachyon::YON_FILTER_QUALITY, atof(optarg), tachyon::YON_CMP_LESS_EQUAL);
			break;
		case 'm':
			filters.Add(tachyon::YON_FILTER_NUMBER_ALT_ALLELES, atoi(optarg), tachyon::YON_CMP_GREATER);
			break;
		case 'M':
			filters.Add(tachyon::YON_FILTER_NUMBER_ALT_ALLELES, atoi(optarg), tachyon::YON_CMP_LESS_EQUAL);
			break;
		case 'a':
			filters.Add(tachyon::YON_FILTER_REFERENCE_ALLELE, std::string(optarg), tachyon::YON_CMP_REGEX);
			break;
		case 'A':
			filters.Add(tachyon::YON_FILTER_ALT_ALLELE, std::string(optarg), tachyon::YON_CMP_REGEX);
			break;
		case 'n':
			filters.Add(tachyon::YON_FILTER_NAME, std::string(optarg), tachyon::YON_CMP_REGEX);
			break;
		case 's':
			filters.Add(tachyon::YON_FILTER_ALLELE_COUNT, atoi(optarg), tachyon::YON_CMP_GREATER);;
			break;
		case 'S':
			filters.Add(tachyon::YON_FILTER_ALLELE_COUNT, atoi(optarg), tachyon::YON_CMP_LESS_EQUAL);
			break;
		case 'y':
			filters.Add(tachyon::YON_FILTER_UNIFORM_PHASE, (bool)true, tachyon::YON_CMP_EQUAL);
			break;
		case 'Y':
			filters.Add(tachyon::YON_FILTER_UNIFORM_PHASE, (bool)false, tachyon::YON_CMP_EQUAL);
			break;
		case 'j':
			filters.Add(tachyon::YON_FILTER_MIXED_PHASING, (bool)true, tachyon::YON_CMP_EQUAL);
			break;
		case 'J':
			filters.Add(tachyon::YON_FILTER_MIXED_PHASING, (bool)false, tachyon::YON_CMP_EQUAL);
			break;
		case 'w':
			filters.Add(tachyon::YON_FILTER_MIXED_PLOIDY, (bool)true, tachyon::YON_CMP_EQUAL);
			break;
		case 'W':
			filters.Add(tachyon::YON_FILTER_MIXED_PLOIDY, (bool)false, tachyon::YON_CMP_EQUAL);
			break;
		case 'u':
			filters.Add(tachyon::YON_FILTER_MISSING_GT, 0, tachyon::YON_CMP_GREATER);
			break;
		case 'U':
			filters.Add(tachyon::YON_FILTER_MISSING_GT, 0, tachyon::YON_CMP_EQUAL);
			break;
		case 'z':
			filters.Add(tachyon::YON_FILTER_KNOWN_NOVEL, (bool)true, tachyon::YON_CMP_EQUAL);
			break;
		case 'Z':
			filters.Add(tachyon::YON_FILTER_KNOWN_NOVEL, (bool)false, tachyon::YON_CMP_EQUAL);
			break;
		case 'f':
			interpret_commands.push_back(std::string(optarg));
			break;
		case 'G':
			settings.drop_format = true;
			break;
		case 'r':
			interval_strings.push_back(std::string(optarg));
			break;
		case 'h':
			settings.header_only = true;
			break;
		case 'H':
			settings.show_header = false;
			break;
		case 'O':
			settings.output_type = optarg[0];
			break;
		case 'X':
			settings.annotate_genotypes = true;
			break;
		case 'e':
			filters.Add(tachyon::YON_FILTER_UNSEEN_ALT, (bool)true, tachyon::YON_CMP_EQUAL);
			break;
		case 'E':
			filters.Add(tachyon::YON_FILTER_UNSEEN_ALT, (bool)false, tachyon::YON_CMP_EQUAL);
			break;

		default:
			std::cerr << tachyon::utility::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
			return(1);
		}
	}

	if(settings.input.length() == 0){
		std::cerr << tachyon::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	if(!reader.open(settings.input)){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to open file: " << settings.input << "..." << std::endl;
		return 1;
	}

	if(settings.header_only){
		reader.UpdateHeaderView();
		reader.GetHeader().PrintVcfHeader(std::cout);
		return(0);
	}

	// User provided '-f' string(s)
	if(interpret_commands.size()){
		if(!reader.GetBlockSettings().ParseCommandString(interpret_commands, reader.GetHeader())){
			std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to parse command..." << std::endl;
			return(1);
		}
	} else {
		reader.GetBlockSettings().LoadAll(true);

		if(settings.drop_format){
			reader.GetBlockSettings().LoadGenotypes(false).LoadDisplayWrapper(false, YON_BLK_BV_PPA).LoadDisplayWrapper(false, YON_BLK_BV_FORMAT);
		}
	}

	if(settings.output_type == 'v'){
		settings.use_htslib = false;
	} else if(settings.output_type == 'z'){
		settings.output_type = 'z';
		settings.use_htslib = true;
	} else if(settings.output_type == 'b'){
		settings.output_type = 'b';
		settings.use_htslib = true;
	} else if(settings.output_type == 'u'){
		settings.output_type = 'u';
		settings.use_htslib = true;
	} else if(settings.output_type == 'y'){
		settings.output_type = 'y';
		settings.use_htslib = false;
	} else {
		std::cerr << tachyon::utility::timestamp("ERROR") << "Unrecognised output option: " << settings.output_type << "..." << std::endl;
		return(1);
	}

	// If user is triggering annotation
	if(settings.annotate_genotypes){
		reader.GetBlockSettings().annotate_extra = true;
		reader.GetBlockSettings().LoadGenotypes(true).LoadMinimumVcf(true);
		if(settings.drop_format) reader.GetBlockSettings().DisplayWrapper(false, YON_BLK_BV_GT);
	}

	if(filters.HasRequireGenotypes()){
		reader.GetBlockSettings().LoadGenotypes(true).LoadMinimumVcf(true);
		if(settings.drop_format) reader.GetBlockSettings().DisplayWrapper(false, YON_BLK_BV_GT);
	}

	reader.GetSettings() = settings;

	if(settings.show_header) reader.GetBlockSettings().show_vcf_header = true;
	else reader.GetBlockSettings().show_vcf_header = false;

	if(reader.AddIntervals(interval_strings) == false) return(1);


	reader.OutputRecords();

	return 0;
}
