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

#include "import.h"
#include "stats.h"
#include "program_utils.h"
#include "view.h"

int main(int argc, char** argv) {
	if (tachyon::utility::IsBigEndian()) {
		std::cerr << tachyon::utility::timestamp("ERROR") << "Tachyon does not support big endian systems..." << std::endl;
		return(1);
	}

	if (argc == 1) {
		programHelp();
		return(1);
	}

	// Literal string input line
	tachyon::LITERAL_COMMAND_LINE = tachyon::TACHYON_PROGRAM_NAME;
	for (uint32_t i = 1; i < argc; ++i)
		tachyon::LITERAL_COMMAND_LINE += " " + std::string(&argv[i][0]);

	const std::string subroutine(argv[1]);

	if (strncmp(subroutine.data(), "import", 6) == 0 && subroutine.size() == 6) {
		return(import(argc, argv));
	} else if (strncmp(subroutine.data(), "view", 4) == 0 && subroutine.size() == 4) {
		return(view(argc, argv));
	} else if (strncmp(subroutine.data(), "stats", 5) == 0 && subroutine.size() == 5) {
		return(stats(argc, argv));
		return(0);
	} else if (strncmp(subroutine.data(), "check", 5) == 0 && subroutine.size() == 5) {
		std::cerr << tachyon::utility::timestamp("ERROR") << "Not implemented" << std::endl;
		return(0);
	} else {
		programHelp();
		std::cerr << tachyon::utility::timestamp("ERROR") << "Illegal command: " << subroutine << std::endl;
		return(1);
	}
	return(1);
}
