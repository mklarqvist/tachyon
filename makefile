# ################################################################
# Copyright (C) 2017-present Genome Research Ltd.
# Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
# ################################################################

# Version numbers slices from the source header
LIBVER_MAJOR_SCRIPT:=`sed -n '/const S32 TACHYON_VERSION_MAJOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < lib/support/magic_constants.h`
LIBVER_MINOR_SCRIPT:=`sed -n '/const S32 TACHYON_VERSION_MINOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < lib/support/magic_constants.h`
LIBVER_PATCH_SCRIPT:=`sed -n '/const S32 TACHYON_VERSION_PATCH = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < lib/support/magic_constants.h`
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_MAJOR := $(shell echo $(LIBVER_MAJOR_SCRIPT))
LIBVER_MINOR := $(shell echo $(LIBVER_MINOR_SCRIPT))
LIBVER_PATCH := $(shell echo $(LIBVER_PATCH_SCRIPT))
LIBVER       := $(shell echo $(LIBVER_SCRIPT))

PREFIX := /usr/local

# If you want to build in debug mode then add DEBUG=true to your build command
# make DEBUG=true
ifdef DEBUG
DEBUG_FLAGS := -g -Wall -Wextra -Wcast-qual -Wcast-align -Wshadow \
                  -Wstrict-aliasing=1 -Wswitch-enum -Wdeclaration-after-statement \
                  -Wstrict-prototypes -Wundef -Wpointer-arith -Wformat-security \
                  -Wvla -Wformat=2 -Winit-self -Wfloat-equal -Wwrite-strings \
                  -Wredundant-decls
else
DEBUG_FLAGS := 
endif

# Global build parameters
INCLUDE_PATH = -I./lib/
ZSTD_LIBRARY_PATH = 

# Check if ZSTD is in the current directory
ifneq ("$(wildcard ./zstd/)","")
  INCLUDE_PATH += -I./zstd/lib/ -I./zstd/lib/common/ 
  ZSTD_LIBRARY_PATH = -L./zstd/lib 
else ifneq ("$(wildcard /usr/local/include/)","")
  INCLUDE_PATH += -I/usr/local/include/
  #ZSTD_LIBRARY_PATH = -L/usr/local/lib 
endif

# Try to deduce where OpenSSL is located
OPENSSL_LIBRARY_PATH = 
ifneq ("$(wildcard ./openssl/)","")
  INCLUDE_PATH += -I./openssl/include/ 
  OPENSSL_LIBRARY_PATH = -L./openssl/
else ifneq ("$(wildcard /usr/local/include/openssl/)","")
  INCLUDE_PATH += -I/usr/local/include/
  #OPENSSL_LIBRARY_PATH = -L/usr/local/lib/
else ifneq ("$(wildcard /usr/include/openssl/evp.h)","")
  INCLUDE_PATH += -I/usr/include/
  OPENSSL_LIBRARY_PATH = -L/usr/lib/x86_64-linux-gnu/
endif

LIBRARY_PATHS := $(ZSTD_LIBRARY_PATH) $(OPENSSL_LIBRARY_PATH) -L/usr/local/lib/

OPTFLAGS := -O3 -msse4.2
# Legacy flags used
#OPTFLAGS := -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops

ifdef library
OPTFLAGS += -fPIC
endif

# OS X linker doesn't support -soname, and use different extension
# see : https://developer.apple.com/library/mac/documentation/DeveloperTools/Conceptual/DynamicLibraries/100-Articles/DynamicLibraryDesignGuidelines.html
ifneq ($(shell uname), Darwin)
SHARED_EXT   = so
LD_LIB_FLAGS = -shared -Wl,-rpath,./zstd/lib,-rpath,./openssl/,-soname,libtachyon.$(SHARED_EXT)
else
SHARED_EXT   = dylib
LD_LIB_FLAGS = -dynamiclib -install_name libtachyon.$(SHARED_EXT) -Wl,-rpath,./zstd/lib,-rpath,./openssl/
endif

CXXFLAGS      = -std=c++0x $(OPTFLAGS) $(DEBUG_FLAGS)
CFLAGS        = -std=c99   $(OPTFLAGS) $(DEBUG_FLAGS)
BINARY_RPATHS = '-Wl,-rpath,$$ORIGIN/zstd/lib,-rpath,$$ORIGIN/openssl/'

LIBS := -lzstd -lcrypto
CXX_SOURCE = $(wildcard lib/algorithm/compression/*.cpp) \
			 $(wildcard lib/algorithm/digest/*.cpp) \
			 $(wildcard lib/algorithm/encryption/*.cpp) \
		  	 $(wildcard lib/algorithm/permutation/*.cpp) \
			 $(wildcard lib/containers/*.cpp) \
			 $(wildcard lib/containers/components/*.cpp) \
			 $(wildcard lib/core/header/*.cpp) \
			 $(wildcard lib/core/*.cpp) \
			 $(wildcard lib/io/*.cpp) \
			 $(wildcard lib/io/bcf/*.cpp) \
			 $(wildcard lib/io/compression/*.cpp) \
			 $(wildcard lib/io/vcf/*.cpp) \
			 $(wildcard lib/*.cpp) \
			 $(wildcard lib/math/*.cpp) \
			 $(wildcard lib/support/*.cpp) \
			 $(wildcard lib/utility/*.cpp)

C_SOURCE = \
lib/third_party/xxhash/xxhash.c \
lib/third_party/zlib/adler32.c \
lib/third_party/zlib/crc32.c \
lib/third_party/zlib/deflate.c \
lib/third_party/zlib/infback.c \
lib/third_party/zlib/inffast.c \
lib/third_party/zlib/inflate.c \
lib/third_party/zlib/inftrees.c \
lib/third_party/zlib/trees.c \
lib/third_party/zlib/zutil.c \
lib/third_party/zlib/compress.c \
lib/third_party/zlib/uncompr.c \
lib/third_party/zlib/gzclose.c \
lib/third_party/zlib/gzlib.c \
lib/third_party/zlib/gzread.c \
lib/third_party/zlib/gzwrite.c \

OBJECTS  = $(CXX_SOURCE:.cpp=.o) $(C_SOURCE:.c=.o)
CPP_DEPS = $(CXX_SOURCE:.cpp=.d) $(C_SOURCE:.c=.d)

LIB_INCLUDE_PATH   = -I./lib/
LIB_EXAMPLE_FLAGS  = -L./ -ltachyon '-Wl,-rpath,$$ORIGIN/../,-rpath,$(PWD),-rpath,$$ORIGIN/../zstd/lib,-rpath,$$ORIGIN/../openssl'
LIB_EXAMPLE_SOURCE = $(wildcard lib_example/*.cpp)
LIB_EXAMPLE_OUTPUT = $(LIB_EXAMPLE_SOURCE:.cpp=)

# Inject git information
BRANCH = $(shell git rev-parse --abbrev-ref HEAD)
ifneq ($(BRANCH), master)
GIT_VERSION = $(shell git describe --abbrev=8 --dirty --always --tags)-$(BRANCH)
else
GIT_VERSION = $(shell git describe --abbrev=8 --dirty --always --tags)
endif

# Default target
all: tachyon

# Third party rules
lib/third_party/xxhash/%.o: lib/third_party/xxhash/%.c
	gcc $(CFLAGS) -c -o $@ $<

lib/third_party/zlib/%.o: lib/third_party/zlib/%.c
	gcc $(CFLAGS) -c -o $@ $<

# Generic rules
%.o: %.cpp
	g++ $(CXXFLAGS) $(INCLUDE_PATH) -c -DVERSION=\"$(GIT_VERSION)\" -o $@ $<

tachyon: $(OBJECTS)
	g++ $(BINARY_RPATHS) $(LIBRARY_PATHS) -pthread $(OBJECTS) $(LIBS) -o tachyon
	$(MAKE) cleanmost
	$(MAKE) library library=true
	$(MAKE) examples

library: $(OBJECTS)
	@echo 'Building dynamic library...'
	g++ $(LD_LIB_FLAGS) $(LIBRARY_PATHS) -pthread $(OBJECTS) $(LIBS) -o libtachyon.$(SHARED_EXT).$(LIBVER)
	@echo 'Building static library...'
	ar crs libtachyon.a $(OBJECTS)
	@echo 'Symlinking library...'
	ln -sf libtachyon.$(SHARED_EXT).$(LIBVER) libtachyon.$(SHARED_EXT)
	ln -sf libtachyon.$(SHARED_EXT).$(LIBVER) ltachyon.$(SHARED_EXT)

examples: $(LIB_EXAMPLE_OUTPUT)

lib_example/%: lib_example/%.cpp
	g++ $(CXXFLAGS) $(INCLUDE_PATH) $(LIB_INCLUDE_PATH) $(LIB_EXAMPLE_FLAGS) -DVERSION=\"$(GIT_VERSION)\" -o $@ $<

# Clean procedures
clean_examples:
	rm -f $(LIB_EXAMPLE_OUTPUT)

cleanmost:
	rm -f $(OBJECTS) $(CPP_DEPS)

clean: cleanmost clean_examples
	rm -f tachyon libtachyon.so libtachyon.so.* ltachyon.so

.PHONY: all clean clean_examples cleanmost library examples