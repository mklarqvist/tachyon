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
LIBVER_MAJOR_SCRIPT:=`sed -n '/const S32 TACHYON_VERSION_MAJOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < lib/support/MagicConstants.h`
LIBVER_MINOR_SCRIPT:=`sed -n '/const S32 TACHYON_VERSION_MINOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < lib/support/MagicConstants.h`
LIBVER_PATCH_SCRIPT:=`sed -n '/const S32 TACHYON_VERSION_RELEASE = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < lib/support/MagicConstants.h`
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
INCLUDE_PATH := -I"lib/" -I"../zstd/lib/" -I"../zstd/lib/common/" -I"/usr/local/opt/openssl/lib/" -I"/usr/include/openssl/" -I"/usr/local/include/"
OPTFLAGS := -O3 -msse4.2
# Legacy flags used
#OPTFLAGS := -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops

ifdef library
OPTFLAGS += -fPIC
endif

# OS X linker doesn't support -soname, and use different extension
# see : https://developer.apple.com/library/mac/documentation/DeveloperTools/Conceptual/DynamicLibraries/100-Articles/DynamicLibraryDesignGuidelines.html
ifneq ($(shell uname), Darwin)
SHARED_EXT = so
LD_LIB_FLAGS := -shared -Wl,-rpath,"./",-soname,libtachyon.$(SHARED_EXT)
else
SHARED_EXT = dylib
LD_LIB_FLAGS := -dynamiclib -install_name libtachyon.$(SHARED_EXT)
endif

CXXFLAGS := -std=c++0x $(OPTFLAGS) $(DEBUG_FLAGS)
CFLAGS   := -std=c99 $(OPTFLAGS) $(DEBUG_FLAGS)

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
lib/third_party/zlib/compress.c \
lib/third_party/zlib/crc32.c \
lib/third_party/zlib/deflate.c \
lib/third_party/zlib/gzclose.c \
lib/third_party/zlib/gzlib.c \
lib/third_party/zlib/gzread.c \
lib/third_party/zlib/gzwrite.c \
lib/third_party/zlib/infback.c \
lib/third_party/zlib/inffast.c \
lib/third_party/zlib/inflate.c \
lib/third_party/zlib/inftrees.c \
lib/third_party/zlib/trees.c \
lib/third_party/zlib/uncompr.c \
lib/third_party/zlib/zutil.c 

OBJECTS = $(CXX_SOURCE:.cpp=.o) $(C_SOURCE:.c=.o)
CPP_DEPS = $(CXX_SOURCE:.cpp=.d)

# Inject git information
BRANCH := $(shell git rev-parse --abbrev-ref HEAD)
ifneq ($(BRANCH), master)
GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags)-$(BRANCH)
else
GIT_VERSION := $(shell git describe --abbrev=8 --dirty --always --tags)
endif

# All Target
all: tachyon

# Third party rules
lib/third_party/xxhash/%.o: lib/third_party/xxhash/%.c
	gcc $(CFLAGS) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"

lib/third_party/zlib/%.o: lib/third_party/zlib/%.c
	gcc $(CFLAGS) -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"

# Generic rules
%.o: %.cpp
	g++ $(CXXFLAGS) $(INCLUDE_PATH) -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"

tachyon: $(OBJECTS)
	g++ -Wl,-rpath,"$(PWD)/zstd/lib/" -L"zstd/lib" -pthread -o "tachyon" $(OBJECTS) $(LIBS)
	$(MAKE) cleanmost
	$(MAKE) library library=true

library: $(OBJECTS)
	@echo 'Building with positional independence...'
	g++ $(LD_LIB_FLAGS) -pthread -o libtachyon.$(SHARED_EXT).$(LIBVER) $(OBJECTS) $(LIBS)
	@echo 'Symlinking library...'
	ln -sf libtachyon.$(SHARED_EXT).$(LIBVER) libtachyon.$(SHARED_EXT)
	ln -sf libtachyon.$(SHARED_EXT) ltachyon.$(SHARED_EXT)

cleanmost:
	rm -f $(OBJECTS) $(CPP_DEPS)

clean: cleanmost
	rm -f tachyon libtachyon.so libtachyon.so.* ltachyon.so

.PHONY: all clean cleanmost library