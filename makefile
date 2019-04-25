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
LIBVER_MAJOR_SCRIPT:=`sed -n '/const int32_t TACHYON_VERSION_MAJOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < include/tachyon.h`
LIBVER_MINOR_SCRIPT:=`sed -n '/const int32_t TACHYON_VERSION_MINOR = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < include/tachyon.h`
LIBVER_PATCH_SCRIPT:=`sed -n '/const int32_t TACHYON_VERSION_PATCH = /s/.*[[:blank:]]\([0-9][0-9]*\).*/\1/p' < include/tachyon.h`
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_SCRIPT:= $(LIBVER_MAJOR_SCRIPT).$(LIBVER_MINOR_SCRIPT).$(LIBVER_PATCH_SCRIPT)
LIBVER_MAJOR := $(shell echo $(LIBVER_MAJOR_SCRIPT))
LIBVER_MINOR := $(shell echo $(LIBVER_MINOR_SCRIPT))
LIBVER_PATCH := $(shell echo $(LIBVER_PATCH_SCRIPT))
LIBVER       := $(shell echo $(LIBVER_SCRIPT))

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin
libdir      = $(exec_prefix)/lib
includedir  = $(exec_prefix)/include
#mandir      = $(prefix)/share/man
#man1dir     = $(mandir)/man1

MKDIR_P = mkdir -p
INSTALL = install -p
INSTALL_DATA    = $(INSTALL) -m 644
INSTALL_DIR     = $(MKDIR_P) -m 755
INSTALL_MAN     = $(INSTALL_DATA)
INSTALL_PROGRAM = $(INSTALL)
INSTALL_SCRIPT  = $(INSTALL_PROGRAM)
INSTALL_LIB     = $(INSTALL_DATA)

# If you want to build in debug mode then add DEBUG=true to your build command
# make DEBUG=true
ifdef DEBUG
DEBUG_FLAGS := -g -Wall -Wextra -Wcast-qual -Wcast-align \
                  -Wstrict-aliasing=1 -Wswitch-enum -Wdeclaration-after-statement \
                  -Wstrict-prototypes -Wundef -Wpointer-arith -Wformat-security \
                  -Wvla -Wformat=2 -Winit-self -Wfloat-equal -Wwrite-strings \
                  -Wredundant-decls
else
DEBUG_FLAGS := 
endif

# Global build parameters
INCLUDE_PATH = -I./lib/ -I./include/
ZSTD_LIBRARY_PATH = 

# Check if ZSTD is in the current directory
UNAME_R := $(shell uname -r)
ifneq ("$(wildcard ./zstd/)","")
  INCLUDE_PATH += -I./zstd/lib/ -I./zstd/lib/common/ 
  ZSTD_LIBRARY_PATH = -L./zstd/lib 
else ifneq ("$(wildcard /usr/local/include/zstd.h)","")
  INCLUDE_PATH += -I/usr/local/include/
  #ZSTD_LIBRARY_PATH = -L/usr/local/lib 
else ifneq ("$(wildcard /usr/src/linux-headers-$(UNAME_R)/include/linux/zstd.h)","")
  INCLUDE_PATH += -I/usr/src/linux-headers-$(UNAME_R)/include/linux/
  #ZSTD_LIBRARY_PATH = -L/usr/src/linux-headers-$(uname -r)/lib
else
  INCLUDE_PATH += "-I/usr/src/linux-headers-$(UNAME_R)/include/linux/"
endif


# Try to deduce where OpenSSL is located
OPENSSL_PATH =
ifeq ($(OPENSSL_PATH),)
OPENSSL_LIBRARY_PATH = 
OPENSSL_INCLUDE_PATH =
ifeq ($(OPENSSL_INCLUDE_PATH),)
ifneq ("$(wildcard ./openssl/)","")
  OPENSSL_INCLUDE_PATH = ./openssl/include/
  # INCLUDE_PATH += -I./openssl/include/ 
  OPENSSL_LIBRARY_PATH = -L./openssl/
else ifneq ("$(wildcard /usr/local/include/openssl/)","")
  OPENSSL_INCLUDE_PATH = -I/usr/local/include/
  # INCLUDE_PATH += -I/usr/local/include/
  # OPENSSL_LIBRARY_PATH = -L/usr/local/lib/
else ifneq ("$(wildcard /usr/include/openssl/evp.h)","")
  # INCLUDE_PATH += -I/usr/include/
  OPENSSL_INCLUDE_PATH += -I/usr/include/
  OPENSSL_LIBRARY_PATH = -L/usr/lib/x86_64-linux-gnu/
endif
endif # end if not defined OPENSSL_INCLUDE_PATH
else
  # if OPENSSL_PATH is set
  OPENSSL_INCLUDE_PATH = -I$(OPENSSL_PATH)/include/
  OPENSSL_LIBRARY_PATH = -L$(OPENSSL_PATH)/lib/
endif # end if no OPENSSL_PATH set

INCLUDE_PATH += ${OPENSSL_INCLUDE_PATH}

# Try to deduce where HTSLib is located
HSLIB_LIBRARY_PATH =
ifneq ("$(wildcard ./htslib/)","")
  INCLUDE_PATH += -I./htslib/
  HSLIB_LIBRARY_PATH = -L./htslib/
else ifneq ("$(wildcard /usr/local/include/htslib/)","")
  INCLUDE_PATH += -I/usr/local/include/
  #OPENSSL_LIBRARY_PATH = -L/usr/local/lib/
endif 

# Sort the include_path vector of strings to remove duplicates. This doesn't have
# any functional effect but dedupes the vector.
# Do NOT use equal sign here as the lazy evaluaton will throw a recusion
# warning.
INCLUDE_PATH := $(sort $(INCLUDE_PATH))

# Library paths
LIBRARY_PATHS := $(ZSTD_LIBRARY_PATH) $(OPENSSL_LIBRARY_PATH) $(HSLIB_LIBRARY_PATH) -L/usr/local/lib/

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
LD_LIB_FLAGS = -shared '-Wl,-rpath-link,$$ORIGIN/,-rpath-link,$(PWD),-rpath-link,$$ORIGIN/zstd/lib,-rpath-link,$$ORIGIN/openssl,-rpath-link,$$ORIGIN/htslib,-soname,libtachyon.$(SHARED_EXT)'
else
SHARED_EXT   = dylib
LD_LIB_FLAGS = -dynamiclib -install_name "@rpath/libtachyon.$(SHARED_EXT)" '-Wl,-rpath,@loader_path/,-rpath,$(PWD),-rpath,@loader_path/zstd/lib,-rpath,@loader_path/openssl,-rpath,@loader_path/htslib'
endif

CXXFLAGS      = -std=c++0x $(OPTFLAGS) $(DEBUG_FLAGS)
CFLAGS        = -std=c99   $(OPTFLAGS) $(DEBUG_FLAGS)
CFLAGS_VENDOR = -std=c99   $(OPTFLAGS)

ifneq ($(shell uname), Darwin)
BINARY_RPATHS = '-Wl,-rpath,$$ORIGIN/,-rpath,$(PWD),-rpath,$$ORIGIN/zstd/lib,-rpath,$$ORIGIN/openssl,-rpath,$$ORIGIN/htslib'
else
BINARY_RPATHS = '-Wl,-rpath,@loader_path/,-rpath,$(PWD),-rpath,@loader_path/zstd/lib,-rpath,@loader_path/openssl,-rpath,@loader_path/htslib'
endif

LIBS := -lzstd -lcrypto -lhts
CXX_SOURCE = $(wildcard lib/algorithm/compression/*.cpp) \
			 $(wildcard lib/algorithm/digest/*.cpp) \
			 $(wildcard lib/algorithm/encryption/*.cpp) \
		  	 $(wildcard lib/algorithm/permutation/*.cpp) \
			 $(wildcard lib/containers/*.cpp) \
			 $(wildcard lib/containers/components/*.cpp) \
			 $(wildcard lib/core/header/*.cpp) \
			 $(wildcard lib/core/*.cpp) \
			 $(wildcard lib/index/*.cpp) \
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

OBJECTS  = $(CXX_SOURCE:.cpp=.o) $(C_SOURCE:.c=.o)
CPP_DEPS = $(CXX_SOURCE:.cpp=.d) $(C_SOURCE:.c=.d)

# Inject git information
BRANCH = $(shell git rev-parse --abbrev-ref HEAD)
ifneq ($(BRANCH), master)
GIT_VERSION = $(shell git describe --abbrev=8 --dirty --always --tags)-$(BRANCH)
else
GIT_VERSION = $(shell git describe --abbrev=8 --dirty --always --tags)
endif

# Default target
all: 
	@echo 'Compiling executable...'
	@$(MAKE) tachyon
	@echo 'Compiling library...'
	@$(MAKE) library

# Third party rules
lib/third_party/xxhash/%.o: lib/third_party/xxhash/%.c
	$(CC) $(CFLAGS_VENDOR) -c -o $@ $<

# Generic rules
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -DVERSION=\"$(GIT_VERSION)\" -o $@ $<

tachyon: $(OBJECTS)
	$(CXX) $(BINARY_RPATHS) $(LIBRARY_PATHS) -pthread $(OBJECTS) $(LIBS) -o tachyon
	
library_hack:
	@$(MAKE) cleanmost

library_inner: $(OBJECTS)
	@echo 'Building dynamic library...'
	$(CXX) $(LD_LIB_FLAGS) $(LIBRARY_PATHS) -pthread $(OBJECTS) $(LIBS) -o libtachyon.$(SHARED_EXT).$(LIBVER)
	@echo 'Building static library...'
	$(AR) crs libtachyon.a $(OBJECTS)
	@echo 'Symlinking library...'
	ln -sf libtachyon.$(SHARED_EXT).$(LIBVER) libtachyon.$(SHARED_EXT)
	ln -sf libtachyon.$(SHARED_EXT).$(LIBVER) ltachyon.$(SHARED_EXT)

library:
	@$(MAKE) library_hack
	@$(MAKE) library_inner library=true

examples: $(LIB_EXAMPLE_OUTPUT)

install:
	$(MAKE) tachyon
	$(MAKE) library
	$(INSTALL_DIR) $(DESTDIR)$(bindir) $(DESTDIR)$(includedir)/tachyon $(DESTDIR)$(libdir)
	$(INSTALL_DATA) include/*.h $(DESTDIR)$(includedir)/tachyon
	$(INSTALL_PROGRAM) tachyon $(DESTDIR)$(bindir)
	$(INSTALL_DATA) libtachyon.a $(DESTDIR)$(libdir)/libtachyon.a
	$(INSTALL_LIB) libtachyon.$(SHARED_EXT).$(LIBVER) $(DESTDIR)$(libdir)/libtachyon.$(SHARED_EXT).$(LIBVER)
	ln -sf libtachyon.$(SHARED_EXT).$(LIBVER) $(DESTDIR)$(libdir)/libtachyon.$(SHARED_EXT)
	ln -sf libtachyon.$(SHARED_EXT).$(LIBVER) $(DESTDIR)$(libdir)/ltachyon.$(SHARED_EXT)
  	#$(INSTALL_MAN) doc/tachyon.1 $(DESTDIR)$(man1dir)

# Clean procedures
cleanmost:
	rm -f $(OBJECTS) $(CPP_DEPS)

clean: cleanmost clean_examples
	rm -f tachyon libtachyon.$(SHARED_EXT) libtachyon.$(SHARED_EXT).* ltachyon.$(SHARED_EXT) libtachyon.a

.PHONY: all clean clean_examples cleanmost library library_hack library_inner install