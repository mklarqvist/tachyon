#!/bin/bash
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

function note_build_stage {
  echo "========== [$(date)] Stage '${1}' starting"
}

if [[ "$1" == "local" ]]; then
note_build_stage "Building zstd..."
if [ ! -d zstd ]; then
git clone https://github.com/facebook/zstd
fi
cd zstd
if [ ! -f lib/libzstd.so ]; then
    make -j$(nproc)
else
    echo "ZSTD already built! Skipping..."
fi
cd ..

 note_build_stage "Building OpenSSL..."
if [ ! -d openssl ]; then
git clone https://github.com/openssl/openssl.git
fi
cd openssl
if [ ! -f libssl.so ]; then
    ./config
    make -j$(nproc)
else
    echo "OpenSSL already built! Skipping..."
fi
cd ..

 note_build_stage "Building htslib"
if [ ! -d htslib ]; then
git clone https://github.com/samtools/htslib.git
fi
cd htslib
if [ ! -f htslib.so ]; then
    autoheader && autoconf && ./configure && make -j$(nproc)
else
    echo "htslib already built! Skipping..."
fi
cd ..

else # Install with sudo
if [ "$(uname)" == "Darwin" ]; then
    # Update package list
    ################################################################################
    brew update
    # Install generic dependencies
    ################################################################################
    note_build_stage "Update misc. dependencies"
    brew install openssl
    brew install zstd
else
    # Update package list
    ################################################################################
    note_build_stage "Update package list"
    sudo -H apt-get -qq -y update

    # Install generic dependencies
    ################################################################################
    note_build_stage "Update misc. dependencies"
    sudo -H apt-get -y install pkg-config zip g++ zlib1g-dev unzip curl git lsb-release liblz4-dev zstd

    # Install htslib dependencies
    ################################################################################
    note_build_stage "Install htslib dependencies"
    sudo -H apt-get -y install libssl-dev libcurl4-openssl-dev liblz-dev libbz2-dev liblzma-dev
fi
fi

# Install Tachyon
################################################################################
note_build_stage "Building Tachyon"
make clean; make -j$(nproc)