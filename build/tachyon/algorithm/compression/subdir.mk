################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/algorithm/compression/compression_manager.cpp \
../tachyon/algorithm/compression/genotype_encoder.cpp \
../tachyon/algorithm/compression/libzpaq.cpp \
../tachyon/algorithm/compression/uncompressed_codec.cpp \
../tachyon/algorithm/compression/zstd_codec.cpp 

OBJS += \
./tachyon/algorithm/compression/compression_manager.o \
./tachyon/algorithm/compression/genotype_encoder.o \
./tachyon/algorithm/compression/libzpaq.o \
./tachyon/algorithm/compression/uncompressed_codec.o \
./tachyon/algorithm/compression/zstd_codec.o 

CPP_DEPS += \
./tachyon/algorithm/compression/compression_manager.d \
./tachyon/algorithm/compression/genotype_encoder.d \
./tachyon/algorithm/compression/libzpaq.d \
./tachyon/algorithm/compression/uncompressed_codec.d \
./tachyon/algorithm/compression/zstd_codec.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/algorithm/compression/%.o: ../tachyon/algorithm/compression/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I"../zstd/lib/" -I"../zstd/lib/common/" -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O2 -msse4.2 -g -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


