################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/containers/checksum_container.cpp \
../tachyon/containers/data_container.cpp \
../tachyon/containers/format_container_string.cpp \
../tachyon/containers/genotype_container.cpp \
../tachyon/containers/info_container_string.cpp \
../tachyon/containers/interval_container.cpp \
../tachyon/containers/meta_container.cpp \
../tachyon/containers/primitive_group_container_string.cpp \
../tachyon/containers/variant_block.cpp 

OBJS += \
./tachyon/containers/checksum_container.o \
./tachyon/containers/data_container.o \
./tachyon/containers/format_container_string.o \
./tachyon/containers/genotype_container.o \
./tachyon/containers/info_container_string.o \
./tachyon/containers/interval_container.o \
./tachyon/containers/meta_container.o \
./tachyon/containers/primitive_group_container_string.o \
./tachyon/containers/variant_block.o 

CPP_DEPS += \
./tachyon/containers/checksum_container.d \
./tachyon/containers/data_container.d \
./tachyon/containers/format_container_string.d \
./tachyon/containers/genotype_container.d \
./tachyon/containers/info_container_string.d \
./tachyon/containers/interval_container.d \
./tachyon/containers/meta_container.d \
./tachyon/containers/primitive_group_container_string.d \
./tachyon/containers/variant_block.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/containers/%.o: ../tachyon/containers/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I"../zstd/lib/" -I"../zstd/lib/common/" -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


