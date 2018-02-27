################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/containers/checksum_container.cpp \
../tachyon/containers/datacontainer.cpp \
../tachyon/containers/format_container_string.cpp \
../tachyon/containers/genotype_container.cpp \
../tachyon/containers/info_container_string.cpp \
../tachyon/containers/meta_cold_container.cpp \
../tachyon/containers/meta_container.cpp \
../tachyon/containers/meta_hot_container.cpp \
../tachyon/containers/primitive_group_container_string.cpp \
../tachyon/containers/variantblock.cpp 

OBJS += \
./tachyon/containers/checksum_container.o \
./tachyon/containers/datacontainer.o \
./tachyon/containers/format_container_string.o \
./tachyon/containers/genotype_container.o \
./tachyon/containers/info_container_string.o \
./tachyon/containers/meta_cold_container.o \
./tachyon/containers/meta_container.o \
./tachyon/containers/meta_hot_container.o \
./tachyon/containers/primitive_group_container_string.o \
./tachyon/containers/variantblock.o 

CPP_DEPS += \
./tachyon/containers/checksum_container.d \
./tachyon/containers/datacontainer.d \
./tachyon/containers/format_container_string.d \
./tachyon/containers/genotype_container.d \
./tachyon/containers/info_container_string.d \
./tachyon/containers/meta_cold_container.d \
./tachyon/containers/meta_container.d \
./tachyon/containers/meta_hot_container.d \
./tachyon/containers/primitive_group_container_string.d \
./tachyon/containers/variantblock.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/containers/%.o: ../tachyon/containers/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


