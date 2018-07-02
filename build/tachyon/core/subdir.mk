################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/core/meta_allele.cpp \
../tachyon/core/meta_entry.cpp \
../tachyon/core/variant_import_writer.cpp 

OBJS += \
./tachyon/core/meta_allele.o \
./tachyon/core/meta_entry.o \
./tachyon/core/variant_import_writer.o 

CPP_DEPS += \
./tachyon/core/meta_allele.d \
./tachyon/core/meta_entry.d \
./tachyon/core/variant_import_writer.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/core/%.o: ../tachyon/core/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I"../zstd/lib/" -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


