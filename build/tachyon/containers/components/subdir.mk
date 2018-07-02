################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/containers/components/variant_block_footer.cpp \
../tachyon/containers/components/variant_block_header.cpp 

OBJS += \
./tachyon/containers/components/variant_block_footer.o \
./tachyon/containers/components/variant_block_header.o 

CPP_DEPS += \
./tachyon/containers/components/variant_block_footer.d \
./tachyon/containers/components/variant_block_header.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/containers/components/%.o: ../tachyon/containers/components/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I"../zstd/lib" -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


