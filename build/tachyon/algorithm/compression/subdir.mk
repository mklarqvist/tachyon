################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/algorithm/compression/genotype_encoder.cpp 

OBJS += \
./tachyon/algorithm/compression/genotype_encoder.o 

CPP_DEPS += \
./tachyon/algorithm/compression/genotype_encoder.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/algorithm/compression/%.o: ../tachyon/algorithm/compression/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


