################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/algorithm/digest/variant_digest_manager.cpp 

OBJS += \
./tachyon/algorithm/digest/variant_digest_manager.o 

CPP_DEPS += \
./tachyon/algorithm/digest/variant_digest_manager.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/algorithm/digest/%.o: ../tachyon/algorithm/digest/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

