################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/index/BlockIndexEntry.cpp 

OBJS += \
./src/index/BlockIndexEntry.o 

CPP_DEPS += \
./src/index/BlockIndexEntry.d 


# Each subdirectory must supply rules for building sources it contributes
src/index/%.o: ../src/index/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/Cellar/openssl/1.0.2n/include/ -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


