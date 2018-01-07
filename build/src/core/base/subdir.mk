################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/core/base/MetaCold.cpp \
../src/core/base/MetaEntry.cpp 

OBJS += \
./src/core/base/MetaCold.o \
./src/core/base/MetaEntry.o 

CPP_DEPS += \
./src/core/base/MetaCold.d \
./src/core/base/MetaEntry.d 


# Each subdirectory must supply rules for building sources it contributes
src/core/base/%.o: ../src/core/base/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


