################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/core/iterator/ContainerIterator.cpp 

OBJS += \
./src/core/iterator/ContainerIterator.o 

CPP_DEPS += \
./src/core/iterator/ContainerIterator.d 


# Each subdirectory must supply rules for building sources it contributes
src/core/iterator/%.o: ../src/core/iterator/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


