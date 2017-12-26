################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/core/BlockEntry.cpp \
../src/core/ImportWriter.cpp \
../src/core/Importer.cpp \
../src/core/StreamContainer.cpp 

OBJS += \
./src/core/BlockEntry.o \
./src/core/ImportWriter.o \
./src/core/Importer.o \
./src/core/StreamContainer.o 

CPP_DEPS += \
./src/core/BlockEntry.d \
./src/core/ImportWriter.d \
./src/core/Importer.d \
./src/core/StreamContainer.d 


# Each subdirectory must supply rules for building sources it contributes
src/core/%.o: ../src/core/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


