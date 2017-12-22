################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/index/IndexBlockEntry.cpp \
../src/index/IndexBlockEntryOffsets.cpp 

OBJS += \
./src/index/IndexBlockEntry.o \
./src/index/IndexBlockEntryOffsets.o 

CPP_DEPS += \
./src/index/IndexBlockEntry.d \
./src/index/IndexBlockEntryOffsets.d 


# Each subdirectory must supply rules for building sources it contributes
src/index/%.o: ../src/index/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/include/openssl/ -I/media/klarqv01/43f895f1-ba84-4abf-b6e7-1510057cc713/Repos2/zstd/lib -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


