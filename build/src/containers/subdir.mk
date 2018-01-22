################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/containers/Block.cpp \
../src/containers/Container.cpp \
../src/containers/MetaColdContainer.cpp \
../src/containers/MetaContainer.cpp \
../src/containers/MetaHotContainer.cpp 

OBJS += \
./src/containers/Block.o \
./src/containers/Container.o \
./src/containers/MetaColdContainer.o \
./src/containers/MetaContainer.o \
./src/containers/MetaHotContainer.o 

CPP_DEPS += \
./src/containers/Block.d \
./src/containers/Container.d \
./src/containers/MetaColdContainer.d \
./src/containers/MetaContainer.d \
./src/containers/MetaHotContainer.d 


# Each subdirectory must supply rules for building sources it contributes
src/containers/%.o: ../src/containers/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


