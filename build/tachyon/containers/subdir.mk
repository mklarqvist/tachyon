################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/containers/DataBlock.cpp \
../tachyon/containers/DataContainer.cpp \
../tachyon/containers/MetaColdContainer.cpp \
../tachyon/containers/MetaContainer.cpp \
../tachyon/containers/MetaHotContainer.cpp 

OBJS += \
./tachyon/containers/DataBlock.o \
./tachyon/containers/DataContainer.o \
./tachyon/containers/MetaColdContainer.o \
./tachyon/containers/MetaContainer.o \
./tachyon/containers/MetaHotContainer.o 

CPP_DEPS += \
./tachyon/containers/DataBlock.d \
./tachyon/containers/DataContainer.d \
./tachyon/containers/MetaColdContainer.d \
./tachyon/containers/MetaContainer.d \
./tachyon/containers/MetaHotContainer.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/containers/%.o: ../tachyon/containers/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/Cellar/openssl/1.0.2n/include -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


