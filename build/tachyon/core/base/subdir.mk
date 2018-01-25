################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/core/base/MetaCold.cpp \
../tachyon/core/base/MetaEntry.cpp \
../tachyon/core/base/MetaHot.cpp 

OBJS += \
./tachyon/core/base/MetaCold.o \
./tachyon/core/base/MetaEntry.o \
./tachyon/core/base/MetaHot.o 

CPP_DEPS += \
./tachyon/core/base/MetaCold.d \
./tachyon/core/base/MetaEntry.d \
./tachyon/core/base/MetaHot.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/core/base/%.o: ../tachyon/core/base/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/Cellar/openssl/1.0.2n/include -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


