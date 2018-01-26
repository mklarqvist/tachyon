################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/core/base/meta_cold.cpp \
../tachyon/core/base/meta_entry.cpp \
../tachyon/core/base/meta_hot.cpp 

OBJS += \
./tachyon/core/base/meta_cold.o \
./tachyon/core/base/meta_entry.o \
./tachyon/core/base/meta_hot.o 

CPP_DEPS += \
./tachyon/core/base/meta_cold.d \
./tachyon/core/base/meta_entry.d \
./tachyon/core/base/meta_hot.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/core/base/%.o: ../tachyon/core/base/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/Cellar/openssl/1.0.2n/include -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


