################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/io/compression/BGZFController.cpp \
../tachyon/io/compression/TGZFController.cpp \
../tachyon/io/compression/TGZFControllerStream.cpp 

OBJS += \
./tachyon/io/compression/BGZFController.o \
./tachyon/io/compression/TGZFController.o \
./tachyon/io/compression/TGZFControllerStream.o 

CPP_DEPS += \
./tachyon/io/compression/BGZFController.d \
./tachyon/io/compression/TGZFController.d \
./tachyon/io/compression/TGZFControllerStream.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/io/compression/%.o: ../tachyon/io/compression/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/Cellar/openssl/1.0.2n/include -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


