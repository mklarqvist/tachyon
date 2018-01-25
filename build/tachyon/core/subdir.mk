################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/core/ImportWriter.cpp \
../tachyon/core/Importer.cpp 

OBJS += \
./tachyon/core/ImportWriter.o \
./tachyon/core/Importer.o 

CPP_DEPS += \
./tachyon/core/ImportWriter.d \
./tachyon/core/Importer.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/core/%.o: ../tachyon/core/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/Cellar/openssl/1.0.2n/include -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


