################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/io/bcf/BCFEntry.cpp \
../tachyon/io/bcf/BCFReader.cpp 

OBJS += \
./tachyon/io/bcf/BCFEntry.o \
./tachyon/io/bcf/BCFReader.o 

CPP_DEPS += \
./tachyon/io/bcf/BCFEntry.d \
./tachyon/io/bcf/BCFReader.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/io/bcf/%.o: ../tachyon/io/bcf/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


