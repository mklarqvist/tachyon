################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/main.cpp \
../tachyon/variant_importer.cpp \
../tachyon/variant_reader.cpp 

OBJS += \
./tachyon/main.o \
./tachyon/variant_importer.o \
./tachyon/variant_reader.o 

CPP_DEPS += \
./tachyon/main.d \
./tachyon/variant_importer.d \
./tachyon/variant_reader.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/%.o: ../tachyon/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


