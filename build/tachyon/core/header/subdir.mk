################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/core/header/tachyon_header.cpp 

OBJS += \
./tachyon/core/header/tachyon_header.o 

CPP_DEPS += \
./tachyon/core/header/tachyon_header.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/core/header/%.o: ../tachyon/core/header/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


