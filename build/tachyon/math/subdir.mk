################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/math/fisher_math.cpp 

OBJS += \
./tachyon/math/fisher_math.o 

CPP_DEPS += \
./tachyon/math/fisher_math.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/math/%.o: ../tachyon/math/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I"../zstd/lib/" -I"../zstd/lib/common/" -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


