################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tachyon/algorithm/permutation/permutation_manager.cpp \
../tachyon/algorithm/permutation/radix_sort_gt.cpp 

OBJS += \
./tachyon/algorithm/permutation/permutation_manager.o \
./tachyon/algorithm/permutation/radix_sort_gt.o 

CPP_DEPS += \
./tachyon/algorithm/permutation/permutation_manager.d \
./tachyon/algorithm/permutation/radix_sort_gt.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/algorithm/permutation/%.o: ../tachyon/algorithm/permutation/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I"../zstd/lib/" -I"../zstd/lib/common/" -I/usr/local/opt/openssl/lib -I/usr/include/openssl/ -I/usr/local/include/ -O3 -msse4.2 -g -Wall -c -fmessage-length=0  -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


