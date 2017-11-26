################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/algorithm/encryption/openssl_aes.cpp 

CPP_DEPS += \
./src/algorithm/encryption/openssl_aes.d 


# Each subdirectory must supply rules for building sources it contributes
src/algorithm/encryption/openssl_aes: ../src/algorithm/encryption/openssl_aes.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c99 -I/Users/mk21/homebrew/opt/openssl/include/ -O3 -march=native -mtune=native -g -Wall -c -fmessage-length=0 -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"src/algorithm/encryption/openssl_aes.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


