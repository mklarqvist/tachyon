################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../tachyon/third_party/xxhash/xxhash.c 

OBJS += \
./tachyon/third_party/xxhash/xxhash.o 

C_DEPS += \
./tachyon/third_party/xxhash/xxhash.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/third_party/xxhash/%.o: ../tachyon/third_party/xxhash/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c99 -O3 -msse4.2 -g -Wall -c -fmessage-length=0 -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


