################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/third_party/FiniteStateEntropy/lib/Archives/hufx6_decompress.c 

OBJS += \
./src/third_party/FiniteStateEntropy/lib/Archives/hufx6_decompress.o 

C_DEPS += \
./src/third_party/FiniteStateEntropy/lib/Archives/hufx6_decompress.d 


# Each subdirectory must supply rules for building sources it contributes
src/third_party/FiniteStateEntropy/lib/Archives/%.o: ../src/third_party/FiniteStateEntropy/lib/Archives/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c99 -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0 -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


