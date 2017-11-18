################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
S_UPPER_SRCS += \
../src/third_party/zlib/contrib/amd64/crc32-pclmul_asm.S 

OBJS += \
./src/third_party/zlib/contrib/amd64/crc32-pclmul_asm.o 


# Each subdirectory must supply rules for building sources it contributes
src/third_party/zlib/contrib/amd64/%.o: ../src/third_party/zlib/contrib/amd64/%.S
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Assembler'
	as  -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


