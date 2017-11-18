################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/third_party/zlib/test/example.c \
../src/third_party/zlib/test/infcover.c \
../src/third_party/zlib/test/minigzip.c 

OBJS += \
./src/third_party/zlib/test/example.o \
./src/third_party/zlib/test/infcover.o \
./src/third_party/zlib/test/minigzip.o 

C_DEPS += \
./src/third_party/zlib/test/example.d \
./src/third_party/zlib/test/infcover.d \
./src/third_party/zlib/test/minigzip.d 


# Each subdirectory must supply rules for building sources it contributes
src/third_party/zlib/test/%.o: ../src/third_party/zlib/test/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c99 -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0 -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


