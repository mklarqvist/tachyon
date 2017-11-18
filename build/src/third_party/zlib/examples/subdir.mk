################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/third_party/zlib/examples/enough.c \
../src/third_party/zlib/examples/fitblk.c \
../src/third_party/zlib/examples/gun.c \
../src/third_party/zlib/examples/gzappend.c \
../src/third_party/zlib/examples/gzjoin.c \
../src/third_party/zlib/examples/gzlog.c \
../src/third_party/zlib/examples/zpipe.c \
../src/third_party/zlib/examples/zran.c 

OBJS += \
./src/third_party/zlib/examples/enough.o \
./src/third_party/zlib/examples/fitblk.o \
./src/third_party/zlib/examples/gun.o \
./src/third_party/zlib/examples/gzappend.o \
./src/third_party/zlib/examples/gzjoin.o \
./src/third_party/zlib/examples/gzlog.o \
./src/third_party/zlib/examples/zpipe.o \
./src/third_party/zlib/examples/zran.o 

C_DEPS += \
./src/third_party/zlib/examples/enough.d \
./src/third_party/zlib/examples/fitblk.d \
./src/third_party/zlib/examples/gun.d \
./src/third_party/zlib/examples/gzappend.d \
./src/third_party/zlib/examples/gzjoin.d \
./src/third_party/zlib/examples/gzlog.d \
./src/third_party/zlib/examples/zpipe.d \
./src/third_party/zlib/examples/zran.d 


# Each subdirectory must supply rules for building sources it contributes
src/third_party/zlib/examples/%.o: ../src/third_party/zlib/examples/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c99 -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0 -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


