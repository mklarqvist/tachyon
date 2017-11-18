################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/third_party/FiniteStateEntropy/programs/bench.c \
../src/third_party/FiniteStateEntropy/programs/commandline.c \
../src/third_party/FiniteStateEntropy/programs/fileio.c \
../src/third_party/FiniteStateEntropy/programs/fseDist.c \
../src/third_party/FiniteStateEntropy/programs/fullbench.c \
../src/third_party/FiniteStateEntropy/programs/fuzzer.c \
../src/third_party/FiniteStateEntropy/programs/fuzzerHuff0.c \
../src/third_party/FiniteStateEntropy/programs/fuzzerU16.c \
../src/third_party/FiniteStateEntropy/programs/probaGenerator.c \
../src/third_party/FiniteStateEntropy/programs/xxhash.c \
../src/third_party/FiniteStateEntropy/programs/zlibh.c 

OBJS += \
./src/third_party/FiniteStateEntropy/programs/bench.o \
./src/third_party/FiniteStateEntropy/programs/commandline.o \
./src/third_party/FiniteStateEntropy/programs/fileio.o \
./src/third_party/FiniteStateEntropy/programs/fseDist.o \
./src/third_party/FiniteStateEntropy/programs/fullbench.o \
./src/third_party/FiniteStateEntropy/programs/fuzzer.o \
./src/third_party/FiniteStateEntropy/programs/fuzzerHuff0.o \
./src/third_party/FiniteStateEntropy/programs/fuzzerU16.o \
./src/third_party/FiniteStateEntropy/programs/probaGenerator.o \
./src/third_party/FiniteStateEntropy/programs/xxhash.o \
./src/third_party/FiniteStateEntropy/programs/zlibh.o 

C_DEPS += \
./src/third_party/FiniteStateEntropy/programs/bench.d \
./src/third_party/FiniteStateEntropy/programs/commandline.d \
./src/third_party/FiniteStateEntropy/programs/fileio.d \
./src/third_party/FiniteStateEntropy/programs/fseDist.d \
./src/third_party/FiniteStateEntropy/programs/fullbench.d \
./src/third_party/FiniteStateEntropy/programs/fuzzer.d \
./src/third_party/FiniteStateEntropy/programs/fuzzerHuff0.d \
./src/third_party/FiniteStateEntropy/programs/fuzzerU16.d \
./src/third_party/FiniteStateEntropy/programs/probaGenerator.d \
./src/third_party/FiniteStateEntropy/programs/xxhash.d \
./src/third_party/FiniteStateEntropy/programs/zlibh.d 


# Each subdirectory must supply rules for building sources it contributes
src/third_party/FiniteStateEntropy/programs/%.o: ../src/third_party/FiniteStateEntropy/programs/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c99 -O3 -march=native -mtune=native -ftree-vectorize -pipe -frename-registers -funroll-loops -g -Wall -c -fmessage-length=0 -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


