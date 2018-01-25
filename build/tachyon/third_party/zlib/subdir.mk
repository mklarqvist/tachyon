################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../tachyon/third_party/zlib/adler32.c \
../tachyon/third_party/zlib/compress.c \
../tachyon/third_party/zlib/crc32.c \
../tachyon/third_party/zlib/deflate.c \
../tachyon/third_party/zlib/gzclose.c \
../tachyon/third_party/zlib/gzlib.c \
../tachyon/third_party/zlib/gzread.c \
../tachyon/third_party/zlib/gzwrite.c \
../tachyon/third_party/zlib/infback.c \
../tachyon/third_party/zlib/inffast.c \
../tachyon/third_party/zlib/inflate.c \
../tachyon/third_party/zlib/inftrees.c \
../tachyon/third_party/zlib/trees.c \
../tachyon/third_party/zlib/uncompr.c \
../tachyon/third_party/zlib/zutil.c 

OBJS += \
./tachyon/third_party/zlib/adler32.o \
./tachyon/third_party/zlib/compress.o \
./tachyon/third_party/zlib/crc32.o \
./tachyon/third_party/zlib/deflate.o \
./tachyon/third_party/zlib/gzclose.o \
./tachyon/third_party/zlib/gzlib.o \
./tachyon/third_party/zlib/gzread.o \
./tachyon/third_party/zlib/gzwrite.o \
./tachyon/third_party/zlib/infback.o \
./tachyon/third_party/zlib/inffast.o \
./tachyon/third_party/zlib/inflate.o \
./tachyon/third_party/zlib/inftrees.o \
./tachyon/third_party/zlib/trees.o \
./tachyon/third_party/zlib/uncompr.o \
./tachyon/third_party/zlib/zutil.o 

C_DEPS += \
./tachyon/third_party/zlib/adler32.d \
./tachyon/third_party/zlib/compress.d \
./tachyon/third_party/zlib/crc32.d \
./tachyon/third_party/zlib/deflate.d \
./tachyon/third_party/zlib/gzclose.d \
./tachyon/third_party/zlib/gzlib.d \
./tachyon/third_party/zlib/gzread.d \
./tachyon/third_party/zlib/gzwrite.d \
./tachyon/third_party/zlib/infback.d \
./tachyon/third_party/zlib/inffast.d \
./tachyon/third_party/zlib/inflate.d \
./tachyon/third_party/zlib/inftrees.d \
./tachyon/third_party/zlib/trees.d \
./tachyon/third_party/zlib/uncompr.d \
./tachyon/third_party/zlib/zutil.d 


# Each subdirectory must supply rules for building sources it contributes
tachyon/third_party/zlib/%.o: ../tachyon/third_party/zlib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -std=c99 -O3 -msse4.2 -g -Wall -c -fmessage-length=0 -DVERSION=\"$(GIT_VERSION)\" -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


