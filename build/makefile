################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm

# All of the sources participating in the build are defined here
-include sources.mk
-include tachyon/third_party/zlib/subdir.mk
-include tachyon/third_party/xxhash/subdir.mk
-include tachyon/support/subdir.mk
-include tachyon/io/vcf/subdir.mk
-include tachyon/io/compression/subdir.mk
-include tachyon/io/bcf/subdir.mk
-include tachyon/io/subdir.mk
-include tachyon/core/base/header/subdir.mk
-include tachyon/core/base/subdir.mk
-include tachyon/core/subdir.mk
-include tachyon/containers/core/subdir.mk
-include tachyon/containers/subdir.mk
-include tachyon/algorithm/permutation/subdir.mk
-include tachyon/algorithm/compression/subdir.mk
-include tachyon/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: tachyon2

# Tool invocations
tachyon2: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: Cross G++ Linker'
	g++ -L/usr/local/Cellar/openssl/1.0.2n/lib -pthread -o "tachyon2" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '
	$(MAKE) --no-print-directory post-build

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(OBJS)$(C_UPPER_DEPS)$(CXX_DEPS)$(C_DEPS)$(CPP_DEPS) tachyon2
	-@echo ' '

post-build:
	-mkdir -p ../bin; mv tachyon2 ../bin/tachyon;
	-@echo ' '

.PHONY: all clean dependents post-build

-include ../makefile.targets
