# Project Name
TARGET = DaisyFM

# Sources
CPP_SOURCES = DaisyFM.cpp

# Library Locations
LIBDAISY_DIR = libDaisy
DAISYSP_DIR = DaisySP

BOARD = patch
#APP_TYPE = BOOT_QSPI

# Includes FatFS source files within project.
USE_FATFS = 1

# Size-focused flags
OPT             = -Os
USE_LTO         = 1
EXTRA_CFLAGS   += -flto -ffunction-sections -fdata-sections
EXTRA_CXXFLAGS += -flto -ffunction-sections -fdata-sections \
                  -fno-exceptions -fno-rtti -fno-unwind-tables -fno-asynchronous-unwind-tables -fno-use-cxa-atexit -fno-threadsafe-statics
EXTRA_LDFLAGS  += -flto -Wl,--gc-sections --specs=nano.specs

# Core location, and generic Makefile.
SYSTEM_FILES_DIR = $(LIBDAISY_DIR)/core
include $(SYSTEM_FILES_DIR)/Makefile

