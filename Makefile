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
OPT            = -Os
#USE_LTO        = 1
EXTRA_CFLAGS   += -ffunction-sections -fdata-sections
EXTRA_LDFLAGS  += -Wl,--gc-sections
#USE_NEWLIB_NANO = 1

# Core location, and generic Makefile.
SYSTEM_FILES_DIR = $(LIBDAISY_DIR)/core
include $(SYSTEM_FILES_DIR)/Makefile

