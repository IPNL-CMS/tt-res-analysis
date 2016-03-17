# Check if the installation path of the package is provided
ifeq ($(TTRES_ANALYSIS_INSTALL), )
  $(error Mandatory environment variable TTRES_ANALYSIS_INSTALL is not set)
endif

# Check if the installation path of PECFwk is provided
ifeq ($(PEC_FWK_INSTALL), )
  $(error Mandatory environment variable PEC_FWK_INSTALL is not set)
endif

# Make sure ROOT is available
ifeq ($(shell which root-config), )
  $(error ROOT installation is not found)
endif


# Flags to control compilation and linking
CC = g++
INCLUDE = -I$(TTRES_ANALYSIS_INSTALL)/include -I$(PEC_FWK_INSTALL)/include -I$(shell root-config --incdir)
OPFLAGS = -O2
CFLAGS = -Wall -Wextra -Wno-unused-function -fPIC -std=c++14 $(INCLUDE) $(OPFLAGS)


# Sources, object files, and library
LIB_DIR := lib
LIB_NAME := libTTRes.so
LIB_PATH := $(LIB_DIR)/$(LIB_NAME)

SOURCE_DIR := src
OBJ_DIR := obj
OBJECTS := $(shell for f in `find $(SOURCE_DIR) -regex ".*\.cpp"`; \
	do echo $(OBJ_DIR)/`basename $$f .cpp`.o; done)

vpath %.cpp $(SOURCE_DIR)


.PHONY: clean


# Building rules
all: $(LIB_PATH)


$(LIB_PATH): $(OBJECTS)
	@ mkdir -p $(LIB_DIR)
	@ rm -f $@
	$(CC) -shared -Wl,-soname,$(LIB_NAME) -o $@ $^


$(OBJ_DIR)/%.o: $(SOURCE_DIR)/%.cpp
	@ mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@


clean:
	rm -rf $(OBJ_DIR)
