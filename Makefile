# Check if the installation path of this package is set
ifeq ($(TTRES_ANALYSIS_INSTALL), )
  $(error Mandatory environment variable TTRES_ANALYSIS_INSTALL is not set)
endif

# Check dependencies
ifeq ($(MENSURA_INSTALL), )
  $(error Mandatory environment variable MENSURA_INSTALL is not set)
endif

ifeq ($(shell which root-config), )
  $(error ROOT installation is not found)
endif

ifeq ($(BOOST_ROOT), )
  $(error Mandatory environment variable BOOST_ROOT is not set)
endif

ifeq ($(shell which lhapdf-config), )
  $(error LHAPDF installation is not found)
endif


# Flags to control compilation and linking
CC := g++
INCLUDE := -I$(TTRES_ANALYSIS_INSTALL)/include -I$(MENSURA_INSTALL)/include \
  -I$(shell root-config --incdir) -I$(shell lhapdf-config --incdir) -I$(BOOST_ROOT)/include
OPFLAGS := -O2
CFLAGS := -Wall -Wextra -Wno-unused-function -fPIC -std=c++14 $(INCLUDE) $(OPFLAGS)
LDFLAGS := -L$(MENSURA_INSTALL)/lib -lmensura -lPECReader $(shell root-config --libs) \
  $(shell lhapdf-config --libs)
LDFLAGS_BIN := -L$(TTRES_ANALYSIS_INSTALL)/lib -lTTRes \
  -L$(BOOST_ROOT)/lib -lboost_program_options


# Sources, object files, executables and library
LIB_DIR := lib
LIB_NAME := libTTRes.so
LIB_PATH := $(LIB_DIR)/$(LIB_NAME)

SOURCE_DIR := src
OBJ_DIR := obj
OBJECTS := $(shell for f in `find $(SOURCE_DIR) -regex ".*\.cpp"`; \
	do echo $(OBJ_DIR)/`basename $$f .cpp`.o; done)

BIN_DIR := bin
BIN_SRC_DIR := prog
PROGS := htt-tuples

vpath %.cpp $(SOURCE_DIR) $(BIN_SRC_DIR)


.PHONY: clean


# Building rules
all: $(LIB_PATH) $(addprefix $(BIN_DIR)/,$(PROGS))


$(LIB_PATH): $(OBJECTS)
	@ mkdir -p $(LIB_DIR)
	@ rm -f $@
	$(CC) -shared -Wl,-soname,$(LIB_NAME) -o $@ $^ $(LDFLAGS)


$(BIN_DIR)/%: $(OBJ_DIR)/%.o $(LIB_PATH)
	@ mkdir -p $(BIN_DIR)
	g++ $^ $(CFLAGS) $(LDFLAGS_BIN) -o $@


$(OBJ_DIR)/%.o: %.cpp
	@ mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@


clean:
	rm -rf $(OBJ_DIR)
