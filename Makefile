DEP_DIR:=./deps
SRC_DIR:=src
LIB_DIR:=lib
BIN_DIR:=bin
OBJ_DIR:=obj

INC_DIR:=include
CWD:=$(shell pwd)
CC ?= gcc-11
CXX ?= g++-11
PKG_CONFIG ?= pkg-config

SFX :=
EXE:=bgtf$(SFX)

all: $(BIN_DIR)/$(EXE)

# dependencies
HTSLIB_DIR:=deps/htslib
LIBDEFLATE_DIR:=deps/libdeflate
INCLUDE_FLAGS :=-I$(CWD)/$(INC_DIR) -isystem $(CWD)/$(INC_DIR) -I. -I$(CWD)/$(SRC_DIR) 
LD_LIB_DIR_FLAGS := -L$(CWD)/$(LIB_DIR)
# LD_LIB_FLAGS := 
LD_STATIC_LIB_FLAGS := $(CWD)/$(LIB_DIR)/libhts.a -lz -lbz2 -llzma -lcurl

CFLAGS = -O3 -Werror=return-type -w

LD_STATIC_LIB_DEPS := -lpthread -lm

LIB_DEPS += $(LIB_DIR)/libhts.a
LIB_DEPS += $(LIB_DIR)/libdeflate.a 

# Travis needs -latomic for all builds *but* GCC on Mac
ifeq ($(strip $(shell $(CXX) -latomic /dev/null -o/dev/null 2>&1 | grep latomic | wc -l)), 0)
    # Use -latomic if the compiler doesn't complain about it
    LD_LIB_FLAGS += -latomic
endif

# These libs need to come after libdw if used, because libdw depends on them
LD_LIB_FLAGS += -ldl -llzma -lbz2

# When building statically, we need to tell the linker not to bail if it sees multiple definitions.
# libc on e.g. our Jenkins host does not define malloc as weak, so other mallocs can't override it in a static build.
# TODO: Why did this problem only begin to happen when libvw was added?
STATIC_FLAGS=-static -static-libgcc -Wl,--allow-multiple-definition 

# Use a fake rule to build .d files, so we don't complain if they don't exist.
$(OBJ_DIR)/%.d: ;
# Don't delete them.
.PRECIOUS: $(OBJ_DIR)/%.d 

# Use no implicit rules
.SUFFIXES:

LIB_DEPS =
DEPS = $(LIB_DIR)/libhts.a $(LIB_DIR)/libdeflate.a 

LIBBGTF_DEPS = $(OBJ) $(DEPS)
EXE_DEPS = $(OBJ) $(OBJ_DIR)/main.o $(DEPS)
.PHONY: build clean 

OBJ = $(filter-out $(OBJ_DIR)/main.o, $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(wildcard $(SRC_DIR)/*.c)))

$(OBJ_DIR)/main.o:$(SRC_DIR)/main.c
		$(CC) -c $(CFLAGS) $(INCLUDES) $(INCLUDE_FLAGS) $< -o $@

$(OBJ_DIR)/utils.o:$(SRC_DIR)/utils.c 
		$(CC) -c $(CFLAGS) $(INCLUDES) $(INCLUDE_FLAGS) $< -o $@
$(OBJ_DIR)/bgtf.o: $(SRC_DIR)/bgtf.c
		$(CC) -c $(CFLAGS) $(INCLUDES) $(INCLUDE_FLAGS) $< -o $@
$(OBJ_DIR)/hashtable.o: $(SRC_DIR)/hashtable.c
		$(CC) -c $(CFLAGS) $(INCLUDES) $(INCLUDE_FLAGS) $< -o $@
$(OBJ_DIR)/rtree.o: $(SRC_DIR)/rtree.c
		$(CC) -c $(CFLAGS) $(INCLUDES) $(INCLUDE_FLAGS) $< -o $@
$(OBJ_DIR)/parser.o: $(SRC_DIR)/parser.c
		$(CC) -c $(CFLAGS) $(INCLUDES) $(INCLUDE_FLAGS) $< -o $@
$(OBJ_DIR)/bgtfIO.o: $(SRC_DIR)/bgtfIO.c
		$(CC) -c $(CFLAGS) $(INCLUDES) $(INCLUDE_FLAGS) $< -o $@
$(OBJ_DIR)/tpool.o: $(SRC_DIR)/tpool.c
		$(CC) -c $(CFLAGS) $(INCLUDES) $(INCLUDE_FLAGS) $< -o $@
$(OBJ_DIR)/list.o: $(SRC_DIR)/list.c
		$(CC) -c $(CFLAGS) $(INCLUDES) $(INCLUDE_FLAGS) $< -o $@
$(OBJ_DIR)/molcount.o: $(SRC_DIR)/molcount.c
		$(CC) -c $(CFLAGS) $(INCLUDES) $(INCLUDE_FLAGS) $< -o $@
$(OBJ_DIR)/sparse.o: $(SRC_DIR)/sparse.c
		$(CC) -c $(CFLAGS) $(INCLUDES) $(INCLUDE_FLAGS) $< -o $@

$(LIB_DIR)/libdeflate.a: $(LIBDEFLATE_DIR)/*.h $(LIBDEFLATE_DIR)/lib/*.h $(LIBDEFLATE_DIR)/lib/*/*.h $(LIBDEFLATE_DIR)/lib/*.c $(LIBDEFLATE_DIR)/lib/*/*.c
	+. ./source_me.sh && cd $(LIBDEFLATE_DIR) && V=1 $(MAKE) $(FILTER) && cp libdeflate.a $(CWD)/$(LIB_DIR)

$(LIB_DIR)/libhts.a: $(HTSLIB_DIR)/*.c $(HTSLIB_DIR)/htslib/*.h
	+. ./source_me.sh && cd $(HTSLIB_DIR) && $(MAKE) lib-static
	cp $(CWD)/$(HTSLIB_DIR)/libhts.a $(CWD)/$(LIB_DIR)

$(LIB_DIR)/libbgtf.a: $(LIBBGTF_DEPS)
	ar rs $@ $(OBJ)

$(BIN_DIR)/$(EXE): $(LIB_DIR)/libbgtf.a $(EXE_DEPS)
	$(CC) $(INCLUDE_FLAGS) $(CFLAGS) $(LD_STATIC_LIB_FLAGS) -o $(BIN_DIR)/$(EXE) $(OBJ_DIR)/main.o $(LD_STATIC_LIB_DEPS) $(LIB_DEPS) $(OBJ) $(DEPS)

clean: 
	$(RM) -r $(HTSLIB_DIR)/*.o
	$(RM) -r $(BIN_DIR)/*
	$(RM) -r $(LIB_DIR)/*.a $(LIB_DIR)/*.so $(LIB_DIR)/*.dylib