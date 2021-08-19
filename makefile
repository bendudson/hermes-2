
BOUT_TOP	?= ../..

TARGET = hermes-2

DIRS = atomicpp

SOURCEC		= hermes-2.cxx div_ops.cxx loadmetric.cxx radiation.cxx neutral-model.cxx \
		  diffusion2d.cxx recycling.cxx full-velocity.cxx mixed.cxx \
		  atomicpp/ImpuritySpecies.cxx atomicpp/Prad.cxx \
		  atomicpp/RateCoefficient.cxx atomicpp/sharedFunctions.cxx

OBJ = $(SOURCEC:%.cxx=%.o)
TARGET ?= $(SOURCEC:%.cxx=%)

ifeq (, $(shell which bout-config))
$(error No bout-config in $$PATH ($(PATH)))
endif

# Use the bout-config script to get compiler, flags
CXX:=$(shell bout-config --cxx)

CFLAGS:=$(shell bout-config --cflags)
LD:=$(shell bout-config --ld)
LDFLAGS:=$(shell bout-config --libs)

$(TARGET): makefile $(OBJ)
	@echo "  Linking" $(TARGET)
	$(LD) -o $(TARGET) $(OBJ) $(LDFLAGS) -lbout++ # -lfmt

%.o: %.cxx
	@echo "  Compiling " $<
	@echo "$<"
	@$(CXX) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJ) $(TARGET)

