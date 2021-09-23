
BOUT_TOP	= ../..

TARGET = hermes-2

SOURCEC		= src/hermes-2.cxx src/div_ops.cxx src/loadmetric.cxx \
		  src/radiation.cxx src/neutral-model.cxx src/diffusion2d.cxx \
		  src/recycling.cxx src/full-velocity.cxx src/mixed.cxx  \
		  atomicpp/ImpuritySpecies.cxx atomicpp/Prad.cxx \
		  atomicpp/RateCoefficient.cxx atomicpp/sharedFunctions.cxx

# Capture the git version, to be printed in the outputs
GIT_VERSION := $(shell git describe --abbrev=40 --dirty --always --tags)
CXXFLAGS += -DHERMES_REVISION=\"$(GIT_VERSION)\" -I$(PWD)/include -I$(PWD)/atomicpp

include $(BOUT_TOP)/make.config

# Note: Need to create revision.hxx
# Modifying all, $(TARGET) or hermes-3 targets doesn't work,
# but target depends on makefile
makefile: include/revision.hxx
include/revision.hxx: include/revision.hxx.in
	cp $< $@

clean::
	rm -f include/revision.hxx
