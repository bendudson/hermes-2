
BOUT_TOP	= ../..

BOUT_GIT_COMMIT:=$(shell git --git-dir=$(BOUT_TOP)/.git rev-parse --short HEAD >> BOUT_commit)

TARGET = hermes-2

DIRS = atomicpp

SOURCEC		= hermes-2.cxx div_ops.cxx loadmetric.cxx radiation.cxx neutral-model.cxx diffusion2d.cxx recycling.cxx full-velocity.cxx mixed.cxx

include $(BOUT_TOP)/make.config
