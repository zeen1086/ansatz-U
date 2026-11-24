ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
BUILD ?= /tmp/lanio2_ansatzU_build
FC = gfortran
FFLAGS ?= -O2 -g -fcheck=bounds -ffree-line-length-none
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
LIBS ?= -framework Accelerate
else
LIBS ?= -llapack -lblas
endif
LOCAL_SOURCES = $(ROOT)/common/numerics.f90 $(ROOT)/6_RVR/frozen_io.f90 $(ROOT)/common/npy_io.f90 $(ROOT)/6_RVR/local_vfunc.f90 $(ROOT)/6_RVR/local_rvr.f90
LOCAL_OBJECTS = $(addprefix $(BUILD)/local/,numerics.o frozen_io.o npy_io.o local_vfunc.o local_rvr.o)
.PHONY: all local-rvr clean
all: local-rvr
local-rvr: $(BUILD)/local/rvr $(BUILD)/local/cal_Hsk
$(BUILD)/local/.built: $(LOCAL_SOURCES) | $(BUILD)
	mkdir -p $(BUILD)/local
	cd $(BUILD)/local && $(FC) $(FFLAGS) -c $(LOCAL_SOURCES)
	touch $@
$(BUILD)/local/rvr: $(ROOT)/6_RVR/rvr.f90 $(BUILD)/local/.built
	$(FC) $(FFLAGS) -I$(BUILD)/local $< $(LOCAL_OBJECTS) -o $@ $(LIBS)
$(BUILD)/local/cal_Hsk: $(ROOT)/7_spin/cal_Hsk.f90 $(BUILD)/local/.built
	$(FC) $(FFLAGS) -I$(BUILD)/local $< $(LOCAL_OBJECTS) -o $@ $(LIBS)
$(BUILD):
	@case "$(abspath $(BUILD))/" in "$(ROOT)/"*) echo 'BUILD must be outside the source delivery'; exit 1;; esac
	mkdir -p $(BUILD)
clean:
	@echo 'Build artifacts are external; remove only your owned BUILD directory explicitly.'
