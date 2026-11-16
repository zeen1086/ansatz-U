ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
BUILD ?= /tmp/lanio2_ansatzU_build
FC = gfortran
FFLAGS ?= -O2 -g -fcheck=bounds -ffree-line-length-none
LEGACYFLAGS = $(FFLAGS) -fdefault-real-8 -fdefault-double-8
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
LIBS ?= -framework Accelerate
else
LIBS ?= -llapack -lblas
endif
JOINT = $(ROOT)/common/numerics.f90 $(ROOT)/4_potential/fit/joint_features.f90 $(ROOT)/4_potential/fit/joint_npy.f90 $(ROOT)/4_potential/fit/joint_model.f90 $(ROOT)/4_potential/fit/joint_svd.f90
JOINT_OBJECTS = $(addprefix $(BUILD)/joint/,numerics.o joint_features.o joint_npy.o joint_model.o joint_svd.o)
LOCAL_SOURCES = $(ROOT)/common/numerics.f90 $(ROOT)/6_RVR/frozen_io.f90 $(ROOT)/4_potential/fit/joint_npy.f90 $(ROOT)/6_RVR/local_vfunc.f90 $(ROOT)/6_RVR/local_rvr.f90
LOCAL_OBJECTS = $(addprefix $(BUILD)/local/,numerics.o frozen_io.o joint_npy.o local_vfunc.o local_rvr.o)
.PHONY: all legacy joint local-rvr clean
all: local-rvr joint legacy
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
joint: $(BUILD)/fit_joint $(BUILD)/eval_joint
$(BUILD)/joint/.built: $(JOINT) | $(BUILD)
	mkdir -p $(BUILD)/joint
	cd $(BUILD)/joint && $(FC) $(FFLAGS) -J. -c $(JOINT)
	touch $@
$(BUILD)/fit_joint: $(ROOT)/4_potential/fit/fit_joint.f90 $(BUILD)/joint/.built
	$(FC) $(FFLAGS) -I$(BUILD)/joint $< $(JOINT_OBJECTS) -o $@ $(LIBS)
$(BUILD)/eval_joint: $(ROOT)/4_potential/fit/eval_joint.f90 $(BUILD)/joint/.built
	$(FC) $(FFLAGS) -I$(BUILD)/joint $< $(JOINT_OBJECTS) -o $@ $(LIBS)
# Independent fitting tool; the V1 spin workflow uses local/cal_Hsk.
legacy: $(BUILD)/fit_new
$(BUILD)/fit_new: $(ROOT)/4_potential/fit/fit_mask.f90 $(ROOT)/4_potential/fit/fit_new.f90 | $(BUILD)
	cd $(BUILD) && $(FC) $(LEGACYFLAGS) -J. $^ -o $@
clean:
	@echo 'Build artifacts are external; remove only your owned BUILD directory explicitly.'
