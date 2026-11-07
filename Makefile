ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
BUILD ?= /tmp/lanio2_ansatzU_build
FC = gfortran
MPIFC ?= mpif90
FFLAGS ?= -O2 -g -fcheck=bounds -ffree-line-length-none
LEGACYFLAGS = $(FFLAGS) -fdefault-real-8 -fdefault-double-8
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
LIBS ?= -framework Accelerate
else
LIBS ?= -llapack -lblas
endif
CORE = $(ROOT)/common/numerics.f90 $(ROOT)/6_RVR/frozen_io.f90 $(ROOT)/6_RVR/tensor_rvr.f90 $(ROOT)/7_spin/local_ansatz.f90
.PHONY: all legacy mpi clean
all: $(BUILD)/ansatzU legacy
$(BUILD):
	@case "$(abspath $(BUILD))/" in "$(ROOT)/"*) echo 'BUILD must be outside the source delivery'; exit 1;; esac
	mkdir -p $(BUILD)
$(BUILD)/ansatzU: $(CORE) $(ROOT)/7_spin/run_ansatzU.f90 | $(BUILD)
	cd $(BUILD) && $(FC) $(FFLAGS) -J. $(CORE) $(ROOT)/7_spin/run_ansatzU.f90 -o ansatzU $(LIBS)
legacy: $(BUILD)/fit_new $(BUILD)/rvr_combine $(BUILD)/rearrange_HR $(BUILD)/tensor_Hsk $(BUILD)/tensor_Hsk_AFM $(BUILD)/combine
$(BUILD)/fit_new: $(ROOT)/4_potential/fit/fit_mask.f90 $(ROOT)/4_potential/fit/fit_new.f90 | $(BUILD)
	cd $(BUILD) && $(FC) $(LEGACYFLAGS) -J. $^ -o $@
$(BUILD)/rvr_combine: $(ROOT)/6_RVR/combine.f90 | $(BUILD)
	$(FC) $(LEGACYFLAGS) $< -o $@
$(BUILD)/rearrange_HR: $(ROOT)/6_RVR/rearrange_HR.f90 | $(BUILD)
	$(FC) $(LEGACYFLAGS) $< -o $@ $(LIBS)
$(BUILD)/tensor_Hsk: $(ROOT)/7_spin/cal_Hsk.f90 | $(BUILD)
	$(FC) $(LEGACYFLAGS) $< -o $@
$(BUILD)/tensor_Hsk_AFM: $(ROOT)/7_spin/AFM/cal_Hsk.f90 | $(BUILD)
	$(FC) $(LEGACYFLAGS) $< -o $@
$(BUILD)/combine: $(ROOT)/7_spin/combine.f90 | $(BUILD)
	$(FC) $(LEGACYFLAGS) $< -o $@
mpi: legacy $(BUILD)/cal_RVR_split $(BUILD)/band $(BUILD)/calc_tot_ene
$(BUILD)/cal_RVR_split: $(ROOT)/6_RVR/cal_RVR_split.f90 | $(BUILD)
	$(MPIFC) $(LEGACYFLAGS) $< -o $@
$(BUILD)/band: $(ROOT)/7_spin/band.f90 | $(BUILD)
	$(MPIFC) $(LEGACYFLAGS) $< -o $@ $(LIBS)
$(BUILD)/calc_tot_ene: $(ROOT)/7_spin/calc_tot_ene.f90 | $(BUILD)
	$(MPIFC) $(LEGACYFLAGS) $< -o $@ $(LIBS)
clean:
	@echo 'Build artifacts are external; remove only your owned BUILD directory explicitly.'
