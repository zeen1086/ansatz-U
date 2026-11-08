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
JOINT = $(ROOT)/common/numerics.f90 $(ROOT)/4_potential/fit/joint_features.f90 $(ROOT)/4_potential/fit/joint_npy.f90 $(ROOT)/4_potential/fit/joint_model.f90 $(ROOT)/4_potential/fit/joint_svd.f90
JOINT_OBJECTS = $(addprefix $(BUILD)/joint/,numerics.o joint_features.o joint_npy.o joint_model.o joint_svd.o)
RUNTIME = $(filter-out $(ROOT)/common/numerics.f90,$(CORE)) $(ROOT)/6_RVR/joint_rvr.f90 $(ROOT)/7_spin/joint/cal_Hsk.f90
RUNTIME_OBJECTS = $(addprefix $(BUILD)/joint_runtime/,frozen_io.o tensor_rvr.o local_ansatz.o joint_rvr.o cal_Hsk.o)
.PHONY: all legacy mpi joint joint-rvr clean
all: $(BUILD)/ansatzU legacy joint joint-rvr
$(BUILD):
	@case "$(abspath $(BUILD))/" in "$(ROOT)/"*) echo 'BUILD must be outside the source delivery'; exit 1;; esac
	mkdir -p $(BUILD)
$(BUILD)/ansatzU: $(ROOT)/7_spin/run_ansatzU.f90 $(BUILD)/joint_runtime/.built
	$(FC) $(FFLAGS) -I$(BUILD)/joint_runtime -I$(BUILD)/joint $< $(RUNTIME_OBJECTS) $(JOINT_OBJECTS) -o $@ $(LIBS)
joint: $(BUILD)/fit_joint $(BUILD)/eval_joint
joint-rvr: $(BUILD)/ansatzU $(BUILD)/joint/cal_RVR_split $(BUILD)/joint/combine
$(BUILD)/joint_runtime/.built: $(RUNTIME) $(BUILD)/joint/.built
	mkdir -p $(BUILD)/joint_runtime
	cd $(BUILD)/joint_runtime && $(FC) $(FFLAGS) -I$(BUILD)/joint -J. -c $(RUNTIME)
	touch $@
$(BUILD)/joint/cal_RVR_split: $(ROOT)/6_RVR/joint/cal_RVR_split.f90 $(BUILD)/joint_runtime/.built
	$(FC) $(FFLAGS) -I$(BUILD)/joint_runtime -I$(BUILD)/joint $< $(RUNTIME_OBJECTS) $(JOINT_OBJECTS) -o $@ $(LIBS)
$(BUILD)/joint/combine: $(ROOT)/6_RVR/joint/combine.f90 $(BUILD)/joint_runtime/.built
	$(FC) $(FFLAGS) -I$(BUILD)/joint_runtime -I$(BUILD)/joint $< $(RUNTIME_OBJECTS) $(JOINT_OBJECTS) -o $@ $(LIBS)
$(BUILD)/joint/.built: $(JOINT) | $(BUILD)
	mkdir -p $(BUILD)/joint
	cd $(BUILD)/joint && $(FC) $(FFLAGS) -J. -c $(JOINT)
	touch $@
$(BUILD)/fit_joint: $(ROOT)/4_potential/fit/fit_joint.f90 $(BUILD)/joint/.built
	$(FC) $(FFLAGS) -I$(BUILD)/joint $< $(JOINT_OBJECTS) -o $@ $(LIBS)
$(BUILD)/eval_joint: $(ROOT)/4_potential/fit/eval_joint.f90 $(BUILD)/joint/.built
	$(FC) $(FFLAGS) -I$(BUILD)/joint $< $(JOINT_OBJECTS) -o $@ $(LIBS)
legacy: $(BUILD)/fit_new $(BUILD)/rvr_combine $(BUILD)/rearrange_HR $(BUILD)/tensor_Hsk $(BUILD)/tensor_Hsk_AFM $(BUILD)/combine
$(BUILD)/fit_new: $(ROOT)/4_potential/fit/fit_mask.f90 $(ROOT)/4_potential/fit/fit_new.f90 | $(BUILD)
	cd $(BUILD) && $(FC) $(LEGACYFLAGS) -J. $^ -o $@
$(BUILD)/rvr_combine: $(ROOT)/6_RVR/combine_lanio2.f90 | $(BUILD)
	$(FC) $(LEGACYFLAGS) $< -o $@
$(BUILD)/rearrange_HR: $(ROOT)/6_RVR/rearrange_HR_lanio2.f90 | $(BUILD)
	$(FC) $(LEGACYFLAGS) $< -o $@ $(LIBS)
$(BUILD)/tensor_Hsk: $(ROOT)/7_spin/cal_Hsk_lanio2.f90 | $(BUILD)
	$(FC) $(LEGACYFLAGS) $< -o $@
$(BUILD)/tensor_Hsk_AFM: $(ROOT)/7_spin/AFM/cal_Hsk_lanio2.f90 | $(BUILD)
	$(FC) $(LEGACYFLAGS) $< -o $@
$(BUILD)/combine: $(ROOT)/7_spin/combine_lanio2.f90 | $(BUILD)
	$(FC) $(LEGACYFLAGS) $< -o $@
mpi: legacy $(BUILD)/cal_RVR_split $(BUILD)/band $(BUILD)/calc_tot_ene
$(BUILD)/cal_RVR_split: $(ROOT)/6_RVR/cal_RVR_split_lanio2.f90 | $(BUILD)
	$(MPIFC) $(LEGACYFLAGS) $< -o $@
$(BUILD)/band: $(ROOT)/7_spin/band_lanio2.f90 | $(BUILD)
	$(MPIFC) $(LEGACYFLAGS) $< -o $@ $(LIBS)
$(BUILD)/calc_tot_ene: $(ROOT)/7_spin/calc_tot_ene_lanio2.f90 | $(BUILD)
	$(MPIFC) $(LEGACYFLAGS) $< -o $@ $(LIBS)
clean:
	@echo 'Build artifacts are external; remove only your owned BUILD directory explicitly.'
