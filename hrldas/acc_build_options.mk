# CPU is the default; preserve site user_build_options and NetCDF settings.
# Clean before changing modes because .o/.mod names are shared.
ACCELERATOR ?= cpu
ifeq ($(ACCELERATOR),gpu)
  ACC_DEVICE_FLAGS ?= -O2 -fno-fast-math -ffp-contract=off
  # Device math (e.g. real MOD -> fmodf) must be linked for NVPTX explicitly.
  ACC_LDFLAGS = -fopenacc -foffload=nvptx-none -foffload-options=nvptx-none="$(ACC_DEVICE_FLAGS) -lm"
  override CPPFLAGS += -DNOAHMP_OPENACC
  # Cancel both legacy and modern bounds flags from site compiler options.
  override F90FLAGS += -O2 -frecursive -fstack-arrays -fno-realloc-lhs -fno-bounds-check -fcheck=no-bounds -fno-fast-math -ffp-contract=off $(ACC_LDFLAGS)
  override FREESOURCE += -ffree-form -ffree-line-length-none
else ifeq ($(ACCELERATOR),cpu)
  ACC_LDFLAGS =
else
  $(error ACCELERATOR must be cpu or gpu)
endif
