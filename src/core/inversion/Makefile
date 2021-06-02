TARGET := isola

PARAMETERS_OBJ  := parameters.o
COMMON_VARS_OBJ := common_vars.o
ISOLA15_OBJ     := isola15.o

OBJS  := $(PARAMETERS_OBJ) $(COMMON_VARS_OBJ) $(ISOLA15_OBJ)

ifeq ($(MODEL),openmp)
	FC := gfortran
	FFLAGS := -O3 -fopenmp -foffload=disable
else ifeq ($(MODEL), openacc)
	FC := /opt/nvidia/hpc_sdk/Linux_x86_64/20.9/compilers/bin/nvfortran
	FFLAGS := -O2 -mcmodel=medium -acc=gpu,multicore -gpu=cc35,nofma -Minfo=accel
endif

all: clean $(ISOLA15_OBJ)
	$(FC) $(FFLAGS) $(OBJS) -o $(TARGET)

ifndef MODEL
	$(error "Please specify your preferred programming model: make MODEL=[openacc,openmp]")
endif

$(ISOLA15_OBJ): $(COMMON_VARS_OBJ)

$(COMMON_VARS_OBJ): $(PARAMETERS_OBJ)

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $<

ifndef MODEL
	$(error "Please specify your preferred programming model: make MODEL=[openacc,openmp]")
endif

clean:
	rm -f $(TARGET) $(OBJS) *.mod

