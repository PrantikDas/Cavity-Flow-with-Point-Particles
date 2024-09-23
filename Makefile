# Variables
PYTHON := python3
SCRIPT := src/fast_poisson_parallel.py
FLUID_OUTPUT_DIR := output/fluid_information
PARTICLE_OUTPUT_DIR := output/particle_information
VTK_FILES := $(FLUID_OUTPUT_DIR)/*.vtk
PARTICLE_FILES := $(PARTICLE_OUTPUT_DIR)/*.liggghts

# Default target: Run the Python script
.PHONY: all
all: run

# Run the Python script
run:
	$(PYTHON) $(SCRIPT)

# Clean target: Remove the generated output files
.PHONY: clean
clean:
	rm -f $(VTK_FILES) $(PARTICLE_FILES)

# Clean everything: Remove directories as well
.PHONY: clean-all
clean-all:
	rm -rf $(FLUID_OUTPUT_DIR) $(PARTICLE_OUTPUT_DIR) Residuals.txt
	
# Re-run the simulation after cleaning
rerun: clean run
