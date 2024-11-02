# Define constants
N_POINTS_X = 40
N_POINTS_Y = 40
DOMAIN_SIZE_X = 1.0
DOMAIN_SIZE_Y = 1.0
N_ITERATIONS = 100000
TIME_STEP_LENGTH = 0.001
KINEMATIC_VISCOSITY = 0.1
DENSITY = 1.0
HORIZONTAL_VELOCITY_TOP = 5.0
N_PARTICLES = 10  # Number of particles
N_PRESSURE_POISSON_ITERATIONS = 200
STABILITY_SAFETY_FACTOR = 0.5
particle_radius = 0.001
DOMAIN_HEIGHT = 0.001
vtk_interval = 100  # VTK output interval (in number of iterations)
output_folder_fluid = "output/fluid_information"  # Output folder for files
output_folder_particle = "output/particle_information" # Output folder for files
output_folder_kinetic_energy = "postProcess/Kinetic_Energy"
output_folder_residual = "postProcess/Residual"