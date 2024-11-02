import numpy as np
from parameters import *
import os

def kinetic_energy(particle_velocities):
    kinetic_energy_sum = np.sum(0.5*pow(particle_velocities,2))
    return kinetic_energy_sum

def write_kinetic_energy(filename, simulation_time, kinetic_energy_sum):
    filepath = os.path.join(output_folder_kinetic_energy, filename)
    with open(filepath,"a") as f:
        f.write(f"{simulation_time} {kinetic_energy_sum}\n")   
