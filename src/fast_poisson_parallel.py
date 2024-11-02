import matplotlib.pyplot as plt
import os
import numpy as np
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
from vtk_write import *
from parameters import *
from fluid_discretization_schemes import *
from particle_motion import *
from energy_budget import *

# Ensure the output folder exists, if not these command creates them
os.makedirs(output_folder_fluid, exist_ok=True)
os.makedirs(output_folder_particle, exist_ok=True)
os.makedirs(output_folder_residual, exist_ok=True)
os.makedirs(output_folder_kinetic_energy, exist_ok=True)

def main():
    element_length_x = DOMAIN_SIZE_X / (N_POINTS_X - 1)
    element_length_y = DOMAIN_SIZE_Y / (N_POINTS_Y - 1)
    x = np.linspace(0.0, DOMAIN_SIZE_X, N_POINTS_X)
    y = np.linspace(0.0, DOMAIN_SIZE_Y, N_POINTS_Y)
    X, Y = np.meshgrid(x, y)
    
    u_prev = np.zeros_like(X)
    v_prev = np.zeros_like(X)
    p_prev = np.zeros_like(X)

    particle_positions = np.random.rand(N_PARTICLES, 2) * min(DOMAIN_SIZE_X, DOMAIN_SIZE_Y)
    particle_velocities = np.zeros((N_PARTICLES, 2))  # Initialize velocities

    maximum_possible_time_step_length = 0.5 * (min(element_length_x, element_length_y))**2 / KINEMATIC_VISCOSITY
    if TIME_STEP_LENGTH > STABILITY_SAFETY_FACTOR * maximum_possible_time_step_length:
        raise RuntimeError("Stability is not guaranteed")

    for iteration in tqdm(range(N_ITERATIONS)):
    
        u_0 = u_prev.copy()
        v_0 = v_prev.copy()
        p_0 = p_prev.copy()
        
        d_u_prev__d_x = central_difference_x(u_prev, element_length_x)
        d_u_prev__d_y = central_difference_y(u_prev, element_length_y)
        d_v_prev__d_x = central_difference_x(v_prev, element_length_x)
        d_v_prev__d_y = central_difference_y(v_prev, element_length_y)
        laplace__u_prev = laplace(u_prev, element_length_x)
        laplace__v_prev = laplace(v_prev, element_length_y)

        u_tent = u_prev + TIME_STEP_LENGTH * (- (u_prev * d_u_prev__d_x + v_prev * d_u_prev__d_y) + KINEMATIC_VISCOSITY * laplace__u_prev)
        v_tent = v_prev + TIME_STEP_LENGTH * (- (u_prev * d_v_prev__d_x + v_prev * d_v_prev__d_y) + KINEMATIC_VISCOSITY * laplace__v_prev)

        # Boundary conditions
        u_tent[0, :] = 0.0
        u_tent[:, 0] = 0.0
        u_tent[:, -1] = 0.0
        u_tent[-1, :] = HORIZONTAL_VELOCITY_TOP
        v_tent[0, :] = 0.0
        v_tent[:, 0] = 0.0
        v_tent[:, -1] = 0.0
        v_tent[-1, :] = 0.0

        d_u_tent__d_x = central_difference_x(u_tent, element_length_x)
        d_v_tent__d_y = central_difference_y(v_tent, element_length_y)

        rhs = DENSITY / TIME_STEP_LENGTH * (d_u_tent__d_x + d_v_tent__d_y)

        for _ in range(N_PRESSURE_POISSON_ITERATIONS):
            p_next = np.zeros_like(p_prev)
            prefactor = (pow(element_length_x,2)*pow(element_length_y,2))/(2*(pow(element_length_x,2) + pow(element_length_y,2)))
            p_next[1:-1, 1:-1] =  prefactor * (((p_prev[1:-1, 0:-2] + p_prev[1:-1, 2:])/(pow(element_length_y,2))) + ((p_prev[0:-2, 1:-1] + p_prev[2:, 1:-1]) / pow(element_length_x,2)) - rhs[1:-1, 1:-1])
            p_prev = p_next

        d_p_next__d_x = central_difference_x(p_next, element_length_x)
        d_p_next__d_y = central_difference_y(p_next, element_length_y)

        u_next = u_tent - TIME_STEP_LENGTH / DENSITY * d_p_next__d_x
        v_next = v_tent - TIME_STEP_LENGTH / DENSITY * d_p_next__d_y

        particle_positions, particle_velocities = parallel_particle_update(u_next, v_next, particle_positions, X, Y, element_length_x, element_length_y, TIME_STEP_LENGTH)
        volume_fraction = calculate_cell_volume_fraction(particle_positions, X, Y, element_length_x, element_length_y, particle_radius)
        
        # Correct the residue code part     
        residue_u = np.max(np.abs(u_next - u_0) / np.maximum(u_0,1e-5))
        residue_v = np.max(np.abs(v_next - v_0) / np.maximum(v_0,1e-5))
        residue_p = np.max(np.abs(p_next - p_0) / np.maximum(p_0,1e-5))

        u_prev = u_next
        v_prev = v_next
        p_prev = p_next
        
        simulation_time = iteration * TIME_STEP_LENGTH

        if iteration % vtk_interval == 0:
            write_vtk(f"output_{iteration:04d}.vtk", X, Y, u_prev, v_prev, p_prev, volume_fraction)
            write_particles_liggghts(f"particles_{iteration:04d}.liggghts", iteration, particle_positions, particle_velocities)
            write_residuals(f"Residuals.txt", f"{simulation_time:.6f}", residue_u, residue_v, residue_p)
            write_kinetic_energy(f"Kinetic_Energy.txt", f"{simulation_time:.6f}", kinetic_energy(particle_velocities))

if __name__ == "__main__":
    main()
