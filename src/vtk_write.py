import os
from parameters import *
from energy_budget import *

def write_vtk(filename, X, Y, u, v, p, cell_volume_fractions):
    filepath = os.path.join(output_folder_fluid, filename)
    with open(filepath, "w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("Lid Driven Cavity Simulation Results\n")
        f.write("ASCII\n")
        f.write("DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {N_POINTS_X} {N_POINTS_Y} 1\n")
        f.write("ORIGIN 0 0 0\n")
        f.write(f"SPACING {DOMAIN_SIZE_X / (N_POINTS_X - 1)} {DOMAIN_SIZE_Y / (N_POINTS_Y - 1)} {1.0}\n")
        f.write(f"POINT_DATA {N_POINTS_X * N_POINTS_Y}\n")

        # Write velocity data
        f.write("VECTORS velocity float\n")
        for j in range(N_POINTS_Y):
            for i in range(N_POINTS_X):
                f.write(f"{u[j, i]} {v[j, i]} 0.0\n")

        # Write pressure data
        f.write("SCALARS pressure float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for j in range(N_POINTS_Y):
            for i in range(N_POINTS_X):
                f.write(f"{p[j, i]}\n")

        # Write cell volume fraction data
        f.write("SCALARS volume_fraction float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for j in range(N_POINTS_Y):
            for i in range(N_POINTS_X):
                # Ensure we are within bounds of cell_volume_fractions
                if j < cell_volume_fractions.shape[0] and i < cell_volume_fractions.shape[1]:
                    cell_volume_fraction = cell_volume_fractions[j, i]
                else:
                    cell_volume_fraction = 0.0  # Default to 0.0 if out of bounds
                f.write(f"{cell_volume_fraction}\n")

def write_particles_liggghts(filename, timestep, particle_positions, particle_velocities):
    filepath = os.path.join(output_folder_particle, filename)
    with open(filepath, "w") as f:
        f.write("ITEM: TIMESTEP\n")
        f.write(f"{timestep}\n")
        num_atoms = len(particle_positions)
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(f"{num_atoms}\n")
        f.write("ITEM: BOX BOUNDS ff ff ff\n")
        f.write("0 1\n")
        f.write("0 1\n")
        f.write("0 0\n")
        f.write("ITEM: ATOMS id type x y z vx vy vz\n")
        for i, (pos, vel) in enumerate(zip(particle_positions, particle_velocities)):
            f.write(f"{i + 1} 1 {pos[0]} {pos[1]} 0.0 {vel[0]} {vel[1]} 0.0\n")
            
def write_residuals(filename, simulation_time, residue_u, residue_v, residue_p):
    filepath = os.path.join(output_folder_residual, filename)
    with open(filepath, "a") as f:
        f.write(f"{simulation_time} {residue_u} {residue_v} {residue_p}\n")
