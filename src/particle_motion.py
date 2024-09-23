from concurrent.futures import ThreadPoolExecutor
from fluid_discretization_schemes import *
from parameters import *

def interpolate_velocity_at_particles(u, v, positions, X, Y, element_length_x, element_length_y):
    velocities = np.zeros_like(positions)
    for i, (px, py) in enumerate(positions):
        ix = min(int(px / element_length_x), N_POINTS_X - 2)
        iy = min(int(py / element_length_y), N_POINTS_Y - 2)
        x1, x2 = X[iy, ix], X[iy, ix + 1]
        y1, y2 = Y[iy, ix], Y[iy + 1, ix]
        u_interp = (u[iy, ix] * (x2 - px) * (y2 - py) +
                    u[iy, ix + 1] * (px - x1) * (y2 - py) +
                    u[iy + 1, ix] * (x2 - px) * (py - y1) +
                    u[iy + 1, ix + 1] * (px - x1) * (py - y1)) / ((x2 - x1) * (y2 - y1))
        v_interp = (v[iy, ix] * (x2 - px) * (y2 - py) +
                    v[iy, ix + 1] * (px - x1) * (y2 - py) +
                    v[iy + 1, ix] * (x2 - px) * (py - y1) +
                    v[iy + 1, ix + 1] * (px - x1) * (py - y1)) / ((x2 - x1) * (y2 - y1))
        velocities[i] = [u_interp, v_interp]
    return velocities

def parallel_particle_update(u_next, v_next, particle_positions, X, Y, element_length_x, element_length_y, TIME_STEP_LENGTH):
    with ThreadPoolExecutor() as executor:
        future_velocity = executor.submit(interpolate_velocity_at_particles, u_next, v_next, particle_positions, X, Y, element_length_x, element_length_y)
        particle_velocities = future_velocity.result()
        particle_positions += particle_velocities * TIME_STEP_LENGTH
        # Check if any particle is outside the domain in x or y directions
        mask_x = (particle_positions[:, 0] <= 0) | (particle_positions[:, 0] >= DOMAIN_SIZE_X)
        mask_y = (particle_positions[:, 1] <= 0) | (particle_positions[:, 1] >= DOMAIN_SIZE_Y)

        # Invert velocities for particles that are outside the domain (either in x or y direction)
        mask = mask_x | mask_y  # Logical OR to combine both conditions
        particle_velocities[mask] *= -1
    return particle_positions, particle_velocities
