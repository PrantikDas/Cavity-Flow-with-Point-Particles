import numpy as np
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

def central_difference_x(f, element_length):
    diff = np.zeros_like(f)
    diff[1:-1, 1:-1] = (f[1:-1, 2:] - f[1:-1, 0:-2]) / (2 * element_length)
    return diff

def central_difference_y(f, element_length):
    diff = np.zeros_like(f)
    diff[1:-1, 1:-1] = (f[2:, 1:-1] - f[0:-2, 1:-1]) / (2 * element_length)
    return diff

def laplace(f, element_length):
    diff = np.zeros_like(f)
    diff[1:-1, 1:-1] = (f[1:-1, 0:-2] + f[0:-2, 1:-1] - 4 * f[1:-1, 1:-1] + f[1:-1, 2:] + f[2:, 1:-1]) / (element_length**2)
    return diff
