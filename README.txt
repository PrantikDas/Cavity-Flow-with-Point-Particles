Solves the incompressible Navier Stokes equations in a lid-driven cavity
scenario using Finite Differences, explicit timestepping and Chorin's Projection.

Momentum:           ∂u/∂t + (u ⋅ ∇) u = − 1/ρ ∇p + ν ∇²u + f

Incompressibility:  ∇ ⋅ u = 0


u:  Velocity (2d vector)
p:  Pressure
f:  Forcing (here =0)
ν:  Kinematic Viscosity
ρ:  Density
t:  Time
∇:  Nabla operator (defining nonlinear convection, gradient and divergence)
∇²: Laplace Operator

----

Lid-Driven Cavity Scenario:


                            ------>>>>> u_top

          1 +-------------------------------------------------+
            |                                                 |
            |             *                      *            |
            |          *           *    *    *                |
        0.8 |                                                 |
            |                                 *               |
            |     *       *                                   |
            |                      *     *                    |
        0.6 |                                            *    |
u = 0       |      *                             *            |   u = 0
v = 0       |                             *                   |   v = 0
            |                     *                           |
            |           *                *         *          |
        0.4 |                                                 |
            |                                                 |
            |      *            *             *               |
            |           *                             *       |
        0.2 |                       *           *             |
            |                               *                 |
            |  *          *      *                 *       *  |
            |                            *                    |
          0 +-------------------------------------------------+
            0        0.2       0.4       0.6       0.8        1

                                    u = 0
                                    v = 0

* Velocity and pressure have zero initial condition.
* Homogeneous Dirichlet Boundary Conditions everywhere except for horizontal
  velocity at top. It is driven by an external flow.

-----

Solution strategy:   (Projection Method: Chorin's Splitting)

1. Solve Momentum equation without pressure gradient for tentative velocity
   (with given Boundary Conditions)

    ∂u/∂t + (u ⋅ ∇) u = ν ∇²u

2. Solve pressure poisson equation for pressure at next point in time
   (with homogeneous Neumann Boundary Conditions everywhere except for
   the top, where it is homogeneous Dirichlet)

    ∇²p = ρ/Δt ∇ ⋅ u           

3. Correct the velocities (and again enforce the Velocity Boundary Conditions)

    u ← u − Δt/ρ ∇ p

-----

    Expected Outcome: After some time a swirling motion will take place

          1 +-------------------------------------------------+
            |                                                 |
            |                                                 |
            |                                                 |
        0.8 |                                                 |
            |                      *-->*                      |
            |                ******     ******                |
            |              **                 **              |
        0.6 |             *                     *             |
            |             *                      *            |
            |            *                        *           |
            |            *                       *            |
            |             *                     *             |
        0.4 |             *                     *             |
            |              **                 **              |
            |                ******     ******                |
            |                      *<--*                      |
        0.2 |                                                 |
            |                                                 |
            |                                                 |
            |                                                 |
          0 +-------------------------------------------------+
            0        0.2       0.4       0.6       0.8        1

------

Strategy in index notation

u = [u, v]
x = [x, y]

1. Solve tentative velocity + velocity BC

    ∂u/∂t + u ∂u/∂x + v ∂u/∂y = ν ∂²u/∂x² + ν ∂²u/∂y²

    ∂v/∂t + u ∂v/∂x + v ∂v/∂y = ν ∂²v/∂x² + ν ∂²v/∂y²

2. Solve pressure poisson + pressure BC

    ∂²p/∂x² + ∂²p/∂y² = ρ/Δt (∂u/∂x + ∂v/∂y)

3. Correct velocity + velocity BC

    u ← u − Δt/ρ ∂p/∂x

    v ← v − Δt/ρ ∂p/∂y

------

IMPORTANT: Take care to select a timestep that ensures stability
"""

if parallel run it as : mpiexec -np 4 make run
if serial run it as : make run 
if you want to clean the system use : make clean
