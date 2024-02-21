# Notes per Exercise
- Ex. 1
    - What was done to check correctness of Kernel (also for precomputation)
    - Analytical gradient vs finite difference approximation
    - ffmpeg -r FRAME_RATE <=> Use this!!!
    - CompactNSearch vs Bruteforce (do we need this??)
    - We had to try different smoothing lengths and compare error (e.g. against distance to (0,0,0))
- Ex. 2
    - Density estimation correctness + Accuracy of sampling (e.g. by density estimation after sampling)
    - Boundary Sampling (Square vs Hexagonal) + Obtuse angle + Triangles with circumscribed sphere smaller than half particle distance sampling
    - Without Velocity Smoothing vs With Velocity Smoothing
    - High velocity particle handling (speed cap, drag forces)
    - Tests:
        - Everything disabled => should be a constant box
        - Pressure Stiffness on => Slowly expanding box
        - Fluid viscosity on => ?
        - Boundary Viscosity on => ?
        - No viscosity but particle smoothing
- Ex. 3
    - Reconstruct sphere and torus
    - Optimizations (avoid looping through all edges and cells, only compute stuff that contributes to final mesh)
    - Decoupled reconstruction and simulation
    - Different parameters for reconstruction (c, cell edge length)
- Ex. 4
    - Comparison of WCSPH and PBF (number of particles, timestep size, constraint iterations and pressure stiffness) - Double Dam Break
    - Showcase main differences
- Ex. 5
    - Water Droplet for Cohesion
    - Boundary Wetting for Adhesion
    - Look in paper for other scenarios showcasing cohesion and adhesion

# Features of Simulation to Show/Verify
- Pressure Forces
- Fluid-Viscosity
- Boundary-Viscosity
- Cohesion
- Adhesion
- Squared and Hexagonal Boundary Sampling + Problems with obtuse and very small triangles
- High Velocity Particle Handling (Speedcap, Drag forces)

# Scenarios to Show
- Squared and Hexagonal Boundary Sampling Side by Side
- Pressure Forces
    - Dam Break
- Fluid-Viscosity
    - Dam Break
    - Stair Drop
- Boundary-Viscosity
    - Diagonal Plate
- Cohesion
    - Water Droplet
    - Water fountain
- Adhesion
    - Boundary Wetting (box with hole -> drops onto sphere or monkey or other cool object)
    - Water fountain
- Fancy
    - Waterfall with rocks sticking out (for adhesion)

# Problems to Show
- Boundary Sampling -> Obtuse + Very Small Triangles
- Parameters are dependent on scale (use same parameters for differently scaled simulation -> will be different)
