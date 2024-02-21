# Fluid Simulation Framework 
Practical Coursework at RWTH Aachen University

## Overview

The Fluid Simulation Framework is a comprehensive C++ library designed for the simulation of fluid dynamics using the Smoothed Particle Hydrodynamics (SPH) method. 
This project aims to provide an efficient and scalable solution for simulating both compressible and incompressible fluids, incorporating advanced physics phenomena such as viscosity, incompressibility, and solid-fluid interactions. 
The framework facilitates the numerical solution of the Navier-Stokes equations through spatial discretization and time integration schemes, allowing for the accurate modeling of fluid motion and behavior.
For more detailed information see [here](https://animation.rwth-aachen.de/course/48/).

## Features

- **SPH Implementation**: Core SPH algorithms for fluid dynamics simulation, offering detailed control over particle interactions and fluid properties.
- **Numerical Integration**: Robust time integration schemes to ensure stability and accuracy over the simulation.
- **Physics Extensions**: Modules for simulating complex physical phenomena, including viscosity, surface tension, and interactions between fluids and solids.
- **Visualization Support**: Tools for exporting simulation data for visualization, including particle systems and surface reconstruction with triangular meshes.
- **Modular Architecture**: Designed for extensibility and collaboration, allowing for easy integration of new features and optimizations.

## Getting Started

### Prerequisites

- C/C++ compiler with C++17 support
- Basic knowledge of numerics, algorithms, and data structures
- Optional: External tools for visualization (e.g., Blender, MeshLab, ParaView)