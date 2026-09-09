# Siconos

[![Pipeline Status](https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/badges/main/pipeline.svg)](https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos/-/commits/main)
[![Version](https://img.shields.io/github/release/siconos/siconos.svg)](https://github.com/siconos/siconos/releases/latest)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://github.com/siconos/siconos/blob/main/COPYING)

**A software package for the modeling and simulation of nonsmooth dynamical systems in C++ and Python.**

> 📚 **Examples and Tutorials**: Most examples and tutorials are in the separate **[siconos-tutorials](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials)** repository. Please ensure you use **consistent versions** of siconos and siconos-tutorials (e.g., both `main` branch or both same release version).

> 📖 **Documentation**: check **[getting started, install, user or developer guides, complete (C++ and python) API reference ...](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos)**, getting started, user or developer guides, complete (C++ and python) API reference ...

[Homepage](http://siconos.org) |
[Documentation](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos) |
[Tutorials](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials) |
[GitLab](https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos)

---

## What is Siconos?

Siconos is an open-source scientific software package designed for modeling and simulating **nonsmooth dynamical systems** - systems where the dynamics are not continuous due to phenomena like impacts, friction, switches, or constraints.

### Key Application Areas

| Domain | Examples |
|--------|----------|
| **Mechanical Systems** | Rigid body dynamics with contact and friction, multibody systems, granular materials |
| **Electrical Circuits** | Power converters, rectifiers, PLLs, ADCs with ideal components |
| **Control Systems** | Sliding mode control, hybrid systems |
| **Biology** | Gene regulatory networks |
| **Optimization** | Complementarity systems, variational inequalities |

---

## Quick Start

### Option 1: Try Online (No Installation)

Run Siconos directly in your browser via Binder:

[![Binder](https://plmbinder.math.cnrs.fr/binder/badge_logo.svg)](https://plmbinder.math.cnrs.fr/binder/v2/git/https%3A%2F%2Fgricad-gitlab.univ-grenoble-alpes.fr%2Fnonsmooth%2Fsiconos-tutorials.git/HEAD)

**Best for:** Testing the software, tutorials, lightweight demos

### Option 2: Docker (Recommended for Development)

Prerequisite: [Docker](https://docs.docker.com/get-docker/)

```bash
# Jupyter Lab with Siconos and tutorials
docker run --rm -p 8888:8888 -ti \
  gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/siconoslab-main:latest

# Or command-line only
docker run --rm --entrypoint /bin/bash -ti \
  gricad-registry.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/siconoslab-main:latest
```

### Option 3: Build from Source

```bash
# Clone the repository
git clone https://github.com/siconos/siconos.git
cd siconos

# Configure and build
cmake -S . -B build -DUSER_OPTIONS_FILE=config_samples/default.cmake
cmake --build build -j$(nproc)

# Run tests (optional)
ctest --test-dir build

# Install
cmake --install build
```

📖 [Detailed installation guide](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/install_guide/index.html)

---

## Main Components

Each component can be used from C/C++ or Python.

### [siconos/numerics](numerics) (C)

Low-level numerical algorithms for solving optimization problems:
- **Complementarity problems**: LCP, MLCP, NCP
- **Friction-contact problems**: 2D/3D with Coulomb friction, rolling friction
- **Variational inequalities**: AVI, VI
- **Second-order cone programming (SOCP)**
- **Relay problems**

📖 [Numerics documentation](numerics/README.md)

### [siconos/kernel](kernel) (C++)

Core library for modeling and simulation:
- Dynamical systems: Lagrangian, Newton-Euler, first-order formulations
- Numerical integration: Event-driven and time-stepping schemes
- Nonsmooth laws: Impact, friction, complementarity, relay

### [siconos/mechanics](mechanics) (C++)

Mechanical simulation with contact detection:
- Collision detection (native and [Bullet3](https://github.com/bulletphysics/bullet3))
- CAD integration ([Open CASCADE](https://github.com/tpaviot/oce))

### [siconos/control](control) (C++)

Control systems with sliding mode control and implicit discretization.

### [siconos/io](io) (C++)

Serialization and visualization:
- HDF5 import/export
- VTK visualization

---

## Examples by Domain

The [siconos-tutorials](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials) repository contains comprehensive examples for each application area:

| Domain | Example Topics | Link |
|--------|---------------|------|
| **Mechanics** | Bouncing ball, slider-crank, gear transmission, robotics, granular materials | [Mechanics Tutorials](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials/tutorials/mechanics.html) |
| **Electronics** | Power converters, rectifiers, PLLs, circuit simulation | [Electronics Tutorials](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials/tutorials/electronics.html) |
| **Control** | Sliding mode control, state observers, hybrid systems | [Control Tutorials](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials/tutorials/control.html) |
| **Numerics** | LCP solvers, friction contact, optimization problems | [Numerics Tutorials](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials/tutorials/numerics.html) |

📁 **Browse all examples**: [siconos-tutorials repository](https://gricad-gitlab.univ-grenoble-alpes.fr/nonsmooth/siconos-tutorials/-/tree/main/examples)

---

## Example: Bouncing Ball

```python
import siconos.modeling as sm
import siconos.simulation as ss
import numpy as np

# Parameters
t0, T, h = 0, 10, 0.005    # Time settings
r, g, m = 0.1, 9.81, 1      # Ball properties
e = 0.9                     # Restitution coefficient

# Create dynamical system (ball)
position = [1, 0, 0]
velocity = [0, 0, 0]
mass = np.eye(3) * m
mass[2, 2] = 2./5 * r * r

ball = sm.LagrangianLinearTIDS(position, velocity, mass)
ball.set_constant_fext(np.array([-m * g, 0, 0]))

# Interaction with floor
nslaw = sm.NewtonImpactNSL(e)
relation = sm.LagrangianLinearTIR(np.array([[1, 0, 0]]))
interaction = sm.Interaction(nslaw, relation)

# Model
model = sm.NonSmoothDynamicalSystem(t0, T)
model.insert_dynamical_system(ball)
model.link(interaction, ball)

# Simulation
osi = ss.MoreauJeanOSI(0.5)
td = ss.TimeDiscretisation(t0, h)
osnspb = ss.LCP()
sim = ss.TimeStepping(model, td, osi, osnspb)

# Run simulation
while sim.has_next_event():
    sim.compute_one_step()
    sim.next_step()
```

---

## Documentation & Resources

| Resource | Link |
|----------|------|
| **User Guide** | [https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/users_guide/index.html](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/users_guide/index.html) |
| **Tutorials** | [https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos-tutorials) |
| **API Reference** | [https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/reference/index.html](https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/reference/index.html) |
| **Bug Reports** | [GitHub Issues](https://github.com/siconos/siconos/issues) |

---

## Citation

If you use Siconos in your research, please cite:

```bibtex
@software{siconos,
  title = {Siconos: A Software Package for Modeling and Simulation of Nonsmooth Dynamical Systems},
  author = {{Siconos Team}},
  url = {http://siconos.org},
  year = {2024}
}
```

---

## License

Siconos is distributed under the [Apache License, Version 2.0](COPYING).

---

**Maintained by:** The Siconos Development Team  
**Contact:** [GitHub Issues](https://github.com/siconos/siconos/issues)
