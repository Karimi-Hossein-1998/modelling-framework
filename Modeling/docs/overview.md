# Mathematical Modeling Toolkit: A C++ Library for Scientific Computing

Mathematical Modeling Toolkit is a comprehensive C++ library designed for scientific computing, offering a collection of tools for solving ordinary differential equations (ODEs), delay differential equations (DDEs), network topology generation, and various linear algebra operations. The library is structured to be modular, allowing users to easily integrate specific components into their projects.

## Core Components

The library is organized into several key modules:

- **Solvers**: Provides numerical methods for solving ODEs and DDEs.
- **Initializers**: Contains functions for setting up initial conditions for simulations.
- **Network**: Tools for generating and manipulating network topologies.
- **Models**: Implementations of various mathematical models, such as the Kuramoto model.
- **Utility**: General utility functions, including printing and file writing.
- **Linear Algebra**: Basic linear algebra operations, including matrix decompositions and linear system solvers.
- **Interpolators**: Functions for interpolation, such as Lagrange and Newton interpolation.
- **Typedefs**: Custom type definitions for complex numbers and other common data structures.

## Getting Started

To use Mathematical Modeling Toolkit, include the necessary header files in your C++ project. The main header file, `mm.hpp`, provides a convenient way to include all core components:

```cpp
#include "mm.hpp"
```

For specific functionalities, you can include individual headers, for example (note: these headers may change as the library evolves and you need to make sure the paths are all relative):

```cpp
#include "solvers/dde_solvers.hpp"
#include "models/kuramoto.hpp"
```

## Examples

Refer to the [`Examples/`](../Examples/) directory for various usage examples demonstrating how to utilize different parts of the library.

## Documentation Structure

This documentation is organized into several files to provide a clear and detailed explanation of each module:

- `overview.md`: General introduction and library structure. [Overview](./overview.md)
- `solvers.md`: Detailed documentation for ODE and DDE solvers. [Solvers](./solvers.md)
- `initializers.md`: Documentation for initial condition setup. [Initializers](./initializers.md)
- `network.md`: Documentation for network topology tools. [Network Topology](./network.md)
- `models.md`: Documentation for implemented mathematical models. [Models](./models.md)
- `utility.md`: Documentation for general utility functions. [Utility](./utility.md)
- `linalg.md`: Documentation for linear algebra components. [Linear Algebra](./linalg.md)
- `interpolators.md`: Documentation for interpolation methods. [Interpolators](./interpolators.md)
- `typedefs.md`: Documentation for custom type definitions. [Type Definitions](./typedefs.md)
