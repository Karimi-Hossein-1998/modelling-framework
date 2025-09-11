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

For specific functionalities, you can include individual headers, for example:

```cpp
#include "solvers/dde_solvers.hpp"
#include "models/kuramoto.hpp"
```

## Examples

Refer to the `Examples/` directory for various usage examples demonstrating how to utilize different parts of the library.

## Documentation Structure

This documentation is organized into several files to provide a clear and detailed explanation of each module:

- `overview.md`: General introduction and library structure.
- `solvers.md`: Detailed documentation for ODE and DDE solvers.
- `initializers.md`: Documentation for initial condition setup.
- `network.md`: Documentation for network topology tools.
- `models.md`: Documentation for implemented mathematical models.
- `utility.md`: Documentation for general utility functions.
- `linalg.md`: Documentation for linear algebra components.
- `interpolators.md`: Documentation for interpolation methods.
- `typedefs.md`: Documentation for custom type definitions.

## Development Roadmap

For a detailed plan of future additions and features, please refer to the [Development Roadmap](Modelling/docs/roadmap.md).

# Contributing to the Mathematical Modeling Toolkit

Thank you for your interest in contributing. This project is built on a foundation of **honesty, practicality, theoretical integrity, correctness (rigor),** and **clean code**. Our goal is to create a robust and reliable toolkit for scientific computing, and all contributions are valued based on their technical merit and adherence to these principles.

---

### How Can I Contribute?

There are many ways to contribute to this project, and not all of them require writing code. Here are some options:

* **Report a Bug**: If you find an issue, please report it on the [Issues](https://github.com/Karimi-Hossein-1998/modelling-framework/issues) page. Provide a clear, concise description, including steps to reproduce the bug, so that we can verify and address the problem efficiently.
* **Suggest an Enhancement**: Have an idea for a new feature? We'd be glad to discuss it. Open an issue to propose your idea, explaining its practicality and theoretical basis.
* **Improve Documentation**: Clarity and correctness are paramount. If you find documentation that is unclear, inaccurate, or could be improved, please let us know.
* **Contribute Code**: If you're ready to contribute code, we welcome your work on existing issues or new features.

---

### Submitting a Pull Request

If you are contributing code, please ensure your work aligns with our core principles. Follow these steps:

1.  **Fork the Repository**: Create a copy of the project in your own GitHub account.
2.  **Create a New Branch**: Make a new branch for your feature or bug fix: `git checkout -b feature/your-feature-name`.
3.  **Code with Rigor**: Write code that is clean, efficient, and well-commented. All new code must be theoretically sound and free of sloppy implementation.
4.  **Write Tests**: All new functionality must be accompanied by comprehensive tests to ensure its correctness and integrity.
5.  **Commit Your Changes**: Write clear and descriptive commit messages. Explain **what** you changed and, more importantly, **why** you changed it.
6.  **Push and Submit**: Push your branch to your forked repository and submit a pull request. In the pull request description, reference the issue it addresses and provide a brief summary of your changes.
