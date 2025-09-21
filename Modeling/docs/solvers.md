# Solvers Module

The `solvers` module provides numerical methods for integrating ordinary differential equations (ODEs) and delay differential equations (DDEs).

## Ordinary Differential Equation (ODE) Solvers


The ODE solvers module provides a comprehensive suite of single-step, multistep, and adaptive methods for numerically integrating systems of ordinary differential equations (ODEs). The main entry point is the `solve_ode` function, which selects the appropriate solver based on the `method` parameter. Below, we describe the available methods, their mathematical background, and provide usage examples.

### General Usage

```cpp
#include "includes/solvers/solver.hpp"

// Define the ODE right-hand side: dy/dt = f(t, y)
MyFunc f = [](double t, const dVec& y) {
    dVec dydt(y.size());
    for (size_t i = 0; i < y.size(); ++i)
        dydt[i] = -y[i];
    return dydt;
};

SolverParameters params;
params.derivative = f;
params.initial_conditions = {1.0};
params.t0 = 0.0;
params.t1 = 5.0;
params.dt = 0.1;

// Solve using RK4
SolverResults results = solve_ode(params, "RK4");
```

---

## Single-Step Runge-Kutta Methods

### Euler's Method (RK1)

**Function:** `rk1_solver`  
**Order:** 1 (first-order)  
**Description:** The simplest explicit method, using a single derivative evaluation per step. Suitable for rough estimates or stiff problems with very small step sizes.

**Mathematical Formulation:**
$$ y_{n+1} = y_n + h f(t_n, y_n) $$

**Example:**
```cpp
SolverResults results = solve_ode(params, "RK1");
```

---

### Midpoint Method (RK2)

**Function:** `rk2_solver`  
**Order:** 2 (second-order)  
**Description:** Improves upon Euler's method by using a half-step evaluation for better accuracy.

**Mathematical Formulation:**
$$
k_1 = f(t_n, y_n) \\
k_2 = f(t_n + h/2, y_n + (h/2) k_1) \\
y_{n+1} = y_n + h k_2
$$

**Example:**
```cpp
SolverResults results = solve_ode(params, "RK2");
```

---

### Kutta's Third-Order Method (RK3)

**Function:** `rk3_solver`  
**Order:** 3 (third-order)  
**Description:** A standard third-order Runge-Kutta method, offering a balance between accuracy and computational cost.

**Mathematical Formulation:**
$$
k_1 = f(t_n, y_n) \\
k_2 = f(t_n + h/2, y_n + (h/2) k_1) \\
k_3 = f(t_n + h, y_n - h k_1 + 2h k_2) \\
y_{n+1} = y_n + \frac{h}{6}(k_1 + 4k_2 + k_3)
$$

**Example:**
```cpp
SolverResults results = solve_ode(params, "RK3");
```

---

### Classical Fourth-Order Runge-Kutta (RK4)

**Function:** `rk4_solver`  
**Order:** 4 (fourth-order)  
**Description:** The most widely used explicit method, known for its accuracy and stability. Requires four derivative evaluations per step.

**Mathematical Formulation:**
$$
k_1 = f(t_n, y_n) \\
k_2 = f(t_n + h/2, y_n + (h/2) k_1) \\
k_3 = f(t_n + h/2, y_n + (h/2) k_2) \\
k_4 = f(t_n + h, y_n + h k_3) \\
y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)
$$

**Example:**
```cpp
SolverResults results = solve_ode(params, "RK4");
```

---

### RK4 Variants

- **3/8 Rule RK4 (`RK4 -3/8`)**: Uses a different set of coefficients for the final weighted average, trading some accuracy for different stability properties.
- **Ralston's RK4 (`RK4 -r`)**: Optimized for a smaller error bound with a specific set of coefficients.
- **Gill's RK4 (`RK4 -g`)**: A memory-efficient variant that reduces the number of temporary variables.

**Example:**
```cpp
SolverResults results = solve_ode(params, "RK4 -3/8");
SolverResults results = solve_ode(params, "RK4 -r");
SolverResults results = solve_ode(params, "RK4 -g");
```

---

## Multistep Methods

### Adams-Bashforth (AB)

**Function:** `adams_bashforth_solver`  
**Order:** 1–10 (explicit)  
**Description:** A family of explicit multistep methods that use a history of previous derivative evaluations to improve efficiency. Requires a single-step method (e.g., RK4) for bootstrapping the initial steps.

**Mathematical Formulation (general):**
$$
y_{n+1} = y_n + h \sum_{k=0}^{s-1} b_k f_{n-k}
$$
where $b_k$ are the Adams-Bashforth coefficients for order $s$.

**Example:**
```cpp
params.order = 4; // 4th order AB
SolverResults results = solve_ode(params, "AB");
```

---

### Adams-Bashforth-Moulton (ABM)

**Function:** `adams_bashforth_moulton_solver`  
**Order:** 1–10 (predictor-corrector)  
**Description:** A predictor-corrector method combining Adams-Bashforth (predictor) and Adams-Moulton (corrector) steps. Offers improved accuracy and stability, especially for stiff or oscillatory problems. The number of corrector iterations can be set via `params.iterations`.

**Mathematical Formulation:**
Predictor (AB):
$$
y^*_{n+1} = y_n + h \sum_{k=0}^{s-1} b_k f_{n-k}
$$
Corrector (AM):
$$
y_{n+1} = y_n + h \sum_{k=0}^{s} a_k f_{n+1-k}
$$
where $b_k$ and $a_k$ are the Adams-Bashforth and Adams-Moulton coefficients, respectively.

**Example:**
```cpp
params.order = 5; // 5th order ABM
params.iterations = 2; // Number of corrector iterations
SolverResults results = solve_ode(params, "ABM");
```

---

## Adaptive and Embedded Runge-Kutta Methods

These methods provide error estimation and adaptive step size control for efficient and accurate integration.

- **Heun-Euler (RK1(2))**: `solve_ode(params, "Heun-Euler")` — Embedded first/second-order method.
- **Bogacki-Shampine (RK2(3))**: `solve_ode(params, "Bogachi-Shampine")` — Embedded second/third-order method.
- **Fehlberg (RK3(4))**: `solve_ode(params, "RK34")` — Embedded third/fourth-order method.
- **Dormand-Prince (RK4(5))**: `solve_ode(params, "Dormand-Prince")` or `solve_ode(params, "RK45")` — Embedded fourth/fifth-order method (basis of MATLAB's `ode45`).

**Example:**
```cpp
params.error_estimate = true;
params.variable_steps = true;
SolverResults results = solve_ode(params, "Dormand-Prince");
```

---

## Example: Exponential Decay ODE

```cpp
#include "includes/solvers/solver.hpp"

// dy/dt = -y, y(0) = 1
MyFunc f = [](double t, const dVec& y) {
    dVec dydt(y.size());
    for (size_t i = 0; i < y.size(); ++i)
        dydt[i] = -y[i];
    return dydt;
};
SolverParameters params;
params.derivative = f;
params.initial_conditions = {1.0};
params.t0 = 0.0;
params.t1 = 5.0;
params.dt = 0.1;
SolverResults results = rk4_solver(params);
for (size_t i = 0; i < results.time_points.size(); ++i) {
    std::cout << "t = " << results.time_points[i] << ", y = " << results.solution[i][0] << std::endl;
}
```

---

## Example: Multistep ABM Solver

```cpp
#include "includes/solvers/solver.hpp"

// dy/dt = sin(t), y(0) = 0
MyFunc f = [](double t, const dVec& y) {
    dVec dydt(y.size());
    for (size_t i = 0; i < y.size(); ++i)
        dydt[i] = std::sin(t);
    return dydt;
};
SolverParameters params;
params.derivative = f;
params.initial_conditions = {0.0};
params.t0 = 0.0;
params.t1 = 6.28; // ~2*pi
params.dt = 0.1;
params.order = 5; // Use 5th order ABM
params.iterations = 2; // Number of corrector iterations
SolverResults results = adams_bashforth_moulton_solver(params);
for (size_t i = 0; i < 10 && i < results.time_points.size(); ++i) {
    std::cout << "t = " << results.time_points[i] << ", y = " << results.solution[i][0] << std::endl;
}
```

## Delay Differential Equation (DDE) Solvers

The DDE solvers are designed to handle equations where the derivative at a given time depends on the solution at previous times. The primary function for DDEs is `rk1_dde_solver`.

### `rk1_dde_solver` Function

This function implements a first-order Runge-Kutta method (Euler's method) for solving DDEs. It requires a history function for `t < t0` and handles multiple delay times.

#### Interface

```cpp
SolverResults rk1_dde_solver(const SolverParameters& params);
```

Where `SolverParameters` is a structure containing the following DDE-specific fields:

- `delay_derivative`: A function of type `DelayFunc` representing the DDE: `f(t, y, y_delayed)`. `y_delayed` is a vector of values corresponding to `y(t - τ_i)` for each delay time.
- `delay_times`: A `std::vector<double>` containing the delay times `τ_i`.
- `history_function`: A function representing the solution for `t < t0`. This is crucial for initializing the solver.
- `is_dde`: A boolean flag, which should be set to `true` when solving DDEs.

#### Mathematical Background

Delay differential equations have the general form:

$$\frac{dy}{dt} = f(t, y(t), y(t - \tau_1), y(t - \tau_2), ..., y(t - \tau_n))$$

where $\tau_i$ are the delay times. The solver uses linear interpolation to estimate $y(t - \tau)$ from previously computed solution values.

#### Example Usage

Here are examples demonstrating the use of `rk1_dde_solver`:

##### Simple Linear DDE Example (`dde_example.cpp`)

- **Equation**: `dy/dt = -y(t) + 0.5 * y(t - 1.0)`
- **Initial condition**: `y(0) = 1.0`
- **History function**: `y(t) = 1.0` for `t < 0`
- **Delay time**: `τ = 1.0`
- **Time interval**: `[0, 10]`
- **Step size**: `dt = 0.1`

```cpp
// ... (includes and setup)

// Define the DDE function
auto delay_func = [](double t, double y, const std::vector<double>& y_delayed) {
    return -y + 0.5 * y_delayed[0]; // y_delayed[0] corresponds to y(t - 1.0)
};

// Define the history function
auto history_func = [](double t) {
    return 1.0;
};

SolverParameters params;
params.initial_value = 1.0;
params.start_time = 0.0;
params.end_time = 10.0;
params.step_size = 0.1;
params.delay_derivative = delay_func;
params.delay_times = {1.0};
params.history_function = history_func;
params.is_dde = true;

SolverResults results = rk1_dde_solver(params);

// ... (process results)
```

##### Complex Nonlinear DDE Example (`dde_complex_example.cpp`)

- **Equation**: `dy/dt = -2*y(t) + sin(y(t - 0.5)) + 0.3*y(t - 0.5)`
- **Initial condition**: `y(0) = 0.5`
- **History function**: `y(t) = sin(t)` for `t < 0`
- **Delay time**: `τ = 0.5`
- **Time interval**: `[0, 15]`
- **Step size**: `dt = 0.05`

```cpp
// ... (includes and setup)

// Define the DDE function
auto delay_func_complex = [](double t, double y, const std::vector<double>& y_delayed) {
    return -2 * y + std::sin(y_delayed[0]) + 0.3 * y_delayed[0];
};

// Define the history function
auto history_func_complex = [](double t) {
    return std::sin(t);
};

SolverParameters params_complex;
params_complex.initial_value = 0.5;
params_complex.start_time = 0.0;
params_complex.end_time = 15.0;
params_complex.step_size = 0.05;
params_complex.delay_derivative = delay_func_complex;
params_complex.delay_times = {0.5};
params_complex.history_function = history_func_complex;
params_complex.is_dde = true;

SolverResults results_complex = rk1_dde_solver(params_complex);

// ... (process results)
```
