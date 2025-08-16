# Utility Module

The `utility` module provides a collection of general-purpose functions for tasks such as printing data to the console, writing data to files with various formatting options, and helper functions for adaptive solvers and delay differential equations.

## Printing Utilities (`print.hpp`)

These functions are primarily for debugging and quick inspection of data in the console.

### `print_matrix`

Prints a 2D vector (matrix) to `std::cout`, with values separated by spaces and each row on a new line.

```cpp
template <typename t>
void print_matrix(const Vec<t>& matrix);
```

- `matrix`: The matrix (e.g., `Matrix` or `Vec<dVec>`) to print.

**Example:**

```cpp
Matrix my_matrix = {{1.0, 2.0}, {3.0, 4.0}};
print_matrix(my_matrix);
// Output:
// 1.0 2.0 
// 3.0 4.0 
```

### `print_vector`

Prints a 1D vector to `std::cout`, with values separated by spaces.

```cpp
template <typename t>
void print_vector(const Vec<t>& vector);
```

- `vector`: The vector (e.g., `dVec` or `wVec`) to print.

**Example:**

```cpp
dVec my_vector = {1.0, 2.5, 3.7};
print_vector(my_vector);
// Output:
// 1.0 2.5 3.7 
```

### `print_progress_bar_eta`

Prints a progress bar to the console, including an estimated time of arrival (ETA).

```cpp
template <typename T>
void print_progress_bar_eta(size_t current, size_t total, T eta, size_t bar_width = 50);
```

- `current`: The current progress count.
- `total`: The total count for completion.
- `eta`: The estimated time remaining (can be any type that can be streamed to `std::cout`).
- `bar_width`: The width of the progress bar in characters (default: 50).

**Example:**

```cpp
// In a loop:
// print_progress_bar_eta(i, total_iterations, estimated_time_remaining);
```

### `print_progress_bar`

Prints a progress bar to the console without an ETA.

```cpp
void print_progress_bar(size_t current, size_t total, size_t bar_width = 50);
```

- `current`: The current progress count.
- `total`: The total count for completion.
- `bar_width`: The width of the progress bar in characters (default: 50).

**Example:**

```cpp
// In a loop:
// print_progress_bar(i, total_iterations);
```

## File Writing Utilities (`write.hpp`)

These functions provide flexible options for writing numerical data to files, supporting various formats and configurations.

### `write_matrix_to_file`

Writes a 2D vector (matrix) to a specified file with extensive formatting options.

```cpp
template <typename t>
void write_matrix_to_file(
    const Vec<t>&      matrix,
    const std::string& directory,
    const std::string& filename,
    const std::string& filetype                = "txt",
    const std::string& separator               = " ",
    const std::string& comment                 = "",
    const std::string& precision               = "16",
    const std::string& format                  = "fixed",
    const std::string& alignment               = "",
    const std::string& header                  = "",
    const std::string& footer                  = "",
    bool               append                  = false,
    size_t             col_width               = 0,
    bool               skip_trailing_separator = false,
    bool               binary                  = false
);
```

- `matrix`: Matrix to write.
- `directory`: Output directory path.
- `filename`: Output filename (without extension).
- `filetype`: File extension (default: "txt").
- `separator`: Separator between values (default: " ").
- `comment`: Optional comment string (added after header).
- `precision`: String for floating-point precision (number of digits after decimal, default: "16").
- `format`: String for floating-point format ("fixed" or "scientific", default: "fixed").
- `alignment`: String for column alignment ("left", "right", or "", default: "").
- `header`: Optional header line.
- `footer`: Optional footer line.
- `append`: If `true`, append to file instead of overwriting (default: `false`).
- `col_width`: Column width for alignment (default: `0` = no width).
- `skip_trailing_separator`: If `true`, do not write separator after last value in row (default: `false`).
- `binary`: If `true`, write as raw binary (default: `false`).

**Example:**

```cpp
Matrix data = {{1.234, 5.678}, {9.012, 3.456}};
write_matrix_to_file(data, "./output", "results", "csv", ",", "Simulation Data", "6", "fixed", "", "Time,Value1,Value2");
```

### `write_vector_to_file`

Writes a 1D vector to a specified file with extensive formatting options.

```cpp
template <typename t>
void write_vector_to_file(
    const Vec<t>&      vector,
    const std::string& directory,
    const std::string& filename,
    const std::string& filetype                = "txt",
    const std::string& separator               = " ",
    const std::string& comment                 = "",
    const std::string& precision               = "16",
    const std::string& format                  = "fixed",
    const std::string& alignment               = "",
    const std::string& header                  = "",
    const std::string& footer                  = "",
    bool               append                  = false,
    size_t             col_width               = 0,
    bool               skip_trailing_separator = false,
    bool               binary                  = false
);
```

- Parameters are analogous to `write_matrix_to_file` but for a 1D vector.

**Example:**

```cpp
dVec data_vec = {1.1, 2.2, 3.3, 4.4};
write_vector_to_file(data_vec, "./output", "single_run", "dat", "\t", "Single Run Results");
```

## Adaptive Solvers Utility (`adaptive_solvers_utility.hpp`)

These functions provide helper logic for adaptive step-size numerical solvers.

### `estimate_error`

Estimates the error between a current solution `y` and a reference solution `y_ref` based on specified error estimation methods (weighted or L2 norm).

```cpp
double estimate_error(const dVec& y, const dVec& y_ref, const SolverParameters& Params);
```

- `y`: The current solution vector.
- `y_ref`: The reference solution vector.
- `Params`: `SolverParameters` object containing `weighted_error`, `norm_error`, and `absolute_tol`.

**Returns:** The estimated error.

### `adjust_step_size`

Adjusts the step size `current_dt` based on the estimated `error` and solver parameters. This is typically used in adaptive step-size algorithms to maintain a desired level of accuracy.

```cpp
double adjust_step_size(double current_dt, double error, const SolverParameters& Params, size_t p);
```

- `current_dt`: The current step size.
- `error`: The estimated error from `estimate_error`.
- `Params`: `SolverParameters` object containing `local_tol`, `local_tol_error_ratio`, `min_dt`, `max_dt`, `decrease_factor`, and `increase_factor`.
- `p`: The order of the numerical method.

**Returns:** The adjusted step size.

## DDE Utility (`dde_utility.hpp`)

These functions are specifically designed to assist Delay Differential Equation (DDE) solvers, particularly for handling historical data and multiple delays.

### `interpolate_delay`

Interpolates the solution at a specific past time `t_delayed` using linear interpolation from a set of `time_points` and corresponding `solution` states. If `t_delayed` is before the initial time `t0`, it uses a provided `history_function`.

```cpp
inline dVec interpolate_delay(
    double        t_delayed,             // Time to interpolate at (t - tau)
    const dVec&   time_points,           // Available time points
    const Matrix& solution,              // Solution at those time points
    double        t0,                    // Initial time
    const Func&   history_function       // Pre-history function values
);
```

- `t_delayed`: The time at which to interpolate the solution (e.g., `t - tau`).
- `time_points`: A vector of time points for which solution data is available.
- `solution`: A matrix where each row corresponds to the solution state at a `time_point`.
- `t0`: The initial time of the simulation.
- `history_function`: A function that provides the solution for times `t < t0`.

**Returns:** A `dVec` representing the interpolated solution state at `t_delayed`.

**Example:**

```cpp
// Assuming time_points, solution, t0, and history_func are defined
dVec y_at_delay = interpolate_delay(current_time - delay_tau, time_points, solution_history, initial_t, history_func);
```

### `interpolate_multiple_delays`

Interpolates the solution for multiple delay times. It calls `interpolate_delay` for each specified delay and returns a matrix where each row corresponds to the interpolated state for a given delay time.

```cpp
inline Matrix interpolate_multiple_delays(
    const dVec&   delay_times,           // Vector of delay times τ_i
    double        t,                     // Current time
    const dVec&   time_points,           // Available time points
    const Matrix& solution,              // Solution at those time points
    double        t0,                    // Initial time
    const Func&   history_function       // Pre-history function values
);
```

- `delay_times`: A vector of delay times `τ_i`.
- `t`: The current simulation time.
- Other parameters are the same as `interpolate_delay`.

**Returns:** A `Matrix` where each row `i` contains the interpolated state `y(t - τ_i)`.

**Example:**

```cpp
// Assuming delay_times, current_t, time_points, solution_history, initial_t, and history_func are defined
Matrix delayed_states = interpolate_multiple_delays(delay_times_vec, current_t, time_points, solution_history, initial_t, history_func);
```
