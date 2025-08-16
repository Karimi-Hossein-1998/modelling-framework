# Typedefs Module

The `typedefs` module defines various type aliases and structures used throughout the MyLib library for clarity, convenience, and to manage complex number operations, solver parameters, and sparse data structures.

## Core Type Aliases (`header.hpp`)

These aliases simplify the declaration of common data structures.

- `Vec<T>`: A template alias for `std::vector<T>`. Provides a shorthand for vector types.
- `dVec`: Alias for `Vec<double>`, representing a 1D vector of double-precision floating-point numbers.
- `wVec`: Alias for `Vec<size_t>`, representing a 1D vector of `size_t` integers, typically used for indices or counts.
- `triple`: Alias for `std::array<double, 3>`, useful for 3D vectors or points.
- `Matrix`: Alias for `Vec<dVec>`, representing a 2D matrix (vector of vectors of doubles).
- `func`: Alias for `std::function<double(double)>`, representing a scalar-valued function of a single double argument.
- `Func`: Alias for `std::function<dVec(double)>`, representing a vector-valued function of a single double argument.
- `MyFunc`: Alias for `std::function<dVec(double, const dVec&)>`, representing an Ordinary Differential Equation (ODE) function `f(t, y)`.
- `DelayFunc`: Alias for `std::function<dVec(double, const dVec&, const dVec&)>`, representing a Delay Differential Equation (DDE) function `f(t, y, y_delayed)` for a single delay.
- `DelayFuncMulti`: Alias for `std::function<dVec(double, const dVec&, const Matrix&)>`, representing a DDE function `f(t, y, y_delayed_states)` for multiple delays.
- `interPolater`: Alias for `std::function<func(const dVec&, const dVec&)>`, representing a function that creates a scalar interpolator.
- `InterPolater`: Alias for `std::function<Func(const dVec&, const Matrix&)>`, representing a function that creates a vector-valued interpolator.

## Constants (`header.hpp`)

- `PI`: `constexpr double` representing the mathematical constant $\pi$.
- `max_order`: `constexpr int` for maximum order (e.g., for Adams-Bashforth methods).
- `myexp`: `constexpr double` for a custom exponential constant.

## Solver Structures (`header.hpp`)

### `SolverParameters`

A comprehensive structure for configuring numerical solvers for ODEs and DDEs. This structure provides fine-grained control over the solver's behavior, accuracy, and performance.

#### Error Estimation and Step Size Control

The adaptive solvers use sophisticated error estimation and step size control mechanisms:

1. **Error Estimation Methods**:
   - Local truncation error is estimated using embedded Runge-Kutta pairs
   - For vector systems, component-wise weighted errors or norm-based errors are supported
   - Mixed absolute/relative error control: $err_i = \frac{|y_i - \hat{y}_i|}{\max(|y_i|, \text{absolute\_tol})}$

2. **Step Size Control Algorithm**:
   ```
   1. Compute solution y and error estimate E
   2. Calculate error ratio r = ‖E‖/local_tol
   3. If r > 1 (error too large):
      - New step = current_step * min(decrease_factor * (1/r)^(1/p), 1)
      - Reject step and try again
   4. If r < 1 (error acceptable):
      - Accept step
      - If r < 0.1 (error very small):
        New step = min(current_step * increase_factor, max_dt)
   ```
   where p is the method order

3. **Stability Considerations**:
   - Step size changes limited by decrease_factor and increase_factor
   - Maximum step size prevents missing rapid changes
   - Minimum step size prevents excessive computation in stiff regions
   - Safety factor (local_tol_error_ratio) provides margin against oscillations

4. **Error Control Parameters**:
   - local_tol: Primary accuracy control
   - absolute_tol: Prevents division by zero in relative error
   - weighted_error: Accounts for different scales in system components
   - norm_error: Alternative for oscillatory systems

```cpp
struct SolverParameters {
    // Core Parameters
    MyFunc derivative;                           // The ODE right-hand side function f(t,y) in dy/dt = f(t,y)
    dVec   initial_conditions;                   // Initial state vector y(t0)
    double t0                    = 0.0;         // Initial time
    double t1                    = 10.0;        // Final time
    double dt                    = 1.0 / 16.0;  // Initial step size (can be adapted if variable_steps=true)

    // Multi-step Method Parameters
    int    order                 = 4;           // Order of accuracy for Adams-Bashforth/Moulton methods (1-10)
    int    iterations            = 2;           // Number of corrector iterations in Adams-Moulton method
                                               // More iterations increase accuracy but cost more computation

    // Adaptive Step Size Control
    double min_dt                = 1.0 / 128.0; // Minimum allowed step size for adaptive methods
                                               // Prevents excessive step size reduction in stiff regions
    double max_dt                = 0.5;         // Maximum allowed step size
                                               // Ensures sufficient temporal resolution
    double decrease_factor       = 0.9;         // Factor to decrease step size when error is too large
                                               // Should be < 1 for stability
    double increase_factor       = 2.0;         // Factor to increase step size when error is small
                                               // Usually 1.5-2.0 for efficiency while maintaining stability

    // Error Control Parameters
    double local_tol             = 1.0e-7;      // Target local error tolerance per step
                                               // Primary control for accuracy vs. speed tradeoff
    double absolute_tol          = 1.0e-10;     // Absolute error tolerance
                                               // Prevents division by zero in relative error calculations
                                               // Important for components that pass through zero
    double local_tol_error_ratio = 0.1;         // Safety factor for step size adjustment
                                               // Smaller values are more conservative
    size_t max_trial             = 25;          // Maximum attempts to achieve desired accuracy per step
                                               // Prevents infinite loops in difficult regions

    // Solver Behavior Flags
    bool   error_estimate        = false;       // Enable error estimation using embedded RK pairs
                                               // Required for adaptive step size control
    bool   variable_steps        = false;       // Enable adaptive step size control
                                               // Automatically adjusts step size based on error estimates
    bool   attempts_history      = false;       // Store history of failed step attempts
                                               // Useful for debugging and algorithm analysis
    bool   weighted_error        = true;        // Use weighted error formula for vector-valued ODEs
                                               // Accounts for different scales in system components
    bool   norm_error            = false;       // Use solution norm instead of component-wise error
                                               // Better for some types of oscillatory problems

    // DDE-specific Parameters
    
    // Core DDE Function Types
    DelayFunc delay_derivative   = nullptr;     // Single-delay DDE function f(t, y(t), y(t-τ))
    DelayFuncMulti delay_derivative_multi = nullptr; // Multi-delay DDE function f(t, y(t), {y(t-τᵢ)})
    
    // Delay Configuration
    dVec      delay_times;                      // Vector of delay times [τ₁, τ₂, ..., τₙ]
                                               // Must be positive and in ascending order
                                               // For state-dependent delays, these are initial values
    
    // History Function Management
    Func      history_function;                 // History function y(t) for t < t₀
                                               // Must be continuous and match y(t₀) at t₀
                                               // Used to evaluate delayed terms during initialization
                                               // For multiple delays, must cover [-max(τᵢ), t₀]
    
    // DDE Solver Configuration
    bool      is_dde            = false;        // Flag to use DDE-specific solver components
                                               // Enables delay interpolation and history tracking
    
    // Interpolation Control
    size_t    max_delay_order   = 10;          // Maximum interpolation order for delayed terms
                                               // Higher orders provide better accuracy but:
                                               //   - Require more memory for solution history
                                               //   - Increase computational cost
                                               //   - May introduce oscillations (Runge phenomenon)
                                               // Recommended range: 4-10 for most problems
    
    double    interpolation_tol = 1.0e-6;      // Tolerance for delay term interpolation
                                               // Controls accuracy vs. efficiency tradeoff
                                               // Should generally be similar to local_tol
                                               // Critical for avoiding artificial oscillations
                                               // and maintaining solution stability

    // Advanced DDE Features to be added!
    // - Automatic detection of discontinuity propagation
    // - Proper handling of multiple delays crossing each other
    // - Support for state-dependent delays
    // - Efficient storage of solution history
    // - Specialized interpolation near discontinuities
};
```

### `SolverResults`

A structure to store the results of a numerical simulation.

```cpp
struct SolverResults {
    Matrix      solution;       // The solution matrix.
    dVec         time_points;    // Time points where solution is computed.
    dVec         errors;         // Error estimates.
    dVec         step_sizes;     // Actual step sizes used.
    Matrix      steps_history;  // History of step sizes used (contains failed steps attempted).
    Matrix      errors_history; // History of errors (contains failed steps attempted).
    wVec failed_trials;  // Number of failed trials in each step.
};
```

## Sparse Data Structures (`header.hpp`)

These structures provide efficient storage for sparse vectors and matrices.

### `SparsedVec`

A structure for storing sparse 1D vectors using (index, value) pairs.

```cpp
struct SparsedVec {
    std::vector<std::pair<size_t, double>> values; // (index, value) pairs
    SparsedVec(size_t N) : values(N) {}
    void add(size_t i, double value) { values.emplace_back(i, value); }
};
```

### `SparseMatrix`

A structure for storing sparse 2D matrices, where each row contains (column, value) pairs for non-zero entries.

```cpp
struct SparseMatrix {
    std::vector<std::vector<std::pair<size_t, double>>> rows; // Each row: (col, value) pairs
    SparseMatrix(size_t N) : rows(N) {}
    void add(size_t i, size_t j, double weight) { rows[i].emplace_back(j, weight); }
};
```

## Complex Number Type (`complex.hpp`)

This section defines a custom `Complex` number implementation designed for high-precision mathematical computations. The structure provides comprehensive support for complex arithmetic, transcendental functions, and special functions, with careful attention to numerical stability and branch cuts.

### Mathematical Background

Complex numbers extend the real number system by introducing the imaginary unit $i$, where $i^2 = -1$. A complex number $z$ is represented in rectangular (Cartesian) form as:

$z = x + yi$

where $x$ is the real part and $y$ is the imaginary part. Alternatively, it can be represented in polar form:

$z = r e^{iθ} = r(\cos θ + i\sin θ)$

where:
- $r = \sqrt{x^2 + y^2}$ is the magnitude (absolute value)
- $θ = \arctan2(y, x)$ is the argument (phase angle)

### `Complex` Structure

Implements a complex number using double-precision floating-point components, optimized for both numerical accuracy and computational efficiency. Special attention is paid to proper handling of:
- Branch cuts for multi-valued functions
- Overflow and underflow conditions
- NaN and Infinity propagation
- Numerical stability near the origin and at large magnitudes

```cpp
struct Complex {
    double re, im;

    // Constructors
    Complex(double r = 0.0, double i = 0.0);

    // Operators (arithmetic, compound assignment, comparison)
    // Supports operations with other Complex numbers, doubles, and integers.

    // Unary operators
    Complex operator+() const; // Unary plus
    Complex operator-() const; // Unary minus

    // Accessors
    double& operator[](int idx); // Access real (0) or imaginary (1) part
    const double& operator[](int idx) const;

    // Explicit bool conversion
    explicit operator bool() const; // True if not (0,0)

    // Basic functions
    Complex exp() const; // e^z
    Complex log() const; // Natural logarithm
    Complex pow(const Complex& o) const; // z^w
    Complex pow(double p) const; // z^p
    Complex log(int branch) const; // Logarithm with branch cut
    Complex pow(const Complex& o, int branch) const; // z^w with branch cut
    Complex pow(double p, int branch) const; // z^p with branch cut
    Complex sqrt() const; // Principal square root
    Complex inverse() const; // 1/z
    Complex root(int n) const; // n-th principal root

    // Polar form
    std::string to_polar() const; // String representation in polar form
    static double abs(const Complex& c); // Magnitude
    static double arg(const Complex& c); // Argument (phase angle)
    static inline Complex from_polar(double r, double theta); // Construct from polar coordinates

    // Template support for other floating-point types
    template<typename T>
    static Complex from(const T& c);

    // Trigonometric functions
    Complex sin() const;
    Complex cos() const;
    Complex tan() const;

    // Inverse trigonometric functions
    Complex asin() const;
    Complex acos() const;
    Complex atan() const;

    // Reciprocal trigonometric functions
    Complex csc() const; // 1/sin(z)
    Complex sec() const; // 1/cos(z)
    Complex cot() const; // 1/tan(z)

    // Inverse reciprocal trigonometric functions
    Complex acsc() const;
    Complex asec() const;
    Complex acot() const;

    // Hyperbolic functions
    Complex sinh() const;
    Complex cosh() const;
    Complex tanh() const;

    // Inverse hyperbolic functions
    Complex asinh() const;
    Complex acosh() const;
    Complex atanh() const;

    // Reciprocal hyperbolic functions
    Complex csch() const; // 1/sinh(z)
    Complex sech() const; // 1/cosh(z)
    Complex coth() const; // 1/tanh(z)

    // Inverse reciprocal hyperbolic functions
    Complex acsch() const;
    Complex asech() const;
    Complex acoth() const;

    // Miscellaneous
    Complex conj() const; // Conjugate
    double norm() const; // Magnitude (same as abs())
    double abs() const; // Magnitude
    double arg() const; // Argument (phase angle)
    double real() const; // Real part
    double imag() const; // Imaginary part
    double phase_deg() const; // Phase angle in degrees
    bool is_nan() const; // Check for NaN components
    bool is_inf() const; // Check for Inf components

    // Static constants
    static const Complex I; // Imaginary unit (0, 1)
    static const Complex Zero; // (0, 0)
    static const Complex One; // (1, 0)
};
```

### Core Operations

#### Arithmetic Operations
- Addition and subtraction: Direct component-wise operations
- Multiplication: $(a + bi)(c + di) = (ac - bd) + (ad + bc)i$
- Division: $\frac{a + bi}{c + di} = \frac{(ac + bd) + (bc - ad)i}{c^2 + d^2}$
  - Implements Smith's algorithm for improved numerical stability

#### Essential Functions
- Conjugate: $\overline{z} = x - yi$
- Magnitude: $|z| = \sqrt{x^2 + y^2}$
  - Uses hypot() for better numerical stability
- Argument (phase): $\arg(z) = \arctan2(y, x)$
  - Principal value in $(-π, π]$

### Transcendental Functions

#### Exponential and Logarithm
- Complex exponential: $e^{x+yi} = e^x(\cos y + i\sin y)$
- Complex logarithm: $\ln(z) = \ln|z| + i\arg(z)$
  - Branch cut along negative real axis
  - Principal branch: $-π < \Im(\ln z) ≤ π$

#### Trigonometric Functions
- Sine: $\sin(z) = \frac{e^{iz} - e^{-iz}}{2i}$
- Cosine: $\cos(z) = \frac{e^{iz} + e^{-iz}}{2}$
- Tangent: $\tan(z) = \frac{\sin(z)}{\cos(z)}$
  - Careful handling near poles at $z = \frac{π}{2} + πn$

#### Inverse Trigonometric Functions
- Arcsine: Principal branch in $[-\frac{π}{2}, \frac{π}{2}]$
- Arccosine: Principal branch in $[0, π]$
- Arctangent: Principal branch in $(-\frac{π}{2}, \frac{π}{2})$

#### Hyperbolic Functions
- Sinh: $\sinh(z) = \frac{e^z - e^{-z}}{2}$
- Cosh: $\cosh(z) = \frac{e^z + e^{-z}}{2}$
- Tanh: $\tanh(z) = \frac{\sinh(z)}{\cosh(z)}$

#### Power Functions
- Integer powers: Optimized multiplication sequence
- Complex powers: $z^w = e^{w\ln(z)}$
  - Branch cut handling for fractional powers
- Square root: Principal branch with cut along negative real axis

### Utility Functions
- String formatting with controllable precision
- Conversion to/from string representations
- NaN and Infinity checking
- Comparisons and equality testing with tolerance

### Stream Operators

- `operator<<`: Overloaded output stream operator for `Complex` numbers, printing in the format `(re+imi)`.
- `operator>>`: Overloaded input stream operator for `Complex` numbers.
