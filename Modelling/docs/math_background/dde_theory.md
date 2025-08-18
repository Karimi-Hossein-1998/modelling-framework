# Delay Differential Equations Mathematical Background

## 1. Basic Theory

### 1.1 General Form

A general delay differential equation (DDE) with multiple delays can be expressed in the following form:

$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}(t), \{\mathbf{y}(t-\tau_i)\}_{i=1}^m), \quad t > t_0$$

where:

- $\mathbf{y}(t) \in \mathbb{R}^n$ is the state vector
- $\mathbf{f}: \mathbb{R} \times \mathbb{R}^n \times (\mathbb{R}^n)^m \rightarrow \mathbb{R}^n$ is the vector field
- $\tau_i > 0$ are the delay times
- $m$ is the number of distinct delays

The initial condition is specified through a history function:
$$\mathbf{y}(t) = \mathbf{\phi}(t), \quad t \in [t_0-r, t_0]$$

where:

- $r = \max_{i}\{\tau_i\}$ is the maximal delay
- $\mathbf{\phi}: [t_0-r, t_0] \rightarrow \mathbb{R}^n$ is the history function

### 1.2 Classification of DDEs

1. **Constant Delay DDEs**:
   $$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}(t), \mathbf{y}(t-\tau))$$
   where $\tau > 0$ is constant.

2. **State-Dependent Delay DDEs**:
   $$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}(t), \mathbf{y}(t-\tau(\mathbf{y}(t))))$$
   where $\tau$ depends on the solution.

3. **Distributed Delay DDEs**:
   $$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}(t), \int_{-\infty}^t K(t-s)\mathbf{y}(s)ds)$$
   where $K$ is a kernel function.

## 2. Theoretical Framework

### 2.1 Solution Space

The natural space for DDEs is $C([-r,0],\mathbb{R}^n)$, where:

- $r = \max_i\{\tau_i\}$ is the maximum delay
- The space consists of continuous functions on $[-r,0]$

### 2.2 Existence and Uniqueness

For the initial value problem:

1. **Existence**: Requires
   - $\mathbf{f}$ is continuous
   - $\mathbf{\phi}$ is continuous
   - Matching conditions at $t_0$

2. **Uniqueness**: Requires
   - Local Lipschitz condition on $\mathbf{f}$
   - Continuity of delays (for state-dependent case)

## 3. Numerical Methods

### 3.1 Method of Steps

1. Start with history function $\mathbf{\phi}(t)$ for $t \leq t_0$
2. Solve as ODE on $[t_0, t_0+\tau]$ using history for delayed terms
3. Extend solution step by step

### 3.2 Implementation of Delay Interpolation

The codebase implements delay interpolation through two main functions:

#### 3.2.1 Single Delay Interpolation

```cpp
inline dVec interpolate_delay(
    double        t_delayed,     // Time to interpolate at (t - tau)
    const dVec&   time_points,   // Available time points
    const Matrix& solution,      // Solution at those time points
    double        t0,           // Initial time
    const Func&   history_function // Pre-history function
)
```

Key features:

1. Linear interpolation between available points
2. History function handling for t < t0
3. Efficient bracketing point search

#### 3.2.2 Multiple Delay Interpolation

```cpp
inline Matrix interpolate_multiple_delays(
    const dVec&   delay_times,   // Vector of delay times τ_i
    double        t,             // Current time
    const dVec&   time_points,   // Available time points
    const Matrix& solution,      // Solution at those time points
    double        t0,            // Initial time
    const Func&   history_function // Pre-history function
)
```

Implementation details:

1. Handles vector of delays efficiently
2. Returns matrix of delayed states
3. Reuses single delay interpolation

### 3.3 Runge-Kutta Methods for DDEs

#### 3.3.1 Implementation Architecture

The framework implements both single-delay and multiple-delay RK methods:

1. **Single Delay Interface**:

   ```cpp
   inline SolverResults rk4_dde_solver(const SolverParameters& Params)
   {
       const auto& f = Params.delay_derivative;        // f(t, y, y_delayed)
       const auto& delay_times = Params.delay_times;   // Single delay
       const auto& history_function = Params.history_function;
   }
   ```

2. **Multiple Delay Interface**:

   ```cpp
   inline SolverResults rk4_dde_general_solver(const SolverParameters& Params)
   {
       const auto& f = Params.delay_derivative_multi;  // f(t, y, y_delayed_states)
       const auto& delay_times = Params.delay_times;   // Vector of delays
       const auto& history_function = Params.history_function;
   }
   ```

#### 3.3.2 General Framework

For a system with multiple delays $\{\tau_i\}_{i=1}^m$, a $s$-stage Runge-Kutta method takes the form:

$$\mathbf{k}_i = \mathbf{f}(t_n + c_ih, \mathbf{Y}_i, \{\mathbf{Z}_i^j\}_{j=1}^m), \quad i=1,\ldots,s$$

where:

- $\mathbf{Y}_i    = \mathbf{y}_n + h\sum_{j=1}^{i-1} a_{ij}\mathbf{k}_j$    are the stage values
- $\mathbf{Z}_i^j  = \mathcal{I}(t_n + c_ih - \tau_j)$                       are the delayed terms
- $\mathcal{I}(t)$ is an interpolant constructed from previous solution values

#### 3.3.3 Implementation Details

1. **Memory Management**:

   ```cpp
   // Pre-allocate working memory
   auto y_temp = dVec(N, 0.0);
   dVec k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
   Matrix y_delayed_states(delay_times.size(), dVec(N, 0.0));
   ```

2. **Stage Evaluations (RK4)**:

   ```cpp
   // Stage 1
   y_delayed_states = interpolate_multiple_delays(delay_times, t, 
       time_points, solution, t0, history_function);
   k1 = f(t, y, y_delayed_states);
   
   // Stage 2
   for (size_t j = 0; j < N; ++j)
       y_temp[j] = y[j] + dt_half * k1[j];
   y_delayed_states = interpolate_multiple_delays(delay_times, 
       t + dt_half, time_points, solution, t0, history_function);
   k2 = f(t + dt_half, y_temp, y_delayed_states);
   ```

3. **Solution Update**:

   ```cpp
   for (size_t j = 0; j < N; ++j)
       y[j] += dt_sixth * (k1[j] + k4[j]) + dt_third * (k2[j] + k3[j]);
   ```

The solution is advanced using:

$$\mathbf{y}_{n+1} = \mathbf{y}_n + h\sum_{i=1}^s b_i\mathbf{k}_i$$

For methods of order $p$, the interpolant $\mathcal{I}(t)$ must satisfy:

1. $\|\mathcal{I}(t) - \mathbf{y}(t)\| = \mathcal{O}(h^p)$ for $t \in [t_n-\tau, t_{n+1}]$
2. $\mathcal{I}(t_j) = \mathbf{y}_j$ at mesh points
3. $\mathcal{I}(t)$ is continuous across step boundaries

#### 3.2.2 Interpolation Implementation

The codebase implements two key interpolation functions in `dde_utility.hpp`:

1. **Single Delay Interpolation**:

   ```cpp
   inline dVec interpolate_delay(
       double t_delayed,        // Time point t - τ
       const dVec& time_points, // Available solution times
       const Matrix& solution,  // Solution history
       double t0,              // Initial time
       const Func& history_function // Pre-history values
   )
   ```

2. **Multiple Delay Interpolation**:

   ```cpp
   inline Matrix interpolate_multiple_delays(
       const dVec& delay_times,  // Delay times {τ₁, ..., τₘ}
       double t,                // Current time
       const dVec& time_points, // Available solution times
       const Matrix& solution,  // Solution history
       double t0,              // Initial time
       const Func& history_function // Pre-history values
   )
   ```

Implementation features:

1. **Pre-history Handling**:

   ```cpp
   if (t_delayed < t0)
   {
       return history_function(t_delayed);
   }
   ```

2. **Linear Interpolation**:

   ```cpp
   const double alpha = (t_delayed - t1) / (t2 - t1);
   for (size_t i = 0; i < N; ++i)
   {
       y_delayed[i] = (1.0 - alpha) * solution[idx][i] + 
                     alpha * solution[idx + 1][i];
   }
   ```

3. **Multiple Delay Efficiency**:
   - Reuses single delay interpolation
   - Returns matrix of delayed states
   - Each row corresponds to state at t - τᵢ

### 3.3 Discontinuity Tracking

#### 3.3.1 Primary Discontinuities

At $t_0$: Junction of history and solution

#### 3.3.2 Propagated Discontinuities

At times $t_k = t_0 + k\tau$:

- Due to propagation of primary discontinuity
- Order typically reduces by 1 at each propagation

### 3.4 Step Size Control

Modified error estimation for DDEs:

$$\text{err} = \max\left(\frac{\|\mathbf{y}_{n+1} - \hat{\mathbf{y}}_{n+1}\|}{\text{tol}}, \frac{\|\mathbf{z}_{n+1} - \hat{\mathbf{z}}_{n+1}\|}{\text{tol}}\right)$$

where $\mathbf{z}$ represents interpolated values.

## 4. Special Considerations

### 4.1 Multiple Delays

For systems with multiple delays $\{\tau_1,\ldots,\tau_m\}$:

1. **Memory Management**:
   - Store solution history efficiently
   - Minimize interpolation overhead

2. **Step Size Selection**:
   - Consider all delay times
   - Avoid very small steps due to delay ratios

### 4.2 State-Dependent Delays

Additional challenges:

1. Location of delayed points requires iteration
2. Error in delay affects solution accuracy
3. Modified error estimation needed

## 5. Stability Analysis

### 5.1 Linear Stability

Consider the scalar linear test equation:

$$\frac{dy}{dt} = ay(t) + by(t-\tau), \quad t > 0$$
$$y(t) = \phi(t), \quad t \in [-\tau, 0]$$

where $a, b \in \mathbb{R}$ and $\tau > 0$.

The characteristic equation is:

$$\Delta(\lambda) := \lambda - a - be^{-\lambda\tau} = 0$$

Properties of the characteristic equation:

1. It has infinitely many complex roots
2. Only finitely many roots have real parts greater than any given number
3. The stability of the zero solution is determined by the rightmost root

#### 5.1.1 Stability Regions

For $\tau > 0$, the stability region in the $(a,b)$-plane is given by:

1. When $a < 0$:
   $$|b| < -a$$

2. When $a > 0$:
   $$|b| < -\frac{\pi}{2\tau}$$

### 5.2 Numerical Stability

The numerical stability of a DDE solver depends on three key components:

1. **Basic Method Stability**:
   - The stability of the underlying RK or multistep method
   - Analysis based on the method's stability function $R(z)$

2. **Interpolation Influence**:
   - The stability properties of the interpolation scheme
   - Effect of interpolation errors on solution growth

3. **Delay Treatment**:
   - Impact of delay ratio $h/\tau$ on stability
   - Modifications needed for multiple delays

#### 5.2.1 Practical Implications

For numerical implementation:

1. Choose step size $h$ satisfying:
   $$h \leq \min\left\{\frac{\tau}{m}, \frac{2}{|a|+|b|}\right\}$$
   where $m \geq 1$ is a safety factor

2. For multiple delays $\{\tau_i\}_{i=1}^m$:
   - Use adaptive mesh refinement near discontinuity points
   - Consider delay ratios in step size selection

## 6. Convergence Analysis and Error Control

### 6.1 Global Error Analysis

For a numerical method of order $p$, the global error can be decomposed as:

$$\mathbf{e}(t) = \mathbf{y}(t) - \mathbf{y}_h(t) = \mathbf{e}_d(t) + \mathbf{e}_c(t)$$

where:

- $\mathbf{e}_d(t)$ is the discretization error ($\mathcal{O}(h^p)$)  
- $\mathbf{e}_c(t)$ is the propagated error from continuous extension  

#### 6.1.1 Error Propagation

The error satisfies a variational equation:

$$\frac{d\mathbf{e}}{dt} = \mathbf{L}(t)\mathbf{e}(t) + \mathbf{M}(t)\mathbf{e}(t-\tau) + \mathbf{\delta}(t)$$

where:

- $\mathbf{L}(t)      = \frac{\partial \mathbf{f}}{\partial \mathbf{y}}(t,\mathbf{y}(t),\mathbf{y}(t-\tau))$  
- $\mathbf{M}(t)      = \frac{\partial \mathbf{f}}{\partial \mathbf{z}}(t,\mathbf{y}(t),\mathbf{y}(t-\tau))$  
- $\mathbf{\delta}(t)$ is the local truncation error

### 6.2 Adaptive Error Control

#### 6.2.1 Step Size Selection

The optimal step size is determined by:

$$h_{\text{new}} = h_{\text{old}}\min\left(\text{facmax},\max\left(\text{facmin},\text{fac}\cdot(\frac{\text{tol}}{\text{err}})^{\frac{1}{p+1}}\right)\right)$$

where:

- $\text{err} = \max(\text{err}_d, \text{err}_i)$
- $\text{err}_d$ is the discrete solution error
- $\text{err}_i$ is the interpolation error

#### 6.2.2 Error Estimation Strategy

For methods with embedded error estimation:

1. Compute main and embedded solutions:
   $$\mathbf{y}_{n+1} = \mathbf{y}_n + h\sum_{i=1}^s b_i\mathbf{k}_i$$
   $$\hat{\mathbf{y}}_{n+1} = \mathbf{y}_n + h\sum_{i=1}^s \hat{b}_i\mathbf{k}_i$$

2. Estimate local error:
   $$\text{err}_d = \|\mathbf{y}_{n+1} - \hat{\mathbf{y}}_{n+1}\|$$

3. Compute interpolation error:
   $$\text{err}_i = \max_{t \in [t_n,t_{n+1}]} \|\mathcal{I}(t) - \hat{\mathcal{I}}(t)\|$$

### 6.3 Implementation Guidelines

1. **Error Control Parameters**:
   - Set $\text{facmin} = 0.1$ (minimum step size reduction)
   - Set $\text{facmax} = 5.0$ (maximum step size increase)
   - Use $\text{fac} = 0.9$ (safety factor)

2. **Discontinuity Handling**:
   - Locate discontinuities to within $\mathcal{O}(h^{p+1})$
   - Restart method at discontinuity points
   - Use special starting procedures after discontinuities

3. **Memory Management**:
   - Store minimal solution history for interpolation
   - Implement efficient interpolation search algorithms
   - Use circular buffers for history storage

## References

1. Bellen, A., & Zennaro, M. (2013). Numerical Methods for Delay Differential Equations.
2. Baker, C. T. H., et al. (1995). Algorithm 769: Software for Delay Differential Equations.
3. Shampine, L. F., & Thompson, S. (2001). Solving DDEs in MATLAB.
4. Paul, C. A. H. (1992). Developing a Delay Differential Equation Solver.
5. Guglielmi, N., & Hairer, E. (2001). Implementing Radau IIA Methods for Delay Differential Equations.
