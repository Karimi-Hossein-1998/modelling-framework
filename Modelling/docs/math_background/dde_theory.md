# Delay Differential Equations Mathematical Background

## 1. Basic Theory

### 1.1 General Form

A general DDE with multiple delays can be written as:

$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}(t), \{\mathbf{y}(t-\tau_i)\}_{i=1}^m)$$

with initial history:
$$\mathbf{y}(t) = \mathbf{\phi}(t), \quad t \leq t_0$$

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

### 3.2 Runge-Kutta Methods for DDEs

#### 3.2.1 General Framework

For a single delay $\tau$:

$$\mathbf{k}_i = \mathbf{f}(t_n + c_ih, \mathbf{Y}_i, \mathbf{Z}_i)$$
$$\mathbf{Y}_i = \mathbf{y}_n + h\sum_{j=1}^{i-1} a_{ij}\mathbf{k}_j$$
$$\mathbf{Z}_i = \text{interpolation}(t_n + c_ih - \tau)$$

#### 3.2.2 Interpolation Requirements

1. **Order**: Interpolant order should match method order
2. **Continuity**: Must be continuous across steps
3. **Efficiency**: Minimal number of evaluations

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

For the test equation:
$$\frac{dy}{dt} = ay(t) + by(t-\tau)$$

Characteristic equation:
$$\lambda - a - be^{-\lambda\tau} = 0$$

### 5.2 Numerical Stability

Modified stability region depends on:

1. Basic method stability
2. Interpolation procedure
3. Treatment of delayed terms

## References

1. Bellen, A., & Zennaro, M. (2013). Numerical Methods for Delay Differential Equations.
2. Baker, C. T. H., et al. (1995). Algorithm 769: Software for Delay Differential Equations.
3. Shampine, L. F., & Thompson, S. (2001). Solving DDEs in MATLAB.
