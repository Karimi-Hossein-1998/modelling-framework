# Numerical Methods Mathematical Background

## 1. One-Step Methods

### 1.1 General Framework

One-step methods for solving ODEs have the general form:

$$\mathbf{y}_{n+1} = \mathbf{y}_n + h\Phi(t_n, \mathbf{y}_n, h)$$

where $\Phi$ is the increment function specific to each method.

### 1.2 Runge-Kutta Methods

General s-stage explicit Runge-Kutta method:

$$\mathbf{k}_i = \mathbf{f}\left(t_n + c_ih, \mathbf{y}_n + h\sum_{j=1}^{i-1} a_{ij}\mathbf{k}_j\right), \quad i=1,\ldots,s$$
$$\mathbf{y}_{n+1} = \mathbf{y}_n + h\sum_{i=1}^s b_i\mathbf{k}_i$$

#### Order Conditions

For order p, all equations up to order p must be satisfied:

1. First order: $\sum_{i=1}^s b_i = 1$
2. Second order: $\sum_{i=1}^s b_ic_i = \frac{1}{2}$
3. Third order:
   - $\sum_{i=1}^s b_ic_i^2 = \frac{1}{3}$
   - $\sum_{i,j=1}^s b_ia_{ij}c_j = \frac{1}{6}$

### 1.3 Specific Methods

#### 1.3.1 Classical RK4

Butcher tableau:

$$\large
\begin{array}{c|c}
\mathbf{c} & \mathbf{A} \\[0.3cm]
\hline %\\[0.1cm]
\begin{matrix}
0                  \\[0.1cm]
\frac{1}{2}       \\[0.1cm]
\frac{1}{2}       \\[0.1cm]
1                  %\\[0.3cm]
\end{matrix}
&
\begin{matrix}
0           & 0           & 0      & 0      \\[0.1cm]
\frac{1}{2} & 0           & 0      & 0      \\[0.1cm]
0           & \frac{1}{2} & 0      & 0      \\[0.1cm]
0           & 0           & 1      & 0      %\\[0.3cm]
\end{matrix} \\[0.3cm]
\hline %\\[0.3cm]
\mathbf{b} &
\begin{matrix}
\frac{1}{6} & \frac{1}{3} & \frac{1}{3} & \frac{1}{6}
\end{matrix}
\end{array}
$$

Local truncation error:
$$\tau_n = \frac{1}{720}h^4\mathbf{y}^{(5)}(\xi_n)$$

#### 1.3.2 Embedded Pairs

For error estimation in adaptive methods:

$$\|\mathbf{e}_n\| \approx \|h\sum_{i=1}^s (b_i - b_i^*)\mathbf{k}_i\|$$

## 2. Linear Multistep Methods

### 2.1 General Form

$$\sum_{j=0}^k \alpha_j\mathbf{y}_{n+j} = h\sum_{j=0}^k \beta_j\mathbf{f}_{n+j}$$

where:

- $\alpha_k = 1$ (normalization)
- $\{\alpha_j\}$ and $\{\beta_j\}$ define the specific method

### 2.2 Adams Methods

#### 2.2.1 Adams-Bashforth (Explicit)

$$\mathbf{y}_{n+1} = \mathbf{y}_n + h\sum_{j=0}^{k-1} \beta_j\mathbf{f}_{n-j}$$

Coefficients derived from integration of interpolating polynomial:
$$\beta_j = (-1)^j\int_0^1 \prod_{i=0,i\neq j}^{k-1} \frac{s+i}{j-i}ds$$

#### 2.2.2 Adams-Moulton (Implicit)

$$\mathbf{y}_{n+1} = \mathbf{y}_n + h\left(\beta_k \mathbf{f}_{n+1} + \sum_{j=0}^{k-1} \beta_j\mathbf{f}_{n-j}\right)$$

### 2.3 Stability Analysis

#### 2.3.1 Linear Stability Analysis

For the test equation $\mathbf{y}' = \lambda\mathbf{y}$, define:

1. **Characteristic Polynomial**:
   $$\rho(\zeta) = \sum_{j=0}^k \alpha_j\zeta^j$$

2. **Stability Polynomial**:
   $$\pi(\zeta,q) = \rho(\zeta) - q\sigma(\zeta)$$
   where $q = h\lambda$ and $\sigma(\zeta) = \sum_{j=0}^k \beta_j\zeta^j$

3. **Root Condition**:
   All roots $\zeta$ of $\pi(\zeta,q) = 0$ must satisfy $|\zeta| \leq 1$,
   with $|\zeta| = 1$ only if $\zeta$ is a simple root.

#### 2.3.2 Stability Regions

1. **Definition**:
   $$S = \{q \in \mathbb{C} : \text{method is stable for } h\lambda = q\}$$

2. **Boundary Locus**:
   $$\partial S = \{q \in \mathbb{C} : \pi(e^{i\theta},q) = 0 \text{ for some } \theta \in [0,2\pi]\}$$

3. **A-Stability**:
   Method is A-stable if $S \supset \mathbb{C}^-$

#### 2.3.3 Practical Stability Considerations

1. **Stiff Problems**:
   - Use A-stable or A($\alpha$)-stable methods
   - Monitor stability function evaluation
   - Consider L-stability for severe stiffness

2. **Variable Step Size**:
   - Recompute stability parameters when step size changes
   - Use stability-based step size control
   - Monitor characteristic roots

## 3. Error Analysis and Control

### 3.1 Local Error Estimation

For embedded RK pairs (p,p-1):
$$\text{err} = \|\mathbf{y}_{n+1} - \hat{\mathbf{y}}_{n+1}\|$$

### 3.2 Step Size Control

New step size calculation:
$$h_{\text{new}} = h_{\text{old}}\min\left(\text{facmax},\max\left(\text{facmin},\text{fac}\cdot(\text{tol}/\text{err})^{1/p}\right)\right)$$

where:

- `facmin`: Minimum reduction factor (e.g., 0.1)
- `facmax`: Maximum increase factor (e.g., 5.0)
- `fac`: Safety factor (e.g., 0.9)
- `p`: Method order

### 3.3 Error Propagation

For linear multistep methods:
$$\mathbf{e}_{n+k} + \sum_{j=0}^{k-1} \alpha_j\mathbf{e}_{n+j} = h\sum_{j=0}^k \beta_j\left(\mathbf{f}(t_{n+j},\mathbf{y}(t_{n+j})) - \mathbf{f}(t_{n+j},\mathbf{y}_{n+j})\right) + h\tau_{n+k}$$

## 4. Implementation Considerations

### 4.1 Method Selection Guidelines

1. **Non-stiff Problems**:
   - Explicit RK methods (RK4, Dormand-Prince)
   - Adams-Bashforth methods
   - Consider embedded pairs for adaptivity

2. **Stiff Problems**:
   - Implicit RK methods (Radau IIA, SDIRK)
   - Adams-Moulton methods
   - BDF methods for severe stiffness

### 4.2 Algorithmic Optimizations

#### 4.2.1 Memory Management

```cpp
struct SolverState {
    vector<Vector> y_history;    // Solution history
    vector<Vector> f_history;    // Function evaluation history
    vector<double> t_history;    // Time points
    
    // Circular buffer implementation
    void push(const Vector& y, const Vector& f, double t) {
        rotate(y_history.begin(), y_history.end()-1, y_history.end());
        y_history[0] = y;
        // Similar for f_history and t_history
    }
};
```

#### 4.2.2 Efficient Implementation Strategies

1. **Stage Value Computations**:
   - Reuse common terms in stage calculations
   - Implement FSAL (First Same As Last) optimization
   - Use workspace arrays effectively

2. **Linear Algebra Operations**:
   - Optimize matrix-vector products
   - Use specialized solvers for linear systems
   - Implement careful scaling strategies

### 4.3 Advanced Error Control

#### 4.3.1 PI Step Size Control

More sophisticated step size control using PI-controller:

$$h_{n+1} = h_n\left(\frac{\text{tol}}{\text{err}_n}\right)^{0.6/p}\left(\frac{\text{err}_{n-1}}{\text{err}_n}\right)^{0.2/p}$$

#### 4.3.2 Error Estimation Strategies

1. **Local Error Analysis**:
   $$\text{err}_{\text{local}} = \|\mathbf{y}_{n+1} - \hat{\mathbf{y}}_{n+1}\|_{\text{scaled}}$$
   where
   $$\|\mathbf{v}\|_{\text{scaled}} = \sqrt{\frac{1}{d}\sum_{i=1}^d\left(\frac{v_i}{\text{atol}_i + \text{rtol}_i|y_i|}\right)^2}$$

2. **Global Error Control**:
   - Use Richardson extrapolation
   - Implement deferred correction
   - Monitor error accumulation

### 4.4 Special Cases

1. **Mass Matrices**:
   $$\mathbf{M}\frac{d\mathbf{y}}{dt} = \mathbf{f}(t,\mathbf{y})$$

   Require special treatment for:
   - Singular mass matrices
   - Time-dependent mass matrices
   - Structure preservation

2. **Conservation Laws**:
   - Monitor conserved quantities
   - Use geometric integrators
   - Implement projection methods

## References

1. Hairer, E., Nørsett, S. P., & Wanner, G. (1993). Solving Ordinary Differential Equations I.
2. Butcher, J. C. (2008). Numerical Methods for Ordinary Differential Equations.
3. Ascher, U. M., & Petzold, L. R. (1998). Computer Methods for Ordinary Differential Equations.
4. Hairer, E., & Wanner, G. (1996). Solving Ordinary Differential Equations II.
5. Gustafsson, K. (1991). Control theoretic techniques for stepsize selection in explicit Runge-Kutta methods.
6. Söderlind, G. (2002). Automatic control and adaptive time-stepping.
