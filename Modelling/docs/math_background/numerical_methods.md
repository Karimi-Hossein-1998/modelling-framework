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
0 \\[0.1cm]
\frac{1}{2} \\[0.1cm]
\frac{1}{2} \\[0.1cm]
1 %\\[0.3cm]
\end{matrix}
&
\begin{matrix}
0      & 0      & 0      & 0      \\[0.1cm]
\frac{1}{2} & 0      & 0      & 0      \\[0.1cm]
0      & \frac{1}{2} & 0      & 0      \\[0.1cm]
0      & 0      & 1      & 0 %\\[0.3cm]
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

$$\mathbf{y}_{n+1} = \mathbf{y}_n + h\sum_{j=-1}^{k-1} \beta_j\mathbf{f}_{n-j}$$

### 2.3 Stability Analysis

#### 2.3.1 Characteristic Polynomial

$$\rho(\zeta) = \sum_{j=0}^k \alpha_j\zeta^j$$

#### 2.3.2 Stability Region

For the test equation $\mathbf{y}' = \lambda\mathbf{y}$:

$$\sum_{j=0}^k (\alpha_j - h\lambda\beta_j)\zeta^j = 0$$

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

## References

1. Hairer, E., Nørsett, S. P., & Wanner, G. (1993). Solving Ordinary Differential Equations I.
2. Butcher, J. C. (2008). Numerical Methods for Ordinary Differential Equations.
3. Ascher, U. M., & Petzold, L. R. (1998). Computer Methods for Ordinary Differential Equations.
