# Mathematical Foundations

This document establishes the mathematical notation and conventions used throughout the MyLib documentation. All modules follow these conventions for consistency.

## Notation Conventions

### General Mathematical Notation

1. **Scalars**:
   - Lower-case letters ($x$, $y$, $t$)
   - Greek letters for parameters ($α$, $β$, $τ$)

2. **Vectors**:
   - Bold lower-case letters ($\mathbf{x}$, $\mathbf{y}$)
   - Components denoted with subscripts ($x_i$)
   - Column vectors by default

3. **Matrices**:
   - Bold upper-case letters ($\mathbf{A}$, $\mathbf{B}$)
   - Elements denoted with double subscripts ($a_{ij}$)

4. **Functions**:
   - Standard font for named functions ($\sin$, $\cos$)
   - Calligraphic letters for operators ($\mathcal{L}$)
   - Arguments in parentheses: $f(x)$, $\mathbf{F}(\mathbf{x})$

### Numerical Analysis Conventions

1. **Discretization**:
   - Time points: $t_n = t_0 + nh$
   - Step size: $h$ (temporal) or $\Delta x$ (spatial)
   - Numerical solution: $y_n \approx y(t_n)$

2. **Error Terms**:
   - Local truncation error: $\tau_n$
   - Global error: $e_n = y(t_n) - y_n$
   - Big-O notation: $\mathcal{O}(h^p)$

3. **Convergence**:
   - Order of convergence: $p$
   - Error constant: $C$
   - Asymptotic behavior: $\|e_n\| \leq Ch^p$

### Linear Algebra Conventions

1. **Vector Spaces**:
   - Real vector space: $\mathbb{R}^n$
   - Complex vector space: $\mathbb{C}^n$
   - Inner product: $\langle \mathbf{x}, \mathbf{y} \rangle$

2. **Matrix Operations**:
   - Transpose: $\mathbf{A}^T$
   - Inverse: $\mathbf{A}^{-1}$
   - Determinant: $\det(\mathbf{A})$ or $|\mathbf{A}|$

3. **Norms**:
   - Vector p-norm: $\|\mathbf{x}\|_p$
   - Matrix norm: $\|\mathbf{A}\|$
   - Frobenius norm: $\|\mathbf{A}\|_F$

### Differential Equations

1. **ODEs**:
   - First-order system: $\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y})$
   - Initial condition: $\mathbf{y}(t_0) = \mathbf{y}_0$

2. **DDEs**:
   - Single delay: $\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}(t), \mathbf{y}(t-τ))$
   - Multiple delays: $\frac{d\mathbf{y}}{dt} = \mathbf{f}(t, \mathbf{y}(t), \{\mathbf{y}(t-τ_i)\}_{i=1}^m)$
   - History function: $\mathbf{y}(t) = \mathbf{\phi}(t)$ for $t \leq t_0$

## Implementation Mappings

This section maps mathematical notation to code implementations:

| Mathematical Notation | Code Implementation | Type |
|---------------------|---------------------|------|
| $x \in \mathbb{R}$ | `double` | Scalar |
| $\mathbf{x} \in \mathbb{R}^n$ | `dVec` | Vector |
| $\mathbf{A} \in \mathbb{R}^{m×n}$ | `Matrix` | Matrix |
| $f: \mathbb{R} \to \mathbb{R}$ | `func` | Scalar function |
| $\mathbf{F}: \mathbb{R} \to \mathbb{R}^n$ | `Func` | Vector function |
| $\frac{d\mathbf{y}}{dt} = \mathbf{f}(t,\mathbf{y})$ | `MyFunc` | ODE right-hand side |
| $\mathbf{f}(t,\mathbf{y},\mathbf{y}_τ)$ | `DelayFunc` | DDE right-hand side |

## Core Mathematical Concepts

### 1. Function Approximation Theory

#### 1.1 Taylor Series Expansion

The foundation of many numerical methods relies on Taylor series expansions:

For a function $f \in C^{n+1}[a,b]$, around a point $x_0$:
$$f(x) = \sum_{k=0}^n \frac{f^{(k)}(x_0)}{k!}(x-x_0)^k + R_n(x)$$

where the remainder term $R_n(x)$ is:
$$R_n(x) = \frac{f^{(n+1)}(\xi)}{(n+1)!}(x-x_0)^{n+1}$$
for some $\xi$ between $x_0$ and $x$.

#### 1.2 Interpolation Theory

For points $(x_i, y_i)$, $i=0,\ldots,n$, the interpolation problem seeks a function $p(x)$ such that:
$$p(x_i) = y_i, \quad i=0,\ldots,n$$

Key forms include:

1. Lagrange Form:
   $$p(x) = \sum_{i=0}^n y_i \prod_{j\neq i}\frac{x-x_j}{x_i-x_j}$$

2. Newton Form:
   $$p(x) = a_0 + a_1(x-x_0) + a_2(x-x_0)(x-x_1) + \cdots$$
   where $a_i$ are divided differences.

### 2. Linear Algebra Foundations

#### 2.1 Vector Spaces

A vector space $V$ over a field $\mathbb{F}$ satisfies:

1. Closure under addition and scalar multiplication
2. Distributive properties
3. Existence of zero vector and additive inverses

Key spaces used:

- $\mathbb{R}^n$: n-dimensional real vector space
- $C[a,b]$: continuous functions on $[a,b]$

#### 2.2 Inner Product Spaces

For vectors $\mathbf{x}, \mathbf{y} \in V$, inner product $\langle \cdot,\cdot \rangle$ satisfies:

1. Symmetry: $\langle \mathbf{x},\mathbf{y} \rangle = \overline{\langle \mathbf{y},\mathbf{x} \rangle}$
2. Linearity: $\langle \alpha\mathbf{x}+\mathbf{y},\mathbf{z} \rangle = \alpha\langle \mathbf{x},\mathbf{z} \rangle + \langle \mathbf{y},\mathbf{z} \rangle$
3. Positive definiteness: $\langle \mathbf{x},\mathbf{x} \rangle > 0$ if $\mathbf{x} \neq \mathbf{0}$

### 3. Differential Equations Theory

#### 3.1 Existence and Uniqueness

For the initial value problem (IVP):
$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t,\mathbf{y}), \quad \mathbf{y}(t_0) = \mathbf{y}_0$$

Picard-Lindelöf Theorem states that a unique solution exists if:

1. $\mathbf{f}$ is Lipschitz continuous in $\mathbf{y}$
2. $\mathbf{f}$ is continuous in $t$

For DDEs:
$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t,\mathbf{y}(t),\mathbf{y}(t-\tau))$$
Additional conditions on the history function $\phi(t)$ are needed.

#### 3.2 Stability Theory

1. Linear Stability Analysis:
   For linear systems $\frac{d\mathbf{y}}{dt} = \mathbf{A}\mathbf{y}$:
   - Solution stability depends on eigenvalues of $\mathbf{A}$
   - Characteristic equation: $\det(\mathbf{A}-\lambda\mathbf{I}) = 0$

2. Numerical Stability:
   - A-stability: Method stable for all $\lambda h$ in left half-plane
   - L-stability: Additional damping for stiff problems

### 4. Numerical Methods Framework

#### 4.1 Discretization Principles

1. Time Discretization:
   $$t_n = t_0 + nh, \quad n=0,1,\ldots,N$$
   where $h$ is step size.

2. Solution Approximation:
   $$\mathbf{y}_n \approx \mathbf{y}(t_n)$$

#### 4.2 Error Analysis

1. Local Truncation Error:
   $$\tau_n = \frac{\mathbf{y}(t_{n+1}) - \mathbf{y}_n - h\Phi(t_n,\mathbf{y}_n,h)}{h}$$
   where $\Phi$ is the increment function.

2. Global Error:
   $$\mathbf{e}_n = \mathbf{y}(t_n) - \mathbf{y}_n$$
   Relation: $\|\mathbf{e}_n\| \leq Ch^p$

### 5. Matrix Decompositions

#### 5.1 LU Decomposition

For a matrix $\mathbf{A}$:
$$\mathbf{A} = \mathbf{L}\mathbf{U}$$
where:

- $\mathbf{L}$ is lower triangular
- $\mathbf{U}$ is upper triangular

#### 5.2 QR Decomposition

$$\mathbf{A} = \mathbf{Q}\mathbf{R}$$
where:

- $\mathbf{Q}$ is orthogonal ($\mathbf{Q}^T\mathbf{Q} = \mathbf{I}$)
- $\mathbf{R}$ is upper triangular

#### 5.3 Cholesky Decomposition

For symmetric positive definite $\mathbf{A}$:
$$\mathbf{A} = \mathbf{L}\mathbf{L}^T$$
where $\mathbf{L}$ is lower triangular.

### 6. Advanced Topics

#### 6.1 Eigenvalue Problems

For matrix $\mathbf{A}$:
$$\mathbf{A}\mathbf{x} = \lambda\mathbf{x}$$
where:

- $\lambda$ is an eigenvalue
- $\mathbf{x}$ is an eigenvector

Methods include:

1. Power Method:
   $$\mathbf{x}_{k+1} = \frac{\mathbf{A}\mathbf{x}_k}{\|\mathbf{A}\mathbf{x}_k\|}$$

2. QR Algorithm:
   $$\mathbf{A}_k = \mathbf{Q}_k\mathbf{R}_k$$
   $$\mathbf{A}_{k+1} = \mathbf{R}_k\mathbf{Q}_k$$

#### 6.2 Delay Differential Equations

For DDEs with multiple delays:
$$\frac{d\mathbf{y}}{dt} = \mathbf{f}(t,\mathbf{y}(t),\mathbf{y}(t-\tau_1),\ldots,\mathbf{y}(t-\tau_m))$$

Key considerations:

1. History Function Requirements:
   - Continuity at $t_0$
   - Sufficient smoothness

2. Numerical Methods:
   - Interpolation of history
   - Treatment of discontinuity propagation
   - Step size selection near delay points

## References

1. Burden, R. L., & Faires, J. D. (2015). Numerical Analysis (10th ed.). Cengage Learning.
2. Golub, G. H., & Van Loan, C. F. (2013). Matrix Computations (4th ed.). Johns Hopkins University Press.
3. Hale, J. K., & Verduyn Lunel, S. M. (1993). Introduction to Functional Differential Equations. Springer-Verlag.
4. Hairer, E., Nørsett, S. P., & Wanner, G. (1993). Solving Ordinary Differential Equations I: Nonstiff Problems (2nd ed.). Springer.
5. Hairer, E., & Wanner, G. (1996). Solving Ordinary Differential Equations II: Stiff and Differential-Algebraic Problems (2nd ed.). Springer.
