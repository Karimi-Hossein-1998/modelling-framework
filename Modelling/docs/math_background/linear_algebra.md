# Linear Algebra Mathematical Background

## 1. Vector Spaces and Inner Products

### 1.1 Vector Spaces

A vector space $V$ over field $\mathbb{F}$ satisfies eight axioms, divided into:

Vector Addition Axioms:

- Closure: $\mathbf{x}, \mathbf{y} \in V \implies \mathbf{x} + \mathbf{y} \in V$
- Commutativity: $\mathbf{x} + \mathbf{y} = \mathbf{y} + \mathbf{x}$
- Associativity: $(\mathbf{x} + \mathbf{y}) + \mathbf{z} = \mathbf{x} + (\mathbf{y} + \mathbf{z})$
- Zero vector: $\exists \mathbf{0} \in V: \mathbf{x} + \mathbf{0} = \mathbf{x}$
- Additive inverse: $\forall \mathbf{x} \in V, \exists -\mathbf{x}: \mathbf{x} + (-\mathbf{x}) = \mathbf{0}$

Scalar Multiplication Axioms:

- Closure: $\alpha \in \mathbb{F}, \mathbf{x} \in V \implies \alpha\mathbf{x} \in V$
- Distributivity: $\alpha(\mathbf{x} + \mathbf{y}) = \alpha\mathbf{x} + \alpha\mathbf{y}$
- Scalar distributivity: $(\alpha + \beta)\mathbf{x} = \alpha\mathbf{x} + \beta\mathbf{x}$

### 1.2 Inner Product Spaces

An inner product $\langle \cdot,\cdot \rangle$ on $V$ satisfies:

1. Conjugate symmetry:
   $$\langle \mathbf{x},\mathbf{y} \rangle = \overline{\langle \mathbf{y},\mathbf{x} \rangle}$$

2. Linearity in first argument:
   $$\langle \alpha\mathbf{x}+\mathbf{y},\mathbf{z} \rangle = \alpha\langle \mathbf{x},\mathbf{z} \rangle + \langle \mathbf{y},\mathbf{z} \rangle$$

3. Positive definiteness:
   $$\langle \mathbf{x},\mathbf{x} \rangle > 0 \text{ if } \mathbf{x} \neq \mathbf{0}$$

## 2. Matrix Theory

### 2.1 Matrix Properties

#### 2.1.1 Basic Definitions

For $\mathbf{A} \in \mathbb{R}^{m \times n}$:

1. Transpose:
   $$(\mathbf{A}^T)_{ij} = a_{ji}$$

2. Trace (for square matrices):
   $$\text{tr}(\mathbf{A}) = \sum_{i=1}^n a_{ii}$$

3. Determinant:
   $$\det(\mathbf{A}) = \sum_{\sigma \in S_n} \text{sgn}(\sigma)\prod_{i=1}^n a_{i\sigma(i)}$$

#### 2.1.2 Special Matrices

1. Symmetric: $\mathbf{A} = \mathbf{A}^T$
2. Orthogonal: $\mathbf{A}\mathbf{A}^T = \mathbf{A}^T\mathbf{A} = \mathbf{I}$
3. Positive Definite: $\mathbf{x}^T\mathbf{A}\mathbf{x} > 0$ for all $\mathbf{x} \neq \mathbf{0}$

### 2.2 Matrix Decompositions

#### 2.2.1 LU Decomposition

For a nonsingular matrix $\mathbf{A}$:
$$\mathbf{A} = \mathbf{L}\mathbf{U}$$

where:

- $\mathbf{L}$ is lower triangular
- $\mathbf{U}$ is upper triangular

Algorithm (Doolittle):

```text
for k = 1 to n
    for i = k to n
        for j = 1 to k-1
            a[i,k] = a[i,k] - a[i,j]*a[j,k]
    for i = k+1 to n
        a[i,k] = a[i,k]/a[k,k]
        for j = k+1 to n
            a[i,j] = a[i,j] - a[i,k]*a[k,j]
```

#### 2.2.2 Cholesky Decomposition

For symmetric positive definite $\mathbf{A}$:
$$\mathbf{A} = \mathbf{L}\mathbf{L}^T$$

Algorithm:

```text
for j = 1 to n
    for i = j to n
        sum = a[i,j]
        for k = 1 to j-1
            sum = sum - a[i,k]*a[j,k]
        if i == j
            a[j,j] = sqrt(sum)
        else
            a[i,j] = sum/a[j,j]
```

#### 2.2.3 QR Decomposition

$$\mathbf{A} = \mathbf{Q}\mathbf{R}$$

where:

- $\mathbf{Q}$ is orthogonal
- $\mathbf{R}$ is upper triangular

Modified Gram-Schmidt Process:

```text
for j = 1 to n
    for i = 1 to j-1
        r[i,j] = dot_product(q[:,i], a[:,j])
        a[:,j] = a[:,j] - r[i,j]*q[:,i]
    r[j,j] = norm(a[:,j])
    q[:,j] = a[:,j]/r[j,j]
```

## 3. Eigenvalue Theory

### 3.1 Basic Concepts

For square matrix $\mathbf{A}$:
$$\mathbf{A}\mathbf{x} = \lambda\mathbf{x}$$

where:

- $\lambda$ is an eigenvalue
- $\mathbf{x}$ is the corresponding eigenvector

Characteristic equation:
$$\det(\mathbf{A} - \lambda\mathbf{I}) = 0$$

### 3.2 Power Method

For finding dominant eigenvalue:

1. Start with unit vector $\mathbf{x}_0$
2. Iterate: $\mathbf{x}_{k+1} = \mathbf{A}\mathbf{x}_k/\|\mathbf{A}\mathbf{x}_k\|$
3. Eigenvalue estimate: $\lambda_k = \mathbf{x}_k^T\mathbf{A}\mathbf{x}_k$

Convergence rate:
$$|\lambda_1| > |\lambda_2| \implies \|\mathbf{e}_k\| \sim \left|\frac{\lambda_2}{\lambda_1}\right|^k$$

### 3.3 QR Algorithm

Basic algorithm:

1. Start with $\mathbf{A}_0 = \mathbf{A}$
2. For k = 0,1,...:
   - Factor $\mathbf{A}_k = \mathbf{Q}_k\mathbf{R}_k$
   - Set $\mathbf{A}_{k+1} = \mathbf{R}_k\mathbf{Q}_k$

With shifts:
$$\mathbf{A}_k - \mu_k\mathbf{I} = \mathbf{Q}_k\mathbf{R}_k$$
$$\mathbf{A}_{k+1} = \mathbf{R}_k\mathbf{Q}_k + \mu_k\mathbf{I}$$

## 4. Linear Systems

### 4.1 Direct Methods

#### 4.1.1 Gaussian Elimination

For system $\mathbf{A}\mathbf{x} = \mathbf{b}$:

1. Forward elimination:
   - Create upper triangular system
   - Update RHS accordingly

2. Back substitution:
   $$x_i = \frac{1}{a_{ii}}\left(b_i - \sum_{j=i+1}^n a_{ij}x_j\right)$$

#### 4.1.2 Band Systems

For bandwidth $p$:
$$a_{ij} = 0 \text{ if } |i-j| > p$$

Storage requirement: $O(np)$ instead of $O(n^2)$

### 4.2 Iterative Methods

#### 4.2.1 Jacobi Method

$$x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j\neq i} a_{ij}x_j^{(k)}\right)$$

Convergence: Requires diagonal dominance

#### 4.2.2 Gauss-Seidel Method

$$x_i^{(k+1)} = \frac{1}{a_{ii}}\left(b_i - \sum_{j<i} a_{ij}x_j^{(k+1)} - \sum_{j>i} a_{ij}x_j^{(k)}\right)$$

## References

1. Golub, G. H., & Van Loan, C. F. (2013). Matrix Computations (4th ed.).
2. Trefethen, L. N., & Bau III, D. (1997). Numerical Linear Algebra.
3. Horn, R. A., & Johnson, C. R. (2012). Matrix Analysis.
