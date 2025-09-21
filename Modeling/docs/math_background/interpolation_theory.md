# Interpolation Theory

## 1. Polynomial Interpolation

### 1.1 Problem Statement

Let $\{(x_i, y_i)\}_{i=0}^n$ be a set of $n+1$ distinct data points where:

- $x_i \in [a,b] \subset \mathbb{R}$ are nodes with $x_i < x_{i+1}$
- $y_i = f(x_i)$ are function values for some $f: [a,b] \to \mathbb{R}$

The interpolation problem is to find a polynomial $p(x)$ of minimal degree such that:

$$
p(x_i) = y_i, \quad i = 0,\ldots,n
$$

It can be shown that such a polynomial has degree $\leq n$ and is unique.

### 1.2 Existence and Uniqueness

**Theorem 1.1** (Uniqueness of Interpolation Polynomial)
> There exists a unique polynomial $p(x)$ of degree $\leq n$ that interpolates the given data points.

*Proof*: Consider two interpolating polynomials $p(x)$ and $q(x)$. Their difference $d(x)=p(x)-q(x)$ has degree $\leq n$ and has $n+1$ zeros. By the Fundamental Theorem of Algebra, $d(x)$ must be identically zero.

## 2. Lagrange Interpolation

### 2.1 Direct Implementation Approach

The implementation uses the direct evaluation formula without storing basis polynomials:

For scalar-valued data:

```cpp
double result = 0.0;
for (size_t i = 0; i < n; ++i)
{
    double term = y_points[i];
    for (size_t j = 0; j < n; ++j)
    {
        if (i != j)
            term *= (x - x_points[j]) / (x_points[i] - x_points[j]);
    }
    result += term;
}
```

For vector-valued data, each component is interpolated independently using the same approach.

### 2.2 Lagrange Basis Polynomials

The mathematical foundation uses basis polynomials defined as:

$$
\ell_j(x) = \prod_{i=0,i\neq j}^n \frac{x-x_i}{x_j-x_i}
$$

**Properties**:

1. $\ell_j(x_k) = \delta_{jk}$ (Kronecker delta)
2. $\sum_{j=0}^n \ell_j(x) = 1$ (Partition of unity)

### 2.2 Lagrange Form

The interpolation polynomial in Lagrange form:

$$
p(x) = \sum_{j=0}^n y_j\ell_j(x)
$$

### 2.3 Error Analysis

Let $f \in C^{n+1}[a,b]$ be the function being interpolated. The interpolation error at any point $x \in [a,b]$ is given by:

$$f(x) - p(x) = \frac{f^{(n+1)}(\xi_x)}{(n+1)!}\prod_{i=0}^n (x-x_i)$$

where $\xi_x \in [a,b]$ depends on $x$.

#### 2.3.1 Error Bounds

The maximum interpolation error can be bounded by:

$$\|f-p\|_{\infty} \leq \frac{M_{n+1}}{4(n+1)}\left(\frac{b-a}{n}\right)^{n+1}$$

where $M_{n+1} = \max_{x\in[a,b]} |f^{(n+1)}(x)|$.

#### 2.3.2 Lebesgue Constant

The stability of polynomial interpolation is governed by the Lebesgue constant $\Lambda_n$:

$$\Lambda_n = \max_{x\in[a,b]} \sum_{j=0}^n |\ell_j(x)|$$

For equally spaced points:
$$\Lambda_n \sim \frac{2^{n+1}}{e\,n\,\log(n)}$$

For Chebyshev points:
$$\Lambda_n \leq \frac{2}{\pi}\log(n) + 1$$

## 3. Newton Interpolation

### 3.1 Implementation Strategy

The implementation uses two main components:

1. **Divided Difference Table Computation**:

   ```cpp
   // For scalar case:
   dVec coef = y_points;
   for (size_t j = 1; j < n; ++j)
   {
       for (size_t i = n - 1; i >= j; --i)
       {
           coef[i] = (coef[i] - coef[i - 1]) / 
                     (x_points[i] - x_points[i - j]);
       }
   }
   ```

2. **Horner's Scheme Evaluation**:

   ```cpp
   double result = coef.back();
   for (int i = static_cast<int>(coef.size()) - 2; i >= 0; --i)
   {
       result = result * (x - x_points[i]) + coef[i];
   }
   ```

### 3.2 Mathematical Foundation

The divided differences are defined recursively:

$$
\begin{align*}
f[x_i] &= f(x_i) \\[1em]
f[x_i,\ldots,x_{i+k}] &= \frac{f[x_{i+1},\ldots,x_{i+k}] - f[x_i,\ldots,x_{i+k-1}]}{x_{i+k}-x_i}
\end{align*}
$$

Leading to the Newton form:

$$
p(x) = f[x_0] + \sum_{k=1}^n f[x_0,\ldots,x_k]\prod_{i=0}^{k-1}(x-x_i)
$$

### 3.3 Vector-Valued Extension

For vector-valued data:

- Compute divided differences for each component
- Store coefficients in a matrix
- Apply Horner's scheme component-wise

### 3.3 Practical Implementation

#### 3.3.1 Divided Difference Table

Compute and store divided differences in a lower triangular matrix:

$$
\begin{bmatrix}
f[x_0]            & & & \\
f[x_0,x_1]        & f[x_1] & & \\
f[x_0,x_1,x_2]    & f[x_1,x_2] & f[x_2] & \\
\vdots            & \vdots      & \vdots & \ddots
\end{bmatrix}
$$

#### 3.3.2 Efficient Evaluation

Using Horner's scheme for evaluation:

```cpp
double evaluate(double x, const vector<double>& coef, const vector<double>& nodes)
{
    double result = coef[n];
    for(int i = n-1; i >= 0; --i)
    {
        result = result * (x - nodes[i]) + coef[i];
    }
    return result;
}
```

#### 3.3.3 Advantages

1. **Incremental Construction**: New points can be added without recalculating previous terms
2. **Error Estimation**: Easy computation of error bounds via divided differences
3. **Numerical Stability**: Better than Lagrange form for many cases
4. **Memory Efficient**: Only requires storage of coefficients and nodes

## 4. Barycentric Interpolation

### 4.1 Implementation Strategy

The implementation uses the second (true) barycentric form with precomputed weights:

1. **Weight Computation**:

   ```cpp
   dVec w(n, 1.0);
   for (size_t j = 0; j < n; ++j)
   {
       for (size_t k = 0; k < n; ++k)
       {
           if (j != k)
               w[j] /= (x_points[j] - x_points[k]);
       }
   }
   ```

2. **Evaluation Optimization**:
   - Direct return for exact node matches
   - Numerically stable division for interpolation
   - Separate handling for vector-valued data

### 4.2 Mathematical Foundation

The barycentric weights are:

$$
w_j = \frac{1}{\prod_{i=0,i\neq j}^n (x_j-x_i)}
$$

The true barycentric form used in implementation:

$$
p(x) = \frac{\displaystyle\sum_{j=0}^n \frac{w_j}{x-x_j}y_j}{\displaystyle\sum_{j=0}^n \frac{w_j}{x-x_j}}
$$

### 4.3 Implementation Details

For scalar-valued data:

```cpp
if (x == x_points[j])
    return y_points[j];  // Exact node match

double temp = w[j] / (x - x_points[j]);
num += temp * y_points[j];
denom += temp;
```

For vector-valued data:

- Compute weights once
- Apply formula component-wise
- Handle exact matches efficiently

### 4.3 Advantages

1. **Numerical Stability**: $\mathcal{O}(1)$ condition number
2. **Efficient Evaluation**: $\mathcal{O}(n)$ operations per point
3. **Node Independence**: Weights can be precomputed

## 5. Implementation Notes

### 5.1 Computational Complexity

| Operation | Lagrange | Newton | Barycentric |
|-----------|----------|--------|-------------|
| Setup     | $\mathcal{O}(n)$ | $\mathcal{O}(n^2)$ | $\mathcal{O}(n^2)$ |
| Evaluation| $\mathcal{O}(n^2)$ | $\mathcal{O}(n)$ | $\mathcal{O}(n)$ |
| Add Point | $\mathcal{O}(n^2)$ | $\mathcal{O}(n)$ | $\mathcal{O}(n^2)$ |

### 5.2 Error Bounds

For $f \in C^{n+1}[a,b]$:

$$
\|f-p\|_\infty \leq \frac{\|f^{(n+1)}\|_\infty}{(n+1)!}\|\omega\|_\infty
$$

where $\omega(x) = \prod_{i=0}^n (x-x_i)$

## 6. Node Selection and Computational Considerations

### 6.1 Optimal Node Distribution

#### 6.1.1 Chebyshev Nodes

The Chebyshev nodes on $[-1,1]$ are given by:

$$
x_j = \cos\left(\frac{2j+1}{2n+2}\pi\right), \quad j = 0,\ldots,n
$$

For general interval $[a,b]$:

$$
x_j = \frac{b+a}{2} + \frac{b-a}{2}\cos\left(\frac{2j+1}{2n+2}\pi\right)
$$

#### 6.1.2 Properties of Node Distributions

1. **Equally Spaced Nodes**:
   - Simple to implement
   - Poor convergence for high degrees
   - Runge phenomenon at endpoints

2. **Chebyshev Nodes**:
   - Near-optimal for minimizing maximum error
   - Excellent convergence properties
   - Minimal Runge phenomenon

### 6.2 Implementation Guidelines

1. **Choice of Method**:
   - Use Newton for incremental construction
   - Use Barycentric for repeated evaluation
   - Use Lagrange for theoretical analysis

2. **Error Control**:
   - Monitor condition number via Lebesgue constant
   - Use adaptive degree selection
   - Implement error estimators

3. **Numerical Considerations**:
   - Scale input domain to $[-1,1]$ for stability
   - Use Chebyshev nodes for high degrees
   - Implement careful overflow protection

### 6.3 Performance Optimization

1. **Precomputation Strategies**:
   - Cache barycentric weights
   - Store divided difference tables
   - Precompute basis functions when possible

2. **Vectorization Opportunities**:
   - Parallelize node evaluations
   - Use SIMD instructions for basis functions
   - Batch interpolation for multiple points

## References

1. Berrut, J. P., & Trefethen, L. N. (2004). Barycentric Lagrange Interpolation.
2. Burden, R. L., & Faires, J. D. (2015). Numerical Analysis.
3. Trefethen, L. N. (2013). Approximation Theory and Approximation Practice.
4. Higham, N. J. (2004). The numerical stability of barycentric Lagrange interpolation.
5. Rack, H. J., & Reimer, M. (1982). The numerical stability of evaluation schemes for polynomials based on the Lagrange interpolation form.
