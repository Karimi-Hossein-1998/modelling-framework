# Linear Algebra Module

The `linalg` module provides fundamental linear algebra operations, including matrix decompositions, basic matrix and vector arithmetic, and solvers for linear systems.

## Basic Operations (`basics/matrix_operations.hpp`)

This section covers common operations on vectors and matrices.

### `dot_product`

Computes the dot product of two vectors.

```cpp
inline double dot_product(const dVec& a, const dVec& b);
```

- `a`, `b`: Input vectors.

**Returns:** The scalar dot product.

**Example:**

```cpp
dVec v1 = {1.0, 2.0, 3.0};
dVec v2 = {4.0, 5.0, 6.0};
double dp = dot_product(v1, v2); // dp will be 32.0
```

### `cross_product` (3D)

Computes the cross product of two 3D vectors.

```cpp
inline triple cross_product(const triple& a, const triple& b);
```

- `a`, `b`: Input 3D vectors (represented as `triple`).

**Returns:** A `triple` representing the cross product vector.

**Example:**

```cpp
triple t1 = {1.0, 0.0, 0.0};
triple t2 = {0.0, 1.0, 0.0};
triple cp = cross_product(t1, t2); // cp will be {0.0, 0.0, 1.0}
```

### `matrix_multiply` (Matrix-Vector)

Computes the product of a matrix and a vector ($A \mathbf{x}$).

```cpp
inline dVec matrix_multiply(const Matrix& A, const dVec& x);
```

- `A`: Input matrix.
- `x`: Input vector.

**Returns:** A `dVec` representing the resulting vector.

**Example:**

```cpp
Matrix mat = {{1, 2}, {3, 4}};
dVec vec = {1, 2};
dVec result = matrix_multiply(mat, vec); // result will be {5, 11}
```

### `matrix_multiply` (Matrix-Matrix)

Computes the product of two matrices ($AB$).

```cpp
inline Matrix matrix_multiply(const Matrix& A, const Matrix& B);
```

- `A`, `B`: Input matrices.

**Returns:** A `Matrix` representing the resulting product matrix.

**Example:**

```cpp
Matrix mat1 = {{1, 2}, {3, 4}};
Matrix mat2 = {{5, 6}, {7, 8}};
Matrix result = matrix_multiply(mat1, mat2); // result will be {{19, 22}, {43, 50}}
```

### `elementwise_multiply` (Matrix)

Computes the element-wise product of two matrices (Hadamard product).

```cpp
inline Matrix elementwise_multiply(const Matrix& A, const Matrix& B);
```

- `A`, `B`: Input matrices of the same dimensions.

**Returns:** A `Matrix` with elements $C_{ij} = A_{ij} \cdot B_{ij}$.

### `elementwise_multiply` (Vector)

Computes the element-wise product of two vectors.

```cpp
inline dVec elementwise_multiply(const dVec& a, const dVec& b);
```

- `a`, `b`: Input vectors of the same size.

**Returns:** A `dVec` with elements $c_i = a_i \cdot b_i$.

### `add` (Matrix)

Computes the sum of two matrices.

```cpp
inline Matrix add(const Matrix& A, const Matrix& B);
```

- `A`, `B`: Input matrices of the same dimensions.

**Returns:** A `Matrix` with elements $C_{ij} = A_{ij} + B_{ij}$.

### `add` (Vector)

Computes the sum of two vectors.

```cpp
inline dVec add(const dVec& a, const dVec& b);
```

- `a`, `b`: Input vectors of the same size.

**Returns:** A `dVec` with elements $c_i = a_i + b_i$.

### `subtract` (Matrix)

Computes the difference of two matrices.

```cpp
inline Matrix subtract(const Matrix& A, const Matrix& B);
```

- `A`, `B`: Input matrices of the same dimensions.

**Returns:** A `Matrix` with elements $C_{ij} = A_{ij} - B_{ij}$.

### `subtract` (Vector)

Computes the difference of two vectors.

```cpp
inline dVec subtract(const dVec& a, const dVec& b);
```

- `a`, `b`: Input vectors of the same size.

**Returns:** A `dVec` with elements $c_i = a_i - b_i$.

### `scalar_multiply` (Matrix)

Computes the scalar product of a matrix ($sA$).

```cpp
inline Matrix scalar_multiply(double s, const Matrix& A);
```

- `s`: Scalar value.
- `A`: Input matrix.

**Returns:** A `Matrix` with elements $C_{ij} = s \cdot A_{ij}$.

### `scalar_multiply` (Vector)

Computes the scalar product of a vector ($s\mathbf{v}$). 

```cpp
inline dVec scalar_multiply(double s, const dVec& v);
```

- `s`: Scalar value.
- `v`: Input vector.

**Returns:** A `dVec` with elements $c_i = s \cdot v_i$.

### `vector_transpose_multiply`

Computes the product of a transposed vector and a matrix ($\mathbf{x}^T A$).

```cpp
inline dVec vector_transpose_multiply(const dVec& x, const Matrix& A);
```

- `x`: Input vector.
- `A`: Input matrix.

**Returns:** A `dVec` representing the resulting row vector.

### `outer_product`

Computes the outer product of two vectors ($\mathbf{x} \otimes \mathbf{y}$ or $|\mathbf{x}\rangle\langle\mathbf{y}|$).

```cpp
inline Matrix outer_product(const dVec& x, const dVec& y);
```

- `x`, `y`: Input vectors.

**Returns:** A `Matrix` where $C_{ij} = x_i \cdot y_j$.

### `vector_transpose`

Returns the transpose of a vector as a 1xN matrix.

```cpp
inline Matrix vector_transpose(const dVec& x);
```

- `x`: Input vector.

**Returns:** A `Matrix` with one row, containing the elements of `x`.

### `matrix_outer_product` (Kronecker product)

Computes the outer product (Kronecker product) of two matrices.

```cpp
inline Matrix matrix_outer_product(const Matrix& A, const Matrix& B);
```

- `A`, `B`: Input matrices.

**Returns:** The Kronecker product matrix.

### `cross_product_2d`

Computes the 2D cross product (scalar result) of two 2D vectors.

```cpp
inline double cross_product_2d(const dVec& a, const dVec& b);
```

- `a`, `b`: Input 2D vectors.

**Returns:** A scalar value.

### `perpendicular_vector`

Returns a vector perpendicular to two input vectors in n-dimensional space. For $n > 3$, this is not unique and returns one possible perpendicular vector.

```cpp
inline dVec perpendicular_vector(const dVec& a, const dVec& b);
```

- `a`, `b`: Input vectors of the same size (at least 2).

**Returns:** A `dVec` representing a perpendicular vector.

### `generalized_cross_product`

Computes the generalized cross product for n-dimensional vectors. Requires $n-1$ input vectors of size $n$.

```cpp
inline dVec generalized_cross_product(const Matrix& vectors);
```

- `vectors`: A `Matrix` where each row is an input vector. There should be $n-1$ rows, and each row should have $n$ columns.

**Returns:** A `dVec` representing the generalized cross product vector.

## Matrix Properties (`basics/matrix_properties.hpp`)

This section provides functions to query properties of matrices.

### `trace`

Returns the trace of a square matrix (sum of diagonal elements).

```cpp
inline double trace(const Matrix& A);
```

- `A`: Input square matrix.

**Returns:** The trace of the matrix.

### `determinant`

Computes the determinant of a square matrix. Uses optimized methods (direct for small matrices, LU decomposition for larger ones).

```cpp
inline double determinant(const Matrix& A);
```

- `A`: Input square matrix.

**Returns:** The determinant of the matrix.

### `minor`

Computes the minor of a matrix at a specified row and column.

```cpp
inline double minor(const Matrix& A, size_t row, size_t col);
```

- `A`: Input matrix.
- `row`, `col`: The row and column for which to compute the minor.

**Returns:** The minor value.

### `transpose`

Returns the transpose of a matrix.

```cpp
inline Matrix transpose(const Matrix& A);
```

- `A`: Input matrix.

**Returns:** The transposed matrix.

### `cofactor_matrix`

Computes the cofactor matrix of a square matrix.

```cpp
inline Matrix cofactor_matrix(const Matrix& A);
```

- `A`: Input square matrix.

**Returns:** The cofactor matrix.

### `adjugate`

Returns the adjugate (classical adjoint) of a square matrix.

```cpp
inline Matrix adjugate(const Matrix& A);
```

- `A`: Input square matrix.

**Returns:** The adjugate matrix.

### `inverse`

Computes the inverse of a square matrix. Uses optimized methods (adjugate for small matrices, LU decomposition for larger ones).

```cpp
inline Matrix inverse(const Matrix& A);
```

- `A`: Input square matrix.

**Returns:** The inverse matrix.

### `is_symmetric`

Checks if a matrix is symmetric.

```cpp
inline bool is_symmetric(const Matrix& A);
```

- `A`: Input matrix.

**Returns:** `true` if the matrix is symmetric, `false` otherwise.

### `is_identity`

Checks if a matrix is an identity matrix.

```cpp
inline bool is_identity(const Matrix& A);
```

- `A`: Input matrix.

**Returns:** `true` if the matrix is an identity matrix, `false` otherwise.

### `is_diagonal`

Checks if a matrix is a diagonal matrix.

```cpp
inline bool is_diagonal(const Matrix& A);
```

- `A`: Input matrix.

**Returns:** `true` if the matrix is diagonal, `false` otherwise.

### `is_upper_triangular`

Checks if a matrix is upper triangular.

```cpp
inline bool is_upper_triangular(const Matrix& A);
```

- `A`: Input matrix.

**Returns:** `true` if the matrix is upper triangular, `false` otherwise.

### `is_lower_triangular`

Checks if a matrix is lower triangular.

```cpp
inline bool is_lower_triangular(const Matrix& A);
```

- `A`: Input matrix.

**Returns:** `true` if the matrix is lower triangular, `false` otherwise.

### `frobenius_norm`

Computes the Frobenius norm of a matrix.

```cpp
inline double frobenius_norm(const Matrix& A);
```

- `A`: Input matrix.

**Returns:** The Frobenius norm.

### `max_norm` (Infinity Norm)

Computes the maximum absolute row sum norm (infinity norm) of a matrix.

```cpp
inline double max_norm(const Matrix& A);
```

- `A`: Input matrix.

**Returns:** The infinity norm.

### `one_norm`

Computes the maximum absolute column sum norm (1-norm) of a matrix.

```cpp
inline double one_norm(const Matrix& A);
```

- `A`: Input matrix.

**Returns:** The 1-norm.

### `rank`

Computes the rank of a matrix using row echelon form (suitable for small matrices).

```cpp
inline size_t rank(const Matrix& A);
```

- `A`: Input matrix.

**Returns:** The rank of the matrix.

## Eigenvalue Decomposition (`basics/eigen.hpp`)

This section provides functions for computing eigenvalues and eigenvectors.

### `power_method`

Finds the dominant eigenvalue and its corresponding eigenvector of a square matrix using the power method.

```cpp
inline std::pair<double, dVec> power_method(
    const Matrix& A,
    size_t max_iter = 100,
    double tol = 1e-10
);
```

- `A`: Input square matrix.
- `max_iter`: Maximum number of iterations (default: 100).
- `tol`: Tolerance for convergence (default: `1e-10`).

**Returns:** A `std::pair` containing the dominant eigenvalue and its corresponding eigenvector.

**Throws:** `std::runtime_error` if the method does not converge.

### `to_hessenberg`

Reduces a matrix to upper Hessenberg form using Householder reflections. This is a preliminary step for more robust eigenvalue algorithms like QR iteration.

```cpp
inline Matrix to_hessenberg(Matrix A);
```

- `A`: Input matrix.

**Returns:** The matrix in upper Hessenberg form.

### `eigen_qr`

Computes all eigenvalues and eigenvectors of a square matrix using the QR iteration with shifts. This is a more general and robust method than the power method.

```cpp
inline std::pair<dVec, Matrix> eigen_qr(
    Matrix A,
    size_t max_iter = 100,
    double tol = 1e-10
);
```

- `A`: Input square matrix.
- `max_iter`: Maximum number of iterations (default: 100).
- `tol`: Tolerance for convergence (default: `1e-10`).

**Returns:** A `std::pair` where the first element is a `dVec` of eigenvalues, and the second element is a `Matrix` where columns are the corresponding eigenvectors.

**Throws:** `std::runtime_error` if the QR iteration does not converge.

### `eigen_decomposition`

Returns the Eigen decomposition (eigenvalues, eigenvectors) for a square matrix. This is a wrapper around `eigen_qr`.

```cpp
inline std::pair<dVec, Matrix> eigen_decomposition(const Matrix& A);
```

- `A`: Input square matrix.

**Returns:** A `std::pair` containing eigenvalues and eigenvectors.

## Matrix Decompositions (`basics/decomposers.hpp`)

This section provides functions for various matrix decompositions.

### `lu_decomposition`

Performs LU decomposition with partial pivoting ($PA = LU$).

```cpp
inline std::tuple<Matrix, Matrix, wVec> lu_decomposition(Matrix A);
```

- `A`: Input square matrix.

**Returns:** A `std::tuple` containing the lower triangular matrix `L`, the upper triangular matrix `U`, and the permutation vector `P`.

**Throws:** `std::invalid_argument` if the matrix is not square, or `std::runtime_error` if the matrix is singular or nearly singular.

### `qr_decomposition`

Performs QR decomposition using modified Gram-Schmidt for an $m \times n$ matrix ($A = QR$).

```cpp
inline std::pair<Matrix, Matrix> qr_decomposition(const Matrix& A);
```

- `A`: Input matrix.

**Returns:** A `std::pair` containing the orthogonal matrix `Q` and the upper triangular matrix `R`.

**Throws:** `std::runtime_error` if the matrix is rank deficient.

### `cholesky_decomposition`

Performs Cholesky decomposition ($A = LL^T$) for a symmetric positive definite matrix.

```cpp
inline Matrix cholesky_decomposition(const Matrix& A);
```

- `A`: Input symmetric positive definite matrix.

**Returns:** The lower triangular matrix `L`.

**Throws:** `std::invalid_argument` if the matrix is not square, or `std::runtime_error` if the matrix is not positive definite or has a zero pivot.

### `ldlt_decomposition`

Performs LDLT decomposition ($A = LDL^T$) for a symmetric matrix.

```cpp
inline std::pair<Matrix, dVec> ldlt_decomposition(const Matrix& A);
```

- `A`: Input symmetric matrix.

**Returns:** A `std::pair` containing the lower triangular matrix `L` and the diagonal vector `D`.

**Throws:** `std::invalid_argument` if the matrix is not square, or `std::runtime_error` if a zero pivot is encountered.

### `svd`

Performs Singular Value Decomposition ($A = USV^T$) for a square matrix.

```cpp
inline std::tuple<Matrix, dVec, Matrix> svd(const Matrix& A);
```

- `A`: Input square matrix.

**Returns:** A `std::tuple` containing the orthogonal matrix `U`, a `dVec` of singular values `S`, and the orthogonal matrix `V`.

## Linear System Solvers (`linsyssolvers/`)

This section provides functions for solving linear systems of equations of the form $A\mathbf{x} = \mathbf{b}$.

### `banded_solve`

Solves a linear system $A\mathbf{x} = \mathbf{b}$ where $A$ is a banded matrix, using compact banded storage and Gaussian elimination.

```cpp
inline dVec banded_solve(Matrix ab, int kl, int ku, dVec b);
```

- `ab`: Banded matrix storage (2*kl+ku+1) x n, where `ab[kl+ku+i-j][j] == A[i][j]`.
- `kl`: Lower bandwidth.
- `ku`: Upper bandwidth.
- `b`: Right-hand side vector.

**Returns:** The solution vector $\mathbf{x}$.

**Throws:** `std::invalid_argument` for invalid matrix/vector sizes, or `std::runtime_error` for zero pivot.

### `gaussian_elimination`

Solves a linear system $A\mathbf{x} = \mathbf{b}$ using Gaussian elimination with partial pivoting.

```cpp
inline dVec gaussian_elimination(Matrix A, dVec b);
```

- `A`: Input square matrix.
- `b`: Right-hand side vector.

**Returns:** The solution vector $\mathbf{x}$.

**Throws:** `std::invalid_argument` for dimension mismatch, or `std::runtime_error` if the matrix is singular.

### `ldlt_solve`

Solves a linear system $A\mathbf{x} = \mathbf{b}$ for a symmetric matrix $A$ using LDLT decomposition.

```cpp
inline dVec ldlt_solve(Matrix A, dVec b);
```

- `A`: Input symmetric square matrix.
- `b`: Right-hand side vector.

**Returns:** The solution vector $\mathbf{x}$.

**Throws:** `std::invalid_argument` for dimension mismatch.

### `llt_solve`

Solves a linear system $A\mathbf{x} = \mathbf{b}$ for a symmetric positive definite matrix $A$ using Cholesky (LLT) decomposition.

```cpp
inline dVec llt_solve(Matrix A, dVec b);
```

- `A`: Input symmetric positive definite square matrix.
- `b`: Right-hand side vector.

**Returns:** The solution vector $\mathbf{x}$.

**Throws:** `std::invalid_argument` for dimension mismatch, or `std::runtime_error` if the matrix is not positive definite.

### `lu_solve`

Solves a linear system $A\mathbf{x} = \mathbf{b}$ using LU factorization with partial pivoting.

```cpp
inline dVec lu_solve(Matrix A, dVec b);
```

- `A`: Input square matrix.
- `b`: Right-hand side vector.

**Returns:** The solution vector $\mathbf{x}$.

**Throws:** `std::invalid_argument` for dimension mismatch, or `std::runtime_error` if the matrix is singular.

### `qr_solve`

Solves a linear system $A\mathbf{x} = \mathbf{b}$ using QR decomposition (Gram-Schmidt).

```cpp
inline dVec qr_solve(Matrix A, dVec b);
```

- `A`: Input square matrix.
- `b`: Right-hand side vector.

**Returns:** The solution vector $\mathbf{x}$.

**Throws:** `std::invalid_argument` for dimension mismatch, or `std::runtime_error` if the matrix is rank deficient.

### `svd_solve`

Solves a linear system $A\mathbf{x} = \mathbf{b}$ using Singular Value Decomposition (SVD). This implementation uses a basic Jacobi SVD and is primarily for demonstration purposes.

```cpp
inline dVec svd_solve(Matrix A, dVec b);
```

- `A`: Input square matrix.
- `b`: Right-hand side vector.

**Returns:** The solution vector $\mathbf{x}$.

**Throws:** `std::invalid_argument` for dimension mismatch.

### `tridiagonal_solve`

Solves a linear system $A\mathbf{x} = \mathbf{d}$ where $A$ is a tridiagonal matrix, using the Thomas algorithm.

```cpp
inline dVec tridiagonal_solve(
    const dVec& a, // Subdiagonal (n-1 elements)
    const dVec& b, // Diagonal (n elements)
    const dVec& c, // Superdiagonal (n-1 elements)
    const dVec& d  // Right-hand side (n elements)
);
```

- `a`: Subdiagonal vector.
- `b`: Diagonal vector.
- `c`: Superdiagonal vector.
- `d`: Right-hand side vector.

**Returns:** The solution vector $\mathbf{x}$.

**Throws:** `std::invalid_argument` for invalid vector sizes.
