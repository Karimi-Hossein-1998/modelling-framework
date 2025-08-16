# Interpolators Module

The `interpolators` module provides functions for constructing interpolation polynomials, allowing you to estimate values between known data points. It includes implementations for both Lagrange and Newton interpolation, supporting scalar and vector-valued data.

## Lagrange Interpolation

Lagrange interpolation constructs a polynomial that passes through a given set of data points. It is defined as:

$$P(x) = \sum_{j=0}^{k} y_j L_j(x)$$

where $L_j(x)$ are the Lagrange basis polynomials.

### `lagrange_interpolator` (Scalar-valued)

Returns a function (lambda) that performs Lagrange interpolation for scalar `y` values given `x_points` and `y_points`.

```cpp
inline func lagrange_interpolator(const dVec& x_points, const dVec& y_points);
```

- `x_points`: A `dVec` of x-coordinates of the data points.
- `y_points`: A `dVec` of y-coordinates (scalar values) of the data points.

**Returns:** A `func` (lambda `double(double)`) that takes an `x` value and returns the interpolated `y` value.

**Example:**

```cpp
dVec x_data = {0.0, 1.0, 2.0};
dVec y_data = {0.0, 1.0, 4.0}; // y = x^2
func interpolate_scalar = lagrange_interpolator(x_data, y_data);
double interpolated_value = interpolate_scalar(0.5); // Should be close to 0.25
```

### `lagrange_interpolator` (Vector-valued)

Returns a function (lambda) that performs Lagrange interpolation for vector-valued `y` data. Each component of the vector is interpolated independently.

```cpp
inline Func lagrange_interpolator(const dVec& x_points, const Matrix& y_points);
```

- `x_points`: A `dVec` of x-coordinates of the data points.
- `y_points`: A `Matrix` where each row corresponds to the vector-valued `y` data at the respective `x_point`.

**Returns:** A `Func` (lambda `dVec(double)`) that takes an `x` value and returns the interpolated `dVec`.

**Example:**

```cpp
dVec x_data_vec = {0.0, 1.0};
Matrix y_data_vec = {{0.0, 0.0}, {1.0, 2.0}}; // y_points[i][k] is k-th component at x_points[i]
Func interpolate_vector = lagrange_interpolator(x_data_vec, y_data_vec);
dVec interpolated_vec = interpolate_vector(0.5); // Should be close to {0.5, 1.0}
```

## Barycentric Lagrange Interpolation

Barycentric Lagrange interpolation is a more numerically stable and efficient form of Lagrange interpolation, especially when evaluating the polynomial at many points.

### `barycentric_lagrange_interpolator` (Scalar-valued)

Returns a function (lambda) that performs barycentric Lagrange interpolation for scalar `y` values.

```cpp
inline func barycentric_lagrange_interpolator(const dVec& x_points, const dVec& y_points);
```

- `x_points`: A `dVec` of x-coordinates of the data points.
- `y_points`: A `dVec` of y-coordinates (scalar values) of the data points.

**Returns:** A `func` (lambda `double(double)`) that takes an `x` value and returns the interpolated `y` value.

**Example:**

```cpp
dVec x_bary = {0.0, 1.0, 2.0};
dVec y_bary = {0.0, 1.0, 4.0};
func interpolate_bary_scalar = barycentric_lagrange_interpolator(x_bary, y_bary);
double interpolated_bary_value = interpolate_bary_scalar(0.5);
```

### `barycentric_lagrange_interpolator` (Vector-valued)

Returns a function (lambda) that performs barycentric Lagrange interpolation for vector-valued `y` data.

```cpp
inline Func barycentric_lagrange_interpolator(const dVec& x_points, const Matrix& y_points);
```

- `x_points`: A `dVec` of x-coordinates of the data points.
- `y_points`: A `Matrix` where each row corresponds to the vector-valued `y` data at the respective `x_point`.

**Returns:** A `Func` (lambda `dVec(double)`) that takes an `x` value and returns the interpolated `dVec`.

**Example:**

```cpp
dVec x_bary_vec = {0.0, 1.0};
Matrix y_bary_vec = {{0.0, 0.0}, {1.0, 2.0}};
Func interpolate_bary_vector = barycentric_lagrange_interpolator(x_bary_vec, y_bary_vec);
dVec interpolated_bary_vec = interpolate_bary_vector(0.5);
```

## Newton Interpolation

Newton interpolation uses divided differences to construct the interpolating polynomial. It is particularly useful when adding new data points, as the existing calculations can be reused.

### `newton_interpolator` (Scalar-valued)

Returns a function (lambda) that performs Newton interpolation for scalar `y` values.

```cpp
inline func newton_interpolator(const dVec& x_points, const dVec& y_points);
```

- `x_points`: A `dVec` of x-coordinates of the data points.
- `y_points`: A `dVec` of y-coordinates (scalar values) of the data points.

**Returns:** A `func` (lambda `double(double)`) that takes an `x` value and returns the interpolated `y` value.

**Example:**

```cpp
dVec x_newton = {0.0, 1.0, 2.0};
dVec y_newton = {0.0, 1.0, 4.0};
func interpolate_newton_scalar = newton_interpolator(x_newton, y_newton);
double interpolated_newton_value = interpolate_newton_scalar(0.5);
```

### `newton_interpolator` (Vector-valued)

Returns a function (lambda) that performs Newton interpolation for vector-valued `y` data.

```cpp
inline Func newton_interpolator(const dVec& x_points, const Matrix& y_points);
```

- `x_points`: A `dVec` of x-coordinates of the data points.
- `y_points`: A `Matrix` where each row corresponds to the vector-valued `y` data at the respective `x_point`.

**Returns:** A `Func` (lambda `dVec(double)`) that takes an `x` value and returns the interpolated `dVec`.

**Example:**

```cpp
dVec x_newton_vec = {0.0, 1.0};
Matrix y_newton_vec = {{0.0, 0.0}, {1.0, 2.0}};
Func interpolate_newton_vector = newton_interpolator(x_newton_vec, y_newton_vec);
dVec interpolated_newton_vec = interpolate_newton_vector(0.5);
```
