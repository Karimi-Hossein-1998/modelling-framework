# Initializers Module

The `initializers` module provides functions for generating initial conditions for simulations, particularly for setting up phases or frequencies based on various statistical distributions.

## Functions

### `random_uniform`

Generates a vector of `N` double values uniformly distributed within a specified range `[min_val, max_val)`.

```cpp
inline dVec random_uniform(
    size_t   N,
    double   min_val,
    double   max_val,
    unsigned seed
);
```

- `N`: The number of values to generate.
- `min_val`: The minimum value of the uniform distribution (inclusive).
- `max_val`: The maximum value of the uniform distribution (exclusive).
- `seed`: The seed for the random number generator.

**Example:**

```cpp
dVec phases = random_uniform(100, 0.0, 2.0 * PI, 123);
```

### `random_normal`

Generates a vector of `N` double values following a normal (Gaussian) distribution with a specified `mean` and `stddev`.

```cpp
inline dVec random_normal(
    size_t      N,
    double      mean,
    double      stddev,
    unsigned    seed
);
```

- `N`: The number of values to generate.
- `mean`: The mean of the normal distribution.
- `stddev`: The standard deviation of the normal distribution.
- `seed`: The seed for the random number generator.

**Example:**

```cpp
dVec frequencies = random_normal(50, 10.0, 2.0, 456);
```

### `random_cauchy`

Generates a vector of `N` double values following a Cauchy (Lorentzian) distribution with a specified `location` and `scale`.

```cpp
inline dVec random_cauchy(
    size_t      N,
    double      location,
    double      scale,
    unsigned    seed
);
```

- `N`: The number of values to generate.
- `location`: The location parameter of the Cauchy distribution.
- `scale`: The scale parameter of the Cauchy distribution.
- `seed`: The seed for the random number generator.

**Example:**

```cpp
dVec values = random_cauchy(75, 0.0, 1.0, 789);
```

### `random_exponential`

Generates a vector of `N` double values following an exponential distribution with a specified `lambda` (rate parameter).

```cpp
inline dVec random_exponential(
    size_t      N,
    double      lambda, // rate parameter
    unsigned    seed
);
```

- `N`: The number of values to generate.
- `lambda`: The rate parameter of the exponential distribution.
- `seed`: The seed for the random number generator.

**Example:**

```cpp
dVec times = random_exponential(60, 0.5, 101);
```

### `random_circle`

Generates a vector of `N` double values uniformly distributed on the unit circle `[0, 2*pi)`.

```cpp
inline dVec random_circle(
    size_t      N,
    unsigned    seed
);
```

- `N`: The number of values to generate.
- `seed`: The seed for the random number generator.

**Example:**

```cpp
dVec phases_circle = random_circle(120, 222);
```

### `splay`

Generates `N` equidistant values (splay phases) around the circle `[0, 2*pi)`.

```cpp
inline dVec splay(
    size_t N
);
```

- `N`: The number of splay phases to generate.

**Example:**

```cpp
dVec splay_phases = splay(8);
```

### `splay_perturbed`

Generates `N` splay phases with a random perturbation applied to each phase. The perturbation is uniformly distributed in `[-amplitude, amplitude]`.

```cpp
inline dVec splay_perturbed(
    size_t   N,
    double   amplitude,
    unsigned seed
);
```

- `N`: The number of splay phases to generate.
- `amplitude`: The maximum perturbation amplitude.
- `seed`: The seed for the random number generator.

**Example:**

```cpp
dVec perturbed_splay_phases = splay_perturbed(16, 0.1, 333);
```

### `module_by_condition`

Generates phases/frequencies for a single module based on a specified distribution type.

```cpp
inline dVec module_by_condition(
    size_t                module_size,
    const std::string&    dist_type,
    double                a,
    double                b,
    unsigned              seed
);
```

- `module_size`: The number of values to generate for the module.
- `dist_type`: A string specifying the distribution type ("uniform", "normal", "cauchy", "exponential", "circle", "splay", "splay_perturbed").
- `a`, `b`: Parameters whose meaning depends on `dist_type` (e.g., min/max for uniform, mean/stddev for normal, location/scale for cauchy, lambda for exponential, amplitude for splay_perturbed).
- `seed`: The seed for the random number generator.

**Example:**

```cpp
dVec module_phases = module_by_condition(50, "normal", 0.0, 1.0, 444);
```

### `identical_modules`

Generates phases/frequencies for one module and then replicates these values across multiple identical modules.

```cpp
inline dVec identical_modules(
    size_t               N_per_module,
    size_t               num_modules,
    const std::string&   dist_type,
    double               a,
    double               b,
    unsigned             seed
);
```

- `N_per_module`: The number of nodes (phases/frequencies) per module.
- `num_modules`: The total number of identical modules.
- `dist_type`: A string specifying the distribution type (same as `module_by_condition`).
- `a`, `b`: Parameters for the distribution (same as `module_by_condition`).
- `seed`: The seed for the random number generator.

**Example:**

```cpp
dVec all_phases = identical_modules(20, 5, "uniform", -PI, PI, 555);
```
