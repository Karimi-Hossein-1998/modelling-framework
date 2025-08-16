# Models Module

The `models` module contains implementations of various mathematical models, primarily focusing on the Kuramoto model and its variations.

## Kuramoto Model

The Kuramoto model describes the synchronization of a population of coupled oscillators. The general form of the Kuramoto model with phase-lag is given by:

$$\frac{d\theta_i}{dt} = \omega_i + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j - \theta_i - \alpha)$$

where:
- $\theta_i$ is the phase of oscillator $i$.
- $\omega_i$ is the natural frequency of oscillator $i$.
- $K$ is the global coupling strength.
- $N$ is the total number of oscillators.
- $A_{ij}$ is the element of the adjacency matrix representing the coupling strength from oscillator $j$ to oscillator $i$.
- $\alpha$ is the phase-lag parameter.

### `kuramoto_general`

Calculates the derivatives of phases for a general Kuramoto model with a dense adjacency matrix and phase-lag.

```cpp
inline dVec kuramoto_general(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              K,
    const Matrix&       adj,
    double              alpha
);
```

- `time`: Current simulation time (unused in this specific model, but common for ODE systems).
- `theta`: Current phases of the oscillators.
- `omega`: Natural frequencies of the oscillators.
- `K`: Global coupling strength.
- `adj`: Dense adjacency matrix (N x N) representing connections between oscillators.
- `alpha`: Phase-lag parameter.

**Returns:** A `dVec` containing the derivatives `d(theta)/dt` for each oscillator.

**Example:**

```cpp
dVec phases = {0.1, 0.2, 0.3};
dVec freqs = {1.0, 1.1, 1.2};
Matrix adj_mat = {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
dVec dtheta = kuramoto_general(0.0, phases, freqs, 10.0, adj_mat, 0.0);
```

### `kuramoto_general_parallel`

Parallelized version of `kuramoto_general` for improved performance on multi-core systems.

```cpp
inline dVec kuramoto_general_parallel(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              K,
    const Matrix&       adj,
    double              alpha
);
```

- Parameters are the same as `kuramoto_general`.

**Returns:** A `dVec` containing the derivatives `d(theta)/dt` for each oscillator.

### `kuramoto_sparse`

Calculates the derivatives of phases for a general Kuramoto model using a sparse adjacency matrix, which is more efficient for sparsely connected networks.

```cpp
inline dVec kuramoto_sparse(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              K,
    const SparseMatrix& sparse_adj,
    double              alpha
);
```

- `sparse_adj`: Sparse adjacency matrix representation.
- Other parameters are the same as `kuramoto_general`.

**Returns:** A `dVec` containing the derivatives `d(theta)/dt` for each oscillator.

**Example:**

```cpp
// Assuming sparse_adj is a SparseMatrix representation
dVec dtheta_sparse = kuramoto_sparse(0.0, phases, freqs, 10.0, sparse_adj, 0.0);
```

### `kuramoto_sparse_parallel`

Parallelized version of `kuramoto_sparse`.

```cpp
inline dVec kuramoto_sparse_parallel(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              K,
    const SparseMatrix& sparse_adj,
    double              alpha
);
```

- Parameters are the same as `kuramoto_sparse`.

**Returns:** A `dVec` containing the derivatives `d(theta)/dt` for each oscillator.

### `kuramoto_special_modular`

Calculates the derivatives for a special modular Kuramoto model, where coupling strengths differ for within-module and between-module connections.

```cpp
inline dVec kuramoto_special_modular(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              intra_K, // Coupling strength within modules
    double              inter_K, // Coupling strength between modules
    double              alpha,
    const size_t        module_size
);
```

- `intra_K`: Coupling strength for connections within the same module.
- `inter_K`: Coupling strength for connections between different modules.
- `module_size`: The number of oscillators in each module.
- Other parameters are similar to `kuramoto_general`.

**Returns:** A `dVec` containing the derivatives `d(theta)/dt` for each oscillator.

**Example:**

```cpp
// Assuming phases, freqs are for a modular network
dVec dtheta_modular = kuramoto_special_modular(0.0, phases, freqs, 10.0, 1.0, 0.0, 20);
```

### `kuramoto_special_modular_parallel`

Parallelized version of `kuramoto_special_modular`.

```cpp
inline dVec kuramoto_special_modular_parallel(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              intra_K,
    double              inter_K,
    double              alpha,
    const size_t        module_size
);
```

- Parameters are the same as `kuramoto_special_modular`.

**Returns:** A `dVec` containing the derivatives `d(theta)/dt` for each oscillator.

## Order Parameter Calculation

The order parameter is a key metric in the Kuramoto model, quantifying the level of synchronization in the system. It is a complex number $r e^{i\psi}$, where $r$ is the coherence and $\psi$ is the average phase.

$$r e^{i\psi} = \frac{1}{N} \sum_{j=1}^{N} e^{i\theta_j}$$

### `calculate_order` (for `dVec` phases)

Calculates the global order parameter (r, mean_sin, mean_cos) for a given set of phases.

```cpp
inline dVec calculate_order(const dVec& phases);
```

- `phases`: A `dVec` of oscillator phases.

**Returns:** A `dVec` of size 3: `[mean_sin, mean_cos, coherence_r]`. `coherence_r` is the magnitude of the order parameter.

**Example:**

```cpp
dVec current_phases = {0.1, 0.2, 0.3, 0.4};
dVec order_params = calculate_order(current_phases);
// order_params[2] will contain the coherence 'r'
```

### `calculate_order` (for `SolverResults`)

Calculates the global order parameter over time for a simulation's results.

```cpp
inline Matrix calculate_order(const SolverResults& results);
```

- `results`: A `SolverResults` object containing the time points and solution (phases over time).

**Returns:** A `Matrix` where each row represents `[mean_sin, mean_cos, coherence_r, time]` for a given time point.

**Example:**

```cpp
// Assuming SolverResults results from a simulation
Matrix global_order_evolution = calculate_order(results);
```

### `calculate_order_per_node` (Dense Matrix)

Calculates the local order parameter for a specific node based on its connections in a dense adjacency matrix.

```cpp
inline dVec calculate_order_per_node(
    const  dVec& phases,
    size_t node_number,
    const  Matrix& adj,
    double tolerance = 1e-12
);
```

- `phases`: Current phases of the oscillators.
- `node_number`: The index of the node for which to calculate the local order parameter.
- `adj`: Dense adjacency matrix.
- `tolerance`: Threshold for considering an edge as existing.

**Returns:** A `dVec` of size 3: `[local_mean_sin, local_mean_cos, local_coherence_r]` for the specified node.

**Example:**

```cpp
dVec node_phases = {0.1, 0.2, 0.3};
Matrix node_adj = {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}};
dVec node0_order = calculate_order_per_node(node_phases, 0, node_adj);
```

### `calculate_order_per_node` (Dense Matrix, all nodes)

Calculates the local order parameter for all nodes in a network with a dense adjacency matrix.

```cpp
inline Matrix calculate_order_per_node(
    const  dVec& phases,
    const  Matrix& adj,
    double tolerance = 1e-12
);
```

- `phases`: Current phases of the oscillators.
- `adj`: Dense adjacency matrix.
- `tolerance`: Threshold for considering an edge as existing.

**Returns:** A `Matrix` where each row is `[local_mean_sin, local_mean_cos, local_coherence_r]` for each node.

### `calculate_order_per_node` (Sparse Matrix)

Calculates the local order parameter for a specific node based on its connections in a sparse adjacency matrix.

```cpp
inline dVec calculate_order_per_node(
    const dVec& phases,
    size_t node,
    const SparseMatrix& sparse_adj
);
```

- `sparse_adj`: Sparse adjacency matrix.
- Other parameters are similar to the dense version.

**Returns:** A `dVec` of size 3: `[local_mean_sin, local_mean_cos, local_coherence_r]` for the specified node.

### `calculate_order_per_node` (Sparse Matrix, all nodes)

Calculates the local order parameter for all nodes in a network with a sparse adjacency matrix.

```cpp
inline Matrix calculate_order_per_node(
    const dVec& phases,
    const SparseMatrix& sparse_adj
);
```

- `sparse_adj`: Sparse adjacency matrix.
- Other parameters are similar to the dense version.

**Returns:** A `Matrix` where each row is `[local_mean_sin, local_mean_cos, local_coherence_r]` for each node.

### `calculate_order_per_module` (single time point)

Calculates the order parameter for each module in a modular network at a single time point.

```cpp
inline dVec calculate_order_per_module(
    const  dVec& phases,
    const  wVec&  module_assignments,
    size_t num_modules
);
```

- `phases`: Current phases of the oscillators.
- `module_assignments`: A vector indicating which module each oscillator belongs to.
- `num_modules`: The total number of modules.

**Returns:** A `dVec` containing `[mean_sin, mean_cos, coherence_r]` for each module, followed by the global order parameter for the entire network.

**Example:**

```cpp
dVec mod_phases = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
wVec assignments = {0, 0, 1, 1, 2, 2}; // 3 modules
dVec module_orders = calculate_order_per_module(mod_phases, assignments, 3);
// module_orders will contain order params for module 0, then module 1, then module 2, then global.
```

### `calculate_order_per_module` (over time)

Calculates the per-module order parameter evolution for a simulation's results.

```cpp
inline Matrix calculate_order_per_module(
    const SolverResults& results,
    const wVec&          module_assignments,
    size_t               num_modules
);
```

- `results`: A `SolverResults` object.
- `module_assignments`: A vector indicating which module each oscillator belongs to.
- `num_modules`: The total number of modules.

**Returns:** A `Matrix` where each row contains the order parameters for all modules and the global order parameter, followed by the time point.

### `calculate_order_hierarchical`

Calculates hierarchical order parameters, aggregating from lower levels to higher levels of a hierarchical network structure.

```cpp
inline Matrix calculate_order_hierarchical(
    const dVec&      phases,
    const Vec<wVec>& level_assignments,
    const wVec&      num_groups_per_level
);
```

- `phases`: Current phases of the oscillators.
- `level_assignments`: A `Vec<wVec>` where each inner `wVec` contains the group assignment for each oscillator at a specific hierarchical level.
- `num_groups_per_level`: A `wVec` indicating the number of groups at each hierarchical level.

**Returns:** A `Matrix` where each row represents a hierarchical level, and columns contain the order parameters `[mean_sin, mean_cos, coherence_r]` for each group at that level.

**Example:**

```cpp
// Assuming phases, level_assignments, and num_groups_per_level are set up for a hierarchical network
Matrix hierarchical_orders = calculate_order_hierarchical(phases, level_assignments, num_groups_per_level);
```
