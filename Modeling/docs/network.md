# Network Module

The `network` module provides a suite of functions for generating various network topologies, converting between dense and sparse matrix representations, and calculating network properties like density and degrees.

## Matrix Conversion and Properties

### `dense_to_sparse`

Converts a dense adjacency matrix to a sparse matrix representation.

```cpp
inline SparseMatrix dense_to_sparse(const Matrix &adj);
```

- `adj`: The dense adjacency matrix (N x N).

**Returns:** A `SparseMatrix` representation of the input matrix.

**Example:**

```cpp
Matrix dense_adj = {{0, 1, 0}, {0, 0, 2}, {3, 0, 0}};
SparseMatrix sparse_adj = dense_to_sparse(dense_adj);
// sparse_adj will represent the connections efficiently
```

### `density`

Computes the density (fraction of nonzero off-diagonal entries) for a dense adjacency matrix.

```cpp
inline double density(
    const Matrix &adj,
    double threshold = 1e-12
);
```

- `adj`: The dense adjacency matrix (N x N).
- `threshold`: Minimum absolute value to consider an entry as nonzero (default: `1e-12`).

**Returns:** The density of the network.

**Example:**

```cpp
Matrix adj_matrix = {{0, 1, 0}, {0, 0, 2}, {0, 0, 0}};
double network_density = density(adj_matrix);
// network_density will be (2 / (3*2)) = 0.333...
```

### `dense_to_sparse_conditional`

Converts a dense matrix to a sparse matrix only if its density is below a specified threshold; otherwise, it just returns the density.

```cpp
inline std::pair<bool, double> dense_to_sparse_conditional(
    const Matrix &adj,
    SparseMatrix &sparse_adj,
    double density_threshold = 0.5,
    double zero_threshold = 1e-12
);
```

- `adj`: Dense adjacency matrix (N x N).
- `sparse_adj`: Output sparse matrix (filled only if density <= `density_threshold`).
- `density_threshold`: Maximum density to allow conversion to sparse (default: `0.5`).
- `zero_threshold`: Minimum absolute value to consider as nonzero (default: `1e-12`).

**Returns:** A `std::pair<bool, double>` where the first element is `true` if `sparse_adj` was filled, and `false` otherwise. The second element is the calculated density.

**Example:**

```cpp
Matrix adj = {{0, 1, 0}, {0, 0, 2}, {3, 0, 0}};
SparseMatrix sparse_result;
auto result = dense_to_sparse_conditional(adj, sparse_result, 0.4);
if (result.first) {
    // sparse_result is now filled
} else {
    // Density was too high: result.second contains the density
}
```

## Network Generation Functions

### `random`

Generates a random dense adjacency matrix where each off-diagonal edge weight is uniformly distributed within `[min_weight, max_weight]`.

```cpp
inline Matrix random(
    size_t   N          = 100,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
);
```

- `N`: Number of nodes.
- `min_weight`: Minimum edge weight.
- `max_weight`: Maximum edge weight.
- `seed`: Random seed.

**Returns:** A randomly generated dense adjacency matrix (N x N).

**Example:**

```cpp
Matrix rand_net = random(50, 0.1, 0.5, 123);
```

### `random_symmetric`

Generates a random symmetric dense adjacency matrix. Edge weights are uniformly distributed within `[min_weight, max_weight]` and `adj[i][j] = adj[j][i]`.

```cpp
inline Matrix random_symmetric(
    size_t   N          = 100,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
);
```

- `N`: Number of nodes.
- `min_weight`: Minimum edge weight.
- `max_weight`: Maximum edge weight.
- `seed`: Random seed.

**Returns:** A randomly generated symmetric dense adjacency matrix (N x N).

**Example:**

```cpp
Matrix sym_rand_net = random_symmetric(50, 0.1, 0.5, 123);
```

### `erdos_renyi`

Generates an Erdos-Renyi random graph where each edge exists with probability `p`. If an edge exists, its weight is uniformly distributed within `[min_weight, max_weight]`.

```cpp
inline Matrix erdos_renyi(
    size_t   N          = 100,
    double   p          = 0.5,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
);
```

- `N`: Number of nodes.
- `p`: Connection probability.
- `min_weight`: Minimum edge weight.
- `max_weight`: Maximum edge weight.
- `seed`: Random seed.

**Returns:** An Erdos-Renyi adjacency matrix (N x N).

**Example:**

```cpp
Matrix er_net = erdos_renyi(100, 0.1, 0.5, 1.0, 456);
```

### `erdos_renyi_uniform`

Generates an Erdos-Renyi random graph with *exactly* `N*(N-1)*p` edges. These edges are chosen uniformly at random from all possible non-self-loop edges, and their weights are uniformly distributed.

```cpp
inline Matrix erdos_renyi_uniform(
    size_t   N          = 100,
    double   p          = 0.5,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
);
```

- `N`: Number of nodes.
- `p`: Connection probability (determines the number of edges).
- `min_weight`: Minimum edge weight.
- `max_weight`: Maximum edge weight.
- `seed`: Random seed.

**Returns:** An Erdos-Renyi adjacency matrix with a fixed number of edges.

**Example:**

```cpp
Matrix er_uniform_net = erdos_renyi_uniform(100, 0.05, 0.5, 1.0, 456);
```

### `erdos_renyi_symmetric`

Generates a symmetric Erdos-Renyi random graph. Each potential undirected edge exists with probability `p`, and if it exists, its weight is uniformly distributed.

```cpp
inline Matrix erdos_renyi_symmetric(
    size_t   N          = 100,
    double   p          = 0.5,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
);
```

- `N`: Number of nodes.
- `p`: Connection probability.
- `min_weight`: Minimum edge weight.
- `max_weight`: Maximum edge weight.
- `seed`: Random seed.

**Returns:** A symmetric Erdos-Renyi adjacency matrix (N x N).

**Example:**

```cpp
Matrix er_sym_net = erdos_renyi_symmetric(100, 0.1, 0.5, 1.0, 456);
```

### `erdos_renyi_symmetric_uniform`

Generates a symmetric Erdos-Renyi random graph with *exactly* `N*(N-1)*p/2` undirected edges. These edges are chosen uniformly at random, and their weights are uniformly distributed.

```cpp
inline Matrix erdos_renyi_symmetric_uniform(
    size_t   N          = 100,
    double   p          = 0.5,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
);
```

- `N`: Number of nodes.
- `p`: Connection probability (determines the number of edges).
- `min_weight`: Minimum edge weight.
- `max_weight`: Maximum edge weight.
- `seed`: Random seed.

**Returns:** A symmetric Erdos-Renyi adjacency matrix with a fixed number of undirected edges.

**Example:**

```cpp
Matrix er_sym_uniform_net = erdos_renyi_symmetric_uniform(100, 0.05, 0.5, 1.0, 456);
```

### `small_world`

Generates a Watts-Strogatz small-world network. Starts with a ring lattice where each node is connected to its `k` nearest neighbors, then rewires each edge with probability `beta`.

```cpp
inline Matrix small_world(
    size_t N      = 100,
    size_t k      = 4,
    double beta   = 0.5,
    double weight = 1.0,
    unsigned seed = 42
);
```

- `N`: Number of nodes.
- `k`: Each node is connected to `k` nearest neighbors in the ring topology.
- `beta`: Rewiring probability (0.0 for pure ring, 1.0 for random).
- `weight`: Edge weight for all connections.
- `seed`: Random seed.

**Returns:** A small-world adjacency matrix (N x N).

**Example:**

```cpp
Matrix sw_net = small_world(100, 6, 0.2, 1.0, 789);
```

### `small_world_sparse`

Generates a Watts-Strogatz small-world network and returns it as a `SparseMatrix`.

```cpp
inline SparseMatrix small_world_sparse(
    size_t   N,
    size_t   k,
    double   beta,
    double   weight,
    unsigned seed = 42
);
```

- Parameters are the same as `small_world`.

**Returns:** A sparse small-world adjacency matrix.

**Example:**

```cpp
SparseMatrix sw_sparse_net = small_world_sparse(100, 6, 0.2, 1.0, 789);
```

### `small_world_directed`

Generates a directed Watts-Strogatz small-world network. Each node has `k` outgoing edges in the initial ring, which are then rewired with probability `beta`.

```cpp
inline Matrix small_world_directed(
    size_t N      = 100,
    size_t k      = 4,
    double beta   = 0.5,
    double weight = 1.0,
    unsigned seed = 42
);
```

- `N`: Number of nodes.
- `k`: Each node has `k` outgoing edges in the ring topology.
- `beta`: Rewiring probability.
- `weight`: Edge weight for all connections.
- `seed`: Random seed.

**Returns:** A directed small-world adjacency matrix (N x N).

**Example:**

```cpp
Matrix sw_directed_net = small_world_directed(100, 4, 0.3, 1.0, 901);
```

### `multilayered`

Combines multiple network layers into a single block-diagonal adjacency matrix, representing a multilayered network where layers are uncoupled.

```cpp
inline Matrix multilayered(const std::vector<Matrix> &layers);
```

- `layers`: A vector of adjacency matrices, where each matrix represents a layer.

**Returns:** A single adjacency matrix representing the combined multilayered network.

**Example:**

```cpp
Matrix layer1 = {{0, 1}, {1, 0}};
Matrix layer2 = {{0, 0.5}, {0.5, 0}};
std::vector<Matrix> layers = {layer1, layer2};
Matrix multi_net = multilayered(layers);
// multi_net will be a 4x4 block-diagonal matrix
```

### `modular`

Generates a modular network where nodes are divided into modules. Connections are dense within modules and sparse between modules.

```cpp
inline Matrix modular(
    size_t   module_size = 100,
    size_t   num_modules = 10,
    double   p_in        = 0.9,
    double   p_out       = 0.1,
    double   in_weight   = 1.0,
    double   out_weight  = 1.0,
    unsigned seed        = 42
);
```

- `module_size`: Number of nodes per module.
- `num_modules`: Number of modules.
- `p_in`: Probability of within-module connection.
- `p_out`: Probability of between-module connection.
- `in_weight`: Weight for within-module connections.
- `out_weight`: Weight for between-module connections.
- `seed`: Random seed.

**Returns:** A modular adjacency matrix (N x N), where `N = module_size * num_modules`.

**Example:**

```cpp
Matrix mod_net = modular(20, 5, 0.8, 0.05, 1.0, 0.5, 1122);
```

### `hierarchical`

Generates a hierarchical network with recursively nested modules and level-dependent connection probabilities.

```cpp
inline Matrix hierarchical(
    size_t   N               = 100, // Base module size at the lowest level
    size_t   levels          = 2, // Number of hierarchical levels (0 returns Erdos-Renyi N x N matrix)
    double   p_in            = 0.9, // Base within-module connection probability (decays with level)
    double   p_out           = 0.1, // Base between-module connection probability (decays with level)
    double   in_weight       = 1.0,
    double   out_weight      = 1.0,
    double   level_decay     = 0.5, // Decay factor for connection probabilities at each level (0 < level_decay <= 1)
    unsigned seed            = 42,
    size_t   base_module_num = 2 // Number of modules in the base level of the hierarchy
);
```

- `N`: Base module size at the lowest level.
- `levels`: Number of hierarchical levels (0 returns an Erdos-Renyi N x N matrix).
- `p_in`: Base within-module connection probability (decays with level).
- `p_out`: Base between-module connection probability (decays with level).
- `in_weight`: Weight for within-module connections.
- `out_weight`: Weight for between-module connections.
- `level_decay`: Decay factor for connection probabilities at each level (0 < `level_decay` <= 1).
- `seed`: Random seed.
- `base_module_num`: Number of modules in the base level of the hierarchy.

**Returns:** A hierarchical adjacency matrix.

**Example:**

```cpp
Matrix hier_net = hierarchical(10, 3, 0.9, 0.01, 1.0, 0.1, 0.7, 3344, 2);
```

### `pick_topology`

A high-level function to create a network topology based on a string identifier and a set of generic parameters.

```cpp
inline Matrix pick_topology(
    const std::string &topology_type = "random",
    size_t   N               = 100,
    double   a               = 0.0,
    double   b               = 1.0,
    double   c               = 1.0,
    unsigned seed            = 42,
    double   d               = 1.0,
    double   e               = 1.0,
    double   f               = 0.5,
    size_t   base_module_num = 2
);
```

- `topology_type`: Type of network to create ("random", "erdos_renyi", "small_world", "modular", "hierarchical", etc.).
- `N`: Number of nodes (interpreted differently based on `topology_type`).
- `a, b, c, d, e, f`: Generic parameters whose meaning depends on `topology_type` (refer to individual topology functions for details).
- `seed`: Random seed.
- `base_module_num`: Specific to hierarchical topology.

**Returns:** The adjacency matrix for the selected topology.

**Example:**

```cpp
// Create an Erdos-Renyi network with N=50, p=0.1, weights between 0.5 and 1.0
Matrix custom_net = pick_topology("erdos_renyi", 50, 0.1, 0.5, 1.0, 567);

// Create a small-world network with N=80, k=6, beta=0.3, weight=1.0
Matrix sw_picked = pick_topology("small_world", 80, 6.0, 0.3, 1.0, 890);
```

### `effective_multiplex`

Calculates the weighted sum of multiple network layers to form a single effective adjacency matrix.

```cpp
inline Matrix effective_multiplex(
    const std::vector<Matrix> &layers,
    const dVec                 &layer_weights
);
```

- `layers`: Vector of adjacency matrices (all must be N x N).
- `layer_weights`: Vector of weights, one for each layer.

**Returns:** The weighted sum of layers as a single adjacency matrix.

**Example:**

```cpp
Matrix layer_A = {{0, 1}, {1, 0}};
Matrix layer_B = {{0, 0}, {1, 0}};
std::vector<Matrix> layers_vec = {layer_A, layer_B};
dVec weights = {0.6, 0.4};
Matrix effective_adj = effective_multiplex(layers_vec, weights);
```

### `generate_multiplex_network_layers`

Generates multiple network layers for a multiplex network based on specified topology types and parameters for each layer.

```cpp
inline std::vector<Matrix> generate_multiplex_network_layers(
    size_t N,
    size_t num_layers,
    const  std::vector<std::string> &layer_generation_types,
    const  Matrix                   &layer_params,
    const  std::vector<unsigned>    &layer_seeds
);
```

- `N`: Number of nodes (consistent across all layers).
- `num_layers`: Number of layers to generate.
- `layer_generation_types`: Vector of topology types for each layer.
- `layer_params`: Matrix of parameters for each layer (one row per layer, parameters correspond to `pick_topology`'s `a, b, c, ...`).
- `layer_seeds`: Vector of random seeds for each layer.

**Returns:** A vector of adjacency matrices, one per layer.

**Example:**

```cpp
size_t N_nodes = 50;
size_t num_L = 2;
std::vector<std::string> types = {"erdos_renyi", "small_world"};
Matrix params_mat = {{0.1, 0.5, 1.0}, {6.0, 0.2, 1.0}};
std::vector<unsigned> seeds = {111, 222};

std::vector<Matrix> multiplex_layers = generate_multiplex_network_layers(
    N_nodes, num_L, types, params_mat, seeds
);
```

## Degree Calculation Functions

### `in_degrees` (Dense Matrix)

Calculates the in-degrees (number of incoming edges) for each node in a dense adjacency matrix.

```cpp
inline wVec in_degrees(const Matrix& adj, double threshold = 1e-12);
```

- `adj`: Dense adjacency matrix (N x N).
- `threshold`: Minimum absolute value to consider as an edge (default: `1e-12`).

**Returns:** A vector of in-degrees for each node.

**Example:**

```cpp
Matrix adj_d = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}};
wVec in_deg = in_degrees(adj_d);
// in_deg will be {1, 1, 1}
```

### `out_degrees` (Dense Matrix)

Calculates the out-degrees (number of outgoing edges) for each node in a dense adjacency matrix.

```cpp
inline wVec out_degrees(const Matrix& adj, double threshold = 1e-12);
```

- `adj`: Dense adjacency matrix (N x N).
- `threshold`: Minimum absolute value to consider as an edge (default: `1e-12`).

**Returns:** A vector of out-degrees for each node.

**Example:**

```cpp
Matrix adj_d = {{0, 1, 0}, {0, 0, 1}, {1, 0, 0}};
wVec out_deg = out_degrees(adj_d);
// out_deg will be {1, 1, 1}
```

### `in_degrees` (Sparse Matrix)

Calculates the in-degrees for each node in a sparse adjacency matrix.

```cpp
inline wVec in_degrees(const SparseMatrix& sparse, double threshold = 1e-12);
```

- `sparse`: Sparse adjacency matrix (N x N).
- `threshold`: Minimum absolute value to consider as an edge (default: `1e-12`).

**Returns:** A vector of in-degrees for each node.

**Example:**

```cpp
// Assuming sparse_adj is a SparseMatrix
wVec in_deg_sparse = in_degrees(sparse_adj);
```

### `out_degrees` (Sparse Matrix)

Calculates the out-degrees for each node in a sparse adjacency matrix.

```cpp
inline wVec out_degrees(const SparseMatrix& sparse, double threshold = 1e-12);
```

- `sparse`: Sparse adjacency matrix (N x N).
- `threshold`: Minimum absolute value to consider as an edge (default: `1e-12`).

**Returns:** A vector of out-degrees for each node.

**Example:**

```cpp
// Assuming sparse_adj is a SparseMatrix
wVec out_deg_sparse = out_degrees(sparse_adj);
```
