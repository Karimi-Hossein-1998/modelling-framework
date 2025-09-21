#pragma once
#include "../typedefs/header.hpp"

// -----------------------------------------------------------------------------
// Convert dense Matrix to sparse Matrix
/*----------------------------------------------------------*/
// Parameters:
//   adj: Dense adjacency matrix (N x N)
// Returns:
//   SparseMatrix representation of the input matrix
inline SparseMatrix dense_to_sparse(const Matrix &adj)
{
    size_t N = adj.size();
    SparseMatrix sparse(N);
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            if (adj[i][j] != 0.0)
            {
                sparse.rows[i].emplace_back(j, adj[i][j]);
            }
        }
    }
    return sparse;
}

// -----------------------------------------------------------------------------
// Compute density (fraction of nonzero off-diagonal entries) for a dense Matrix
/*----------------------------------------------------------*/
// Parameters:
//   adj: Dense adjacency matrix (N x N)
//   threshold: Minimum absolute value to consider an entry as nonzero (default 1e-12)
// Returns:
//   Density (fraction of nonzero off-diagonal entries)
inline double density(
    const Matrix &adj,
    double threshold = 1e-12)
{
    size_t N = adj.size(), nonzero = 0, total = N * (N - 1);
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            if (i != j && std::abs(adj[i][j]) > threshold)
            {
                ++nonzero;
            }
        }
    }
    return total ? static_cast<double>(nonzero) / total : 0.0;
}

// -----------------------------------------------------------------------------
// Convert dense Matrix to sparse Matrix if sparse enough, else return density
/*----------------------------------------------------------*/
// Parameters:
//   adj: Dense adjacency matrix (N x N)
//   sparse_adj: Output sparse matrix (only filled if density <= density_threshold)
//   density_threshold: Maximum density to allow conversion to sparse (default 0.5)
//   zero_threshold: Minimum absolute value to consider as nonzero (default 1e-12)
// Returns:
//   Pair: (true if sparse_adj is filled, false if only density is set; density value)
inline std::pair<bool, double> dense_to_sparse_conditional(
    const Matrix &adj,
    SparseMatrix &sparse_adj,
    double density_threshold = 0.5,
    double zero_threshold = 1e-12)
{
    size_t N = adj.size();
    auto density_measure = density(adj, zero_threshold);
    if (density_measure <= density_threshold)
    {
        sparse_adj = SparseMatrix(N);
        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                if (std::abs(adj[i][j]) > zero_threshold)
                {
                    sparse_adj.rows[i].emplace_back(j, adj[i][j]);
                }
            }
        }
        return {true, density_measure};
    }
    else
    {
        return {false, density_measure};
    }
}

// Random Matrix: each edge weight is random in [min_weight, max_weight]
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated dense adjacency matrix (N x N)
inline Matrix random(
    size_t   N          = 100,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0)
    {
        throw std::invalid_argument("[random_Matrix] Number of nodes N cannot be zero.");
    }
    if (min_weight > max_weight)
    {
        throw std::invalid_argument("[random_Matrix] min_weight (" + std::to_string(min_weight) + ") cannot be greater than max_weight (" + std::to_string(max_weight) + ").");
    }

    // Original logic
    std::mt19937                           rng(seed);
    // uniform_real_distribution handles min_weight == max_weight correctly.
    std::uniform_real_distribution<double> dist(min_weight, max_weight);
    Matrix                                 adj(N, dVec(N));
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            adj[i][j] = (i == j) ? 0.0 : dist(rng);
        }
    }
    return adj;
}

// Random Symmetric Matrix: each edge weight is random in [min_weight, max_weight]
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated dense adjacency matrix (N x N)
inline Matrix random_symmetric(
    size_t   N          = 100,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0)
    {
        throw std::invalid_argument("[random_Matrix] Number of nodes N cannot be zero.");
    }
    if (min_weight > max_weight)
    {
        throw std::invalid_argument("[random_Matrix] min_weight (" + std::to_string(min_weight) + ") cannot be greater than max_weight (" + std::to_string(max_weight) + ").");
    }

    // Original logic
    std::mt19937                           rng(seed);
    // uniform_real_distribution handles min_weight == max_weight correctly.
    std::uniform_real_distribution<double> dist(min_weight, max_weight);
    Matrix                                 adj(N, dVec(N));
    for (size_t i = 1; i < N; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            adj[i][j] = dist(rng);
            adj[j][i] = adj[i][j];
        }
    }
    return adj;
}

// Erdos-Renyi Matrix: each edge exists with probability p, weight random in [min_weight, max_weight]
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   p: Connection probability (default 0.5)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated Erdos-Renyi adjacency matrix (N x N)
inline Matrix erdos_renyi(
    size_t   N          = 100,
    double   p          = 0.5,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0)
    {
        throw std::invalid_argument("[erdos_renyi_Matrix] Number of nodes N cannot be zero.");
    }
    if (p < 0.0 || p > 1.0)
    {
        throw std::invalid_argument("[erdos_renyi_Matrix] Probability p (" + std::to_string(p) + ") must be between 0.0 and 1.0.");
    }
    if (min_weight > max_weight)
    {
        throw std::invalid_argument("[erdos_renyi_Matrix] min_weight (" + std::to_string(min_weight) + ") cannot be greater than max_weight (" + std::to_string(max_weight) + ").");
    }

    // Original logic
    std::mt19937                           rng(seed);
    std::uniform_real_distribution<double> weight_dist(min_weight, max_weight); // Renamed dist
    std::bernoulli_distribution            edge_dist(p);
    Matrix                                 adj(N, dVec(N, 0.0));
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            if (i != j && edge_dist(rng))
            {
                adj[i][j] = weight_dist(rng);
            }
        }

    }
    return adj;
}

// Erdos-Renyi Matrix (Uniform): exactly N*(N-1)*p edges with random weights
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   p: Connection probability (default 0.5)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated Erdos-Renyi adjacency matrix (N x N) with exactly N*(N-1)*p edges
inline Matrix erdos_renyi_uniform(
    size_t   N          = 100,
    double   p          = 0.5,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0)
    {
        throw std::invalid_argument("[erdos_renyi_uniform] Number of nodes N cannot be zero.");
    }
    if (p < 0.0 || p > 1.0)
    {
        throw std::invalid_argument("[erdos_renyi_uniform] Probability p (" + std::to_string(p) + ") must be between 0.0 and 1.0.");
    }
    if (min_weight > max_weight)
    {
        throw std::invalid_argument("[erdos_renyi_uniform] min_weight (" + std::to_string(min_weight) + ") cannot be greater than max_weight (" + std::to_string(max_weight) + ").");
    }

    std::mt19937                           rng(seed);
    std::uniform_real_distribution<double> weight_dist(min_weight, max_weight);
    Matrix                                 adj(N, dVec(N, 0.0));

    // Calculate exact number of edges needed
    size_t total_edges = static_cast<size_t>(N * (N - 1) * p);
    
    // Create vector of all possible edges
    std::vector<std::pair<size_t, size_t>> possible_edges;
    possible_edges.reserve(N * (N - 1));  // Reserve space for all possible edges
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            if (i != j)  // No self-loops
            {
                possible_edges.emplace_back(i, j);
            }
        }
    }

    // Shuffle the edges
    std::shuffle(possible_edges.begin(), possible_edges.end(), rng);

    // Create exactly total_edges edges
    for (size_t e = 0; e < total_edges; ++e)
    {
        size_t i  = possible_edges[e].first;
        size_t j  = possible_edges[e].second;
        adj[i][j] = weight_dist(rng);
    }

    return adj;
}

// Symmetric Erdos-Renyi Matrix: each edge exists with probability p, weight random in [min_weight, max_weight]
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   p: Connection probability (default 0.5)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated Erdos-Renyi (Symmetric) adjacency matrix (N x N)
inline Matrix erdos_renyi_symmetric(
    size_t   N          = 100,
    double   p          = 0.5,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0)
    {
        throw std::invalid_argument("[erdos_renyi_symmetric_Matrix] Number of nodes N cannot be zero.");
    }
    if (p < 0.0 || p > 1.0)
    {
        throw std::invalid_argument("[erdos_renyi_symmetric_Matrix] Probability p (" + std::to_string(p) + ") must be between 0.0 and 1.0.");
    }
    if (min_weight > max_weight)
    {
        throw std::invalid_argument("[erdos_renyi_symmetric_Matrix] min_weight (" + std::to_string(min_weight) + ") cannot be greater than max_weight (" + std::to_string(max_weight) + ").");
    }

    // Original logic
    std::mt19937                           rng(seed);
    std::uniform_real_distribution<double> weight_dist(min_weight, max_weight); // Renamed dist
    std::bernoulli_distribution            edge_dist(p);
    Matrix                                 adj(N, dVec(N, 0.0));
    for (size_t i = 1; i < N; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            if (edge_dist(rng))
            {
                adj[i][j] = weight_dist(rng);
                adj[j][i] = adj[i][j];
            }
        }
    }
    return adj;
}

// Symmetric Erdos-Renyi Matrix (Uniform): exactly N*(N-1)*p/2 edges with random weights
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   p: Connection probability (default 0.5)
//   min_weight: Minimum edge weight (default 0.0)
//   max_weight: Maximum edge weight (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Randomly generated Erdos-Renyi (Symmetric) adjacency matrix (N x N) with exactly N*(N-1)*p/2 edges
inline Matrix erdos_renyi_symmetric_uniform(
    size_t   N          = 100,
    double   p          = 0.5,
    double   min_weight = 0.0,
    double   max_weight = 1.0,
    unsigned seed       = 42
)
{
    // Error Handling
    if (N == 0)
    {
        throw std::invalid_argument("[erdos_renyi_symmetric_uniform] Number of nodes N cannot be zero.");
    }
    if (p < 0.0 || p > 1.0)
    {
        throw std::invalid_argument("[erdos_renyi_symmetric_uniform] Probability p (" + std::to_string(p) + ") must be between 0.0 and 1.0.");
    }
    if (min_weight > max_weight)
    {
        throw std::invalid_argument("[erdos_renyi_symmetric_uniform] min_weight (" + std::to_string(min_weight) + ") cannot be greater than max_weight (" + std::to_string(max_weight) + ").");
    }

    std::mt19937                           rng(seed);
    std::uniform_real_distribution<double> weight_dist(min_weight, max_weight);
    Matrix                                 adj(N, dVec(N, 0.0));

    // Calculate exact number of edges needed
    size_t total_edges = static_cast<size_t>(N * (N - 1) * p / 2);
    
    // Create vector of all possible edges (only lower triangle)
    std::vector<std::pair<size_t, size_t>> possible_edges;
    possible_edges.reserve(N * (N - 1) / 2);  // Reserve space for all possible edges
    for (size_t i = 1; i < N; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            possible_edges.emplace_back(i, j);
        }
    }

    // Shuffle the edges
    std::shuffle(possible_edges.begin(), possible_edges.end(), rng);

    // Create exactly total_edges edges
    for (size_t e = 0; e < total_edges; ++e)
    {
        size_t i      = possible_edges[e].first;
        size_t j      = possible_edges[e].second;
        double weight = weight_dist(rng);
        adj[i][j]     = weight;
        adj[j][i]     = weight;  // Make it symmetric
    }

    return adj;
}

// Small-world Matrix (Watts-Strogatz model, ring lattice with rewiring)
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   k: Each node is connected to k nearest neighbors in ring topology (default 4)
//   beta: Rewiring probability (default 0.5)
//   weight: Edge weight for all connections (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Small-world adjacency matrix (N x N)
inline Matrix small_world(
    size_t N      = 100,
    size_t k      = 4,
    double beta   = 0.5,
    double weight = 1.0,
    unsigned seed = 42
)
{
    // Error Handling
    if (N == 0)
    {
        throw std::invalid_argument("[small_world] Number of nodes N cannot be zero.");
    }
    if (k >= N)
    {
        throw std::invalid_argument("[small_world] k (" + std::to_string(k) + ") must be less than N (" + std::to_string(N) + ").");
    }
    if (k == 0)
    {
        throw std::invalid_argument("[small_world] k must be greater than 0.");
    }
    if (beta < 0.0 || beta > 1.0)
    {
        throw std::invalid_argument("[small_world] Rewiring probability beta (" + std::to_string(beta) + ") must be between 0.0 and 1.0.");
    }

    std::mt19937 rng(seed);
    Matrix       adj(N, dVec(N, 0.0));

    // Initial ring lattice - ensure k/2 neighbors on each side
    size_t half_k = k / 2; // Integer division intentional
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 1; j <= half_k; ++j)
        {
            size_t right  = (i + j) % N;
            size_t left   = (i + N - j) % N;
            adj[i][right] = weight;
            adj[i][left]  = weight;
        }
        // Handle odd k by adding one more connection
        if (k % 2 == 1 && N > 2)
        {
            size_t extra  = (i + (k / 2 + 1)) % N;
            adj[i][extra] = weight;
        }
    }
    // Rewiring with improved logic
    std::bernoulli_distribution rewire_dist(beta);
    std::uniform_int_distribution<size_t> node_dist(0, N - 1);

    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 1; j <= half_k; ++j)
        {
            size_t neighbor = (i + j) % N;
            if (rewire_dist(rng))
            {
                size_t attempts     = 0;
                size_t max_attempts = N * static_cast<size_t>(std::sqrt(N)); // Reasonable limit to prevent infinite loops
                bool   found        = false;
                size_t new_neighbor;

                while (!found && attempts < max_attempts)
                {
                    new_neighbor              = node_dist(rng);
                    if (new_neighbor         != i &&
                        new_neighbor         != neighbor &&
                        adj[i][new_neighbor] == 0.0)
                    {
                        found = true;
                    }
                    attempts++;
                }

                if (found)
                {
                    adj[i][neighbor]     = 0.0;
                    adj[neighbor][i]     = 0.0;
                    adj[i][new_neighbor] = weight;
                    adj[new_neighbor][i] = weight;
                }
            }
        }

        // Handle rewiring for the extra connection in case of odd k
        if (k % 2 == 1 && N > 2)
        {
            size_t extra = (i + (k / 2 + 1)) % N;
            if (rewire_dist(rng))
            {
                size_t attempts     = 0;
                size_t max_attempts = N * static_cast<size_t>(std::sqrt(N));
                bool   found        = false;
                size_t new_neighbor;

                while (!found && attempts < max_attempts)
                {
                    new_neighbor              = node_dist(rng);
                    if (new_neighbor         != i &&
                        new_neighbor         != extra &&
                        adj[i][new_neighbor] == 0.0)
                    {
                        found = true;
                    }
                    attempts++;
                }

                if (found)
                {
                    adj[i][extra]        = 0.0;
                    adj[extra][i]        = 0.0;
                    adj[i][new_neighbor] = weight;
                    adj[new_neighbor][i] = weight;
                }
            }
        }
    }
    return adj;
}

// -----------------------------------------------------------------------------
// Sparse small-world Matrix (Watts-Strogatz, returns SparseMatrix)
inline SparseMatrix small_world_sparse(
    size_t   N,
    size_t   k,
    double   beta,
    double   weight,
    unsigned seed = 42
)
{
    auto   dense  = small_world(N, k, beta, weight, seed);
    auto   sparse = dense_to_sparse(dense);
    return sparse;
}

// Directed Small-world Matrix (Watts-Strogatz model, ring lattice with rewiring)
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (default 100)
//   k: Each node has k outgoing edges in ring topology (default 4)
//   beta: Rewiring probability (default 0.5)
//   weight: Edge weight for all connections (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Directed small-world adjacency matrix (N x N)
inline Matrix small_world_directed(
    size_t N      = 100,
    size_t k      = 4,
    double beta   = 0.5,
    double weight = 1.0,
    unsigned seed = 42
)
{
    // Error Handling
    if (N == 0)
    {
        throw std::invalid_argument("[small_world_directed] Number of nodes N cannot be zero.");
    }
    if (k >= N)
    {
        throw std::invalid_argument("[small_world_directed] k (" + std::to_string(k) + ") must be less than N (" + std::to_string(N) + ").");
    }
    if (k == 0)
    {
        throw std::invalid_argument("[small_world_directed] k must be greater than 0.");
    }
    if (beta < 0.0 || beta > 1.0)
    {
        throw std::invalid_argument("[small_world_directed] Rewiring probability beta (" + std::to_string(beta) + ") must be between 0.0 and 1.0.");
    }

    std::mt19937 rng(seed);
    Matrix       adj(N, dVec(N, 0.0));

    // Initial ring lattice - create k outgoing edges for each node
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 1; j <= k; ++j)
        {
            size_t target  = (i + j) % N;
            adj[target][i] = weight;
        }
    }

    // Rewiring with improved logic
    std::bernoulli_distribution           rewire_dist(beta);
    std::uniform_int_distribution<size_t> node_dist(0, N - 1);

    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 1; j <= k; ++j)
        {
            size_t neighbor = (i + j) % N;
            if (rewire_dist(rng))
            {
                size_t attempts     = 0;
                size_t max_attempts = N * static_cast<size_t>(std::sqrt(N)); // Reasonable limit to prevent infinite loops
                bool   found        = false;
                size_t new_neighbor;

                while (!found && attempts < max_attempts)
                {
                    new_neighbor = node_dist(rng);
                    if (new_neighbor != i && 
                        new_neighbor != neighbor && 
                        adj[i][new_neighbor] == 0.0)
                    {
                        found = true;
                    }
                    attempts++;
                }

                if (found)
                {
                    adj[neighbor][i]     = 0.0;
                    adj[new_neighbor][i] = weight;
                }
            }
        }
    }
    return adj;
}

// Multilayered Matrix: block diagonal matrix, each block is a layer
inline Matrix multilayered(const std::vector<Matrix> &layers)
{
    size_t offset = 0;
    size_t total_N = 0;
    for (const auto &layer : layers) total_N += layer.size();
    Matrix adj(total_N, dVec(total_N, 0.0));
    for (const auto &layer : layers)
    {
        size_t n = layer.size();
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 0; j < n; ++j)
            {
                adj[offset + i][offset + j] = layer[i][j];
            }
        }
        offset += n;
    }
    return adj;
}

// Modular Matrix: nodes are divided into modules, dense within, sparse between
/*----------------------------------------------------------*/
// Parameters:
//   module_size: Number of nodes per module (default 100)
//   num_modules: Number of modules (default 10)
//   p_in: Probability of within-module connection (default 0.9)
//   p_out: Probability of between-module connection (default 0.1)
//   in_weight: Weight for within-module connections (default 1.0)
//   out_weight: Weight for between-module connections (default 1.0)
//   seed: Random seed (default 42)
// Returns:
//   Modular adjacency matrix (N x N), where N = module_size * num_modules
inline Matrix modular(
    size_t   module_size = 100,
    size_t   num_modules = 10,
    double   p_in        = 0.9,
    double   p_out       = 0.1,
    double   in_weight   = 1.0,
    double   out_weight  = 1.0,
    unsigned seed        = 42
)
{
    // Error handling
    if (module_size == 0)
    {
        throw std::invalid_argument("[modular] Module size cannot be zero.");
    }
    if (num_modules == 0)
    {
        throw std::invalid_argument("[modular] Number of modules cannot be zero.");
    }
    if (p_in < 0.0 || p_in > 1.0)
    {
        throw std::invalid_argument("[modular] Within-module connection probability p_in (" +
                                    std::to_string(p_in) + ") must be between 0.0 and 1.0.");
    }
    if (p_out < 0.0 || p_out > 1.0)
    {
        throw std::invalid_argument("[modular] Between-module connection probability p_out (" +
                                    std::to_string(p_out) + ") must be between 0.0 and 1.0.");
    }
    // Calculate total network size
    size_t N = module_size * num_modules;
    if (N < num_modules)
    { // Check for overflow
        throw std::invalid_argument("[modular] Total network size (module_size * num_modules) exceeds maximum size_t value.");
    }

    Matrix                      adj(N, dVec(N, 0.0));
    std::mt19937                rng(seed);
    std::bernoulli_distribution in_dist(p_in), out_dist(p_out);
    for (size_t m = 0; m < num_modules; ++m)
    {
        size_t start = m * module_size;
        size_t end   = (m == num_modules - 1) ? N : (m + 1) * module_size;
        // Within module
        for (size_t i = start; i < end; ++i)
        {
            for (size_t j = start; j < end; ++j)
            {
                if (i != j && in_dist(rng))
                {
                    adj[i][j] = in_weight;
                }
            }
        }
        // Between modules
        for (size_t n = m + 1; n < num_modules; ++n)
        {
            size_t n_start = n * module_size;
            size_t n_end   = (n == num_modules - 1) ? N : (n + 1) * module_size;
            for (size_t i = start; i < end; ++i)
            {
                for (size_t j = n_start; j < n_end; ++j)
                {
                    if (out_dist(rng))
                    {
                        adj[i][j] = out_weight;
                        adj[j][i] = out_weight;
                    }
                }
            }
        }
    }
    return adj;
}

// Hierarchical Matrix: recursively nested modules with level-dependent connection probabilities
/*----------------------------------------------------------*/
// Parameters:
//   N: Base module size at the lowest level (default 100)
//   levels: Number of hierarchical levels (0 returns Erdos-Renyi N x N matrix) (default 2)
//   p_in: Base within-module connection probability (decays with level) (default 0.9)
//   p_out: Base between-module connection probability (decays with level) (default 0.1)
//   in_weight: Weight for within-module connections (default 1.0)
//   out_weight: Weight for between-module connections (default 1.0)
//   level_decay: Decay factor for connection probabilities at each level (0 < level_decay <= 1) (default 0.5)
//   seed: Random seed (default 42)
//   base_module_num: Number of modules in the base level of the hierarchy (default 2)
// Returns:
//   Hierarchical adjacency matrix (size: N * 2^(levels-1) * base_module_num)
inline Matrix hierarchical(
    size_t   N               = 100,
    size_t   levels          = 2,
    double   p_in            = 0.9,
    double   p_out           = 0.1,
    double   in_weight       = 1.0,
    double   out_weight      = 1.0,
    double   level_decay     = 0.5,
    unsigned seed            = 42,
    size_t   base_module_num = 2
)
{
    // Parameter validation and base cases
    if (N == 0)
    {
        throw std::invalid_argument("[hierarchical] Base module size N cannot be zero.");
    }
    if (p_in < 0.0 || p_in > 1.0)
    {
        throw std::invalid_argument("[hierarchical] Within-module probability p_in (" +
                                    std::to_string(p_in) + ") must be between 0.0 and 1.0.");
    }
    if (p_out < 0.0 || p_out > 1.0)
    {
        throw std::invalid_argument("[hierarchical] Between-module probability p_out (" +
                                    std::to_string(p_out) + ") must be between 0.0 and 1.0.");
    }
    if (level_decay <= 0.0 || level_decay > 1.0)
    {
        throw std::invalid_argument("[hierarchical] Level decay factor (" +
                                    std::to_string(level_decay) + ") must be between 0.0 and 1.0.");
    }
    // Check for potential size overflow
    size_t max_nodes = N * (static_cast<size_t>(1) << (levels - 1)) * base_module_num;
    if (max_nodes / (N * base_module_num) != (static_cast<size_t>(1) << (levels - 1)))
    {
        throw std::invalid_argument("[hierarchical] Network size (N * 2^levels) exceeds maximum size_t value.");
    }
    if (base_module_num == 0)
    {
        throw std::invalid_argument("[hierarchical] base_module_num cannot be zero.");
    }
    if (levels == 0)
    {
        return erdos_renyi(N, p_in, in_weight, out_weight, seed);
    }
    double this_p_in  = p_in * std::pow(level_decay, levels - 1);
    double this_p_out = p_out * std::pow(level_decay, levels - 1);
    if (levels == 1)
    {
        return modular(N, base_module_num, this_p_in, this_p_out, in_weight, out_weight, seed);
    }
    Matrix adj(max_nodes, dVec(max_nodes, 0.0));
    // Recursively build left and right submodules
    auto   left  = hierarchical(N, levels - 1, p_in, p_out, in_weight, out_weight, level_decay, seed + 1, base_module_num);
    auto   right = hierarchical(N, levels - 1, p_in, p_out, in_weight, out_weight, level_decay, seed + 2, base_module_num);
    size_t half  = N * std::pow(2, levels - 2) * base_module_num;
    for (size_t i = 0; i < half; ++i)
    {
        for (size_t j = 0; j < half; ++j)
        {
            adj[i][j]               = left[i][j];
            adj[half + i][half + j] = right[i][j];
        }
    }
    std::mt19937                rng(seed);
    std::bernoulli_distribution out_dist(this_p_out);
    for (size_t i = 0; i < half; ++i)
    {
        for (size_t j = half; j < max_nodes; ++j)
        {
            if (out_dist(rng))
            {
                adj[i][j] = out_weight;
                adj[j][i] = out_weight;
            }
        }
    }
    return adj;
}

// pick_topology: Create a network topology based on type and parameters
/*----------------------------------------------------------*/
// Parameters:
//   topology_type: Type of network to create ("random", "erdos_renyi", "small_world", "modular", "hierarchical") (default "random")
//   N: Number of nodes interpreted per topology type (default 100)
//   a, b, c, d, e, f: Parameters whose meaning depends on topology_type (see below)
//   seed: Random seed (default 42)
//   base_module_num: Number of modules in the base level of the hierarchy (default 2)
// Returns:
//   Adjacency matrix for the selected topology
/*----------------------------------------------------------*/
// Note: Error-handling is performed at this high level for us  er-friendliness, so that invalid parameters are caught early with clear messages.
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
)
{
    // Basic N=0 check, specific functions might have more detailed checks.
    // Some types like multilayered might derive N differently, so not checking N for all types here.
    if (N == 0)
    {
        throw std::invalid_argument("[pick_topology] Number of nodes N cannot be zero for type '" + topology_type + "'.");
    }

    if (topology_type == "random")
    {   // a = min_weight, b = max_weight, c is not used
        return random(N, a, b, seed);
    }
    else if (topology_type == "erdos_renyi")
    {   // a = p, b = min_weight, c = max_weight
        return erdos_renyi(N, a, b, c, seed);
    }
    else if (topology_type == "erdos_renyi_uniform")
    {   // a = p, b = min_weight, c = max_weight
        return erdos_renyi_uniform(N, a, b, c, seed);
    }
    else if (topology_type == "erdos_renyi_symmetric")
    {   // a = p, b = min_weight, c = max_weight
        return erdos_renyi_symmetric(N, a, b, c, seed);
    }
    else if (topology_type == "erdos_renyi_symmetric_uniform")
    {   // a = p, b = min_weight, c = max_weight
        return erdos_renyi_symmetric_uniform(N, a, b, c, seed);
    }    
    else if (topology_type == "small_world")
    {   // a = k (num_neighbors_half for ring), b = beta, c = weight
        return small_world(N, static_cast<size_t>(a), b, c, seed);
    }
    else if (topology_type == "small_world_directed")
    {   // a = k (num_neighbors_half for ring), b = beta, c = weight
        return small_world_directed(N, static_cast<size_t>(a), b, c, seed);
    }
    else if (topology_type == "modular")
    {   // Assuming: a = num_modules, b = p_in, c = p_out
        // d = in_weight, e = out_weight, and N is the base module size
        size_t num_modules_val = static_cast<size_t>(a);
        if (num_modules_val == 0)
        {
            throw std::invalid_argument("[pick_topology] For 'modular', num_modules (parameter 'a') cannot be zero if N > 0.");
        }
        return modular(N, num_modules_val, b, c, d, e, seed);
    }
    else if (topology_type == "hierarchical")
    {   // Assuming: a = levels, b = p_in, c = p_out
        // d = in_weight, e = out_weight, f = level_decay, and N is the base module size
        size_t levels_val = static_cast<size_t>(a);
        if (levels_val == 0)
        {
            levels_val = 2;
            std::cout  << "[pick_topology] Warning: levels (parameter 'a') is 0, setting to 2 (default for hierarchical)." << std::endl;
        }
        // hierarchical itself handles levels_val == 0 by returning an empty N x N matrix.
        if (base_module_num == 0)
        {
            throw std::invalid_argument("[pick_topology] For 'hierarchical', base_module_num (parameter 'f') cannot be zero.");
        }
        return hierarchical(N, levels_val, b, c, d, e, f, seed, base_module_num);
    }
    else
    {
        throw std::invalid_argument("[pick_topology] Unknown topology type: '" + topology_type + "'.");
    }
}

// effective_multiplex: Weighted sum of multiple network layers
/*----------------------------------------------------------*/
// Parameters:
//   layers: Vector of adjacency matrices (all must be N x N)
//   layer_weights: Vector of weights (one per layer)
// Returns:
//   Weighted sum of layers as a single adjacency matrix
inline Matrix effective_multiplex(
    const std::vector<Matrix> &layers,
    const dVec                 &layer_weights
)
{
    if (layers.empty())
    {
        throw std::invalid_argument("Input 'layers' cannot be empty.");
    }
    if (layers.size() != layer_weights.size())
    {
        throw std::invalid_argument("Number of layers must match number of layer weights.");
    }

    const size_t num_layers = layers.size();

    const size_t N          = layers[0].size();
    if (N == 0)
    {
        throw std::invalid_argument("Dimension N of layers (number of nodes) cannot be zero.");
    }

    for (size_t l = 0; l < num_layers; ++l)
    {
        if (layers[l].size() != N)
        {
            throw std::invalid_argument("All layers must have the same dimension N (number of nodes). Layer " + std::to_string(l) + " has " + std::to_string(layers[l].size()) + " nodes, expected " + std::to_string(N) + ".");
        }
        for (size_t i = 0; i < N; ++i)
        {
            if (layers[l][i].size() != N)
            {
                throw std::invalid_argument("All layers must be N x N square matrices. Layer " + std::to_string(l) + ", row " + std::to_string(i) + " has " + std::to_string(layers[l][i].size()) + " columns, expected " + std::to_string(N) + ".");
            }
        }
    }

    Matrix adj_out(N, dVec(N, 0.0));
    for (size_t l = 0; l < num_layers; ++l)
    {
        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                adj_out[i][j] += layer_weights[l] * layers[l][i][j];
            }
        }
    }
    return adj_out;
}

// generate_multiplex_network_layers: Generate multiple network layers for a multiplex network
/*----------------------------------------------------------*/
// Parameters:
//   N: Number of nodes (consistent across all layers) (default 100)
//   num_layers: Number of layers to generate (default 2)
//   layer_generation_types: Vector of topology types for each layer (default "erdos_renyi")
//   layer_params: Matrix of parameters for each layer (one row per layer)
//   layer_seeds: Vector of random seeds for each layer (default 42)
// Returns:
//   Vector of adjacency matrices, one per layer
inline std::vector<Matrix> generate_multiplex_network_layers(
    size_t N,
    size_t num_layers,
    const  std::vector<std::string> &layer_generation_types,
    const  Matrix                   &layer_params,
    const  std::vector<unsigned>    &layer_seeds
)
{
    // Basic parameter validation
    if (num_layers == 0)
    {
        throw std::invalid_argument("[generate_multiplex_network_layers] Number of layers cannot be zero.");
    }
    if (N == 0)
    {
        throw std::invalid_argument("[generate_multiplex_network_layers] Number of nodes N cannot be zero.");
    }

    // Check consistency of input vectors/matrix
    if (layer_generation_types.size() != num_layers)
    {
        throw std::invalid_argument("[generate_multiplex_network_layers] Number of layer types (" +
                                    std::to_string(layer_generation_types.size()) + ") does not match num_layers (" +
                                    std::to_string(num_layers) + ").");
    }
    if (layer_params.size() != num_layers)
    {
        throw std::invalid_argument("[generate_multiplex_network_layers] Number of parameter sets (" +
                                    std::to_string(layer_params.size()) + ") does not match num_layers (" +
                                    std::to_string(num_layers) + ").");
    }
    if (layer_seeds.size() != num_layers)
    {
        throw std::invalid_argument("[generate_multiplex_network_layers] Number of seeds (" +
                                    std::to_string(layer_seeds.size()) + ") does not match num_layers (" +
                                    std::to_string(num_layers) + ").");
    }

    // Validate layer types and parameter counts
    for (size_t i = 0; i < num_layers; ++i)
    {
        const auto &type = layer_generation_types[i];
        const auto &params = layer_params[i];

        // Check if topology type is valid
        if (type != "random" && type != "erdos_renyi" && type != "small_world" &&
            type != "modular" && type != "hierarchical" && 
            type != "erdos_renyi_uniform" && type != "erdos_renyi_symmetric"  && 
            type != "erdos_renyi_symmetric_uniform" && type != "small_world_directed")
        {
            throw std::invalid_argument("[generate_multiplex_network_layers] Invalid topology type '" +
                                        type + "' for layer " + std::to_string(i) + ".");
        }

        // Check minimum parameter count for each type
        size_t min_params = 3; // Most types need at least 3 parameters
        if (type == "modular")
            min_params = 5; // Modular needs 5 parameters
        if (type == "hierarchical")
            min_params = 7; // Hierarchical needs 6 parameters

        if (params.size() < min_params)
        {
            throw std::invalid_argument("[generate_multiplex_network_layers] Insufficient parameters for " +
                                        type + " topology in layer " + std::to_string(i) + ". Expected at least " +
                                        std::to_string(min_params) + ", got " + std::to_string(params.size()) + ".");
        }
    }

    std::vector<Matrix> multiplex_layers(num_layers);

    for (size_t l = 0; l < num_layers; ++l)
    {
        const std::string &type                     = layer_generation_types[l];
        const dVec        &current_layer_gen_params = layer_params[l]; // Renamed for clarity
        unsigned          seed                      = layer_seeds[l];
        Matrix            current_layer_adj;

        if (type == "erdos_renyi")
        { // Expected params: p, min_weight, max_weight
            current_layer_adj = erdos_renyi(N, current_layer_gen_params[0], current_layer_gen_params[1], current_layer_gen_params[2], seed);
        }
        else if (type == "small_world")
        { // Expected params: k, beta, weight
            current_layer_adj = small_world(N, static_cast<size_t>(current_layer_gen_params[0]), current_layer_gen_params[1], current_layer_gen_params[2], seed);
        }
        else if (type == "random")
        {
            // Expected params: min_weight, max_weight
            current_layer_adj = random(N, current_layer_gen_params[0], current_layer_gen_params[1], seed);
        }
        else if (type == "modular")
        { // Expected params: num_modules, p_in, p_out, in_weight, out_weight
            current_layer_adj = modular(N, static_cast<size_t>(current_layer_gen_params[0]), current_layer_gen_params[1], current_layer_gen_params[2], current_layer_gen_params[3], current_layer_gen_params[4], seed);
        }
        else if (type == "hierarchical")
        { // Expected params: levels, p_in, p_out, in_weight, out_weight, level_decay
            current_layer_adj = hierarchical(N, static_cast<size_t>(current_layer_gen_params[0]), current_layer_gen_params[1], current_layer_gen_params[2], current_layer_gen_params[3], current_layer_gen_params[4], current_layer_gen_params[5], seed, static_cast<size_t>(current_layer_gen_params[6]));
        }
        multiplex_layers[l] = current_layer_adj;
    }
    return multiplex_layers;
}

// -----------------------------------------------------------------------------
// Calculate in-degrees of a given adjacency matrix
// Parameters:
//   adj: Dense adjacency matrix (N x N)
//   threshold: Minimum absolute value to consider as an edge (default 1e-12)
// Returns:
//   Vector of in-degrees (number of incoming edges for each node)
inline wVec in_degrees(const Matrix& adj, double threshold = 1e-12)
{
    size_t N = adj.size();
    wVec   indeg(N, 0);
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            if (std::abs(adj[i][j]) > threshold)
                ++indeg[i];
        }
    }
    return indeg;
}

// -----------------------------------------------------------------------------
// Calculate out-degrees of a given adjacency matrix
// Parameters:
//   adj: Dense adjacency matrix (N x N)
//   threshold: Minimum absolute value to consider as an edge (default 1e-12)
// Returns:
//   Vector of out-degrees (number of outgoing edges for each node)
inline wVec out_degrees(const Matrix& adj, double threshold = 1e-12)
{
    size_t N = adj.size();
    wVec   outdeg(N, 0);
    for (size_t i = 0; i < N; ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            if (std::abs(adj[i][j]) > threshold)
                ++outdeg[i];
        }
    }
    return outdeg;
}

// -----------------------------------------------------------------------------
// Calculate in-degrees of a given sparse adjacency matrix
// Parameters:
//   sparse: Sparse adjacency matrix (N x N)
//   threshold: Minimum absolute value to consider as an edge (default 1e-12)
// Returns:
//   Vector of in-degrees (number of incoming edges for each node)
inline wVec in_degrees(const SparseMatrix& sparse, double threshold = 1e-12)
{
    size_t N = sparse.rows.size();
    wVec   indeg(N, 0);
    for (size_t i = 0; i < N; ++i)
    {
        for (const auto& entry : sparse.rows[i])
        {
            size_t j     = entry.first;
            double value = entry.second;
            if (std::abs(value) > threshold)
                ++indeg[i];
        }
    }
    return indeg;
}

// -----------------------------------------------------------------------------
// Calculate out-degrees of a given sparse adjacency matrix
// Parameters:
//   sparse: Sparse adjacency matrix (N x N)
//   threshold: Minimum absolute value to consider as an edge (default 1e-12)
// Returns:
//   Vector of out-degrees (number of outgoing edges for each node)
inline wVec out_degrees(const SparseMatrix& sparse, double threshold = 1e-12)
{
    size_t N = sparse.rows.size();
    wVec   outdeg(N, 0);
    for (size_t i = 0; i < N; ++i)
    {
        for (const auto& entry : sparse.rows[i])
        {
            size_t j     = entry.first;
            double value = entry.second;
            if (std::abs(value) > threshold)
                ++outdeg[j];
        }
    }
    return outdeg;
}
 