# Network Theory and Topology

## 1. Basic Definitions

### 1.1 Network Representation

A network (or graph) is represented by an adjacency matrix $A \in \mathbb{R}^{N \times N}$ where:

- $N$ is the number of nodes
- $A_{ij}$ represents the connection weight from node $j$ to node $i$
- $A_{ij} = 0$ indicates no connection
- $A_{ij} \neq 0$ indicates a connection with weight $A_{ij}$

For a sparse representation, only non-zero elements are stored as pairs $(j, A_{ij})$ for each row $i$.

### 1.2 Network Properties

For a given network with adjacency matrix $A$:

1. **In-degree** of node $i$:
   $$
   k_i^{\text{in}} = \sum_{j=1}^N \mathbb{1}(|A_{ij}| > \theta)
   $$
   where $\theta$ is a threshold (typically $10^{-12}$)

2. **Out-degree** of node $j$:
   $$
   k_j^{\text{out}} = \sum_{i=1}^N \mathbb{1}(|A_{ij}| > \theta)
   $$

3. **Density** of the network:
   $$
   \rho = \frac{\sum_{i=1}^N \sum_{j \neq i} \mathbb{1}(|A_{ij}| > \theta)}{N(N-1)}
   $$

## 2. Network Topologies

### 2.1 Random Networks

A random network with $N$ nodes is generated with weights:

$$
A_{ij} = \begin{cases}
w_{ij} \sim U(w_{\text{min}}, w_{\text{max}}) & \text{if } i \neq j \\
0 & \text{if } i = j
\end{cases}
$$

where $U(a,b)$ denotes uniform distribution on $[a,b]$.

### 2.2 Erdős-Rényi Networks

For $N$ nodes and connection probability $p$:

$$
A_{ij} = \begin{cases}
w_{ij} \sim U(w_{\text{min}}, w_{\text{max}}) & \text{with probability } p \text{ if } i \neq j \\
0 & \text{otherwise}
\end{cases}
$$

Variants include:

- Uniform: Exactly $N(N-1)p$ edges
- Symmetric: $A_{ij} = A_{ji}$
- Symmetric Uniform: Exactly $\frac{N(N-1)p}{2}$ undirected edges

### 2.3 Small-World Networks (Watts-Strogatz Model)

1. **Initial Ring Lattice**:
   - Each node connected to $k$ nearest neighbors
   - $A_{ij} = w$ if $\min(|i-j|, N-|i-j|) \leq \frac{k}{2}$

2. **Rewiring Process** with probability $\beta$:
   - For each edge $(i,j)$:
     - With probability $\beta$, rewire to random node $m$
     - Preserve weight: $A_{im} = w$
     - Ensure no self-loops or duplicate edges

### 2.4 Modular Networks

For $M$ modules of size $s$ each ($N = Ms$):

$$
A_{ij} = \begin{cases}
w_{\text{in}} & \text{with prob. } p_{\text{in}} \text{ if in same module} \\
w_{\text{out}} & \text{with prob. } p_{\text{out}} \text{ if in different modules} \\
0 & \text{otherwise}
\end{cases}
$$

### 2.5 Hierarchical Networks

For $L$ levels with base module size $N$:

1. **Connection Probabilities**:
   $$p_l = p_{\text{base}} \cdot \alpha^{l-1}$$
   where $\alpha$ is the level decay factor

2. **Structure**:
   - Level 1: Base modules of size $N$
   - Level $l$: Combines two level $(l-1)$ structures
   - Total size: $N \cdot 2^{L-1}$

### 2.6 Multiplex Networks

A multiplex network consists of $M$ layers sharing the same nodes:

1. **Layer Generation**:
   - Each layer $A^{(l)}$ generated independently
   - Can use different topologies per layer

2. **Effective Combined Network**:
   $$
   A_{\text{eff}} = \sum_{l=1}^M w_l A^{(l)}
   $$
   where $w_l$ are layer weights

## 3. Implementation Details

### 3.1 Sparse Conversion

Convert dense to sparse when:
$$\rho < \rho_{\text{threshold}}$$
typically with $\rho_{\text{threshold}} = 0.5$

### 3.2 Matrix Operations

1. **Density Calculation**:
   $$
   \rho = \frac{|\{(i,j) : |A_{ij}| > \theta, i \neq j\}|}{N(N-1)}
   $$

2. **Weight Distribution**:
   - Uniform on $[w_{\text{min}}, w_{\text{max}}]$ for random weights
   - Constant weights for specific topologies (e.g., small-world)

### 3.3 Error Handling

Key validation criteria:

1. $N > 0$: Number of nodes must be positive
2. $0 \leq p \leq 1$: Probabilities must be in $[0,1]$
3. $w_{\text{min}} \leq w_{\text{max}}$: Weight bounds must be ordered
4. $k < N$: Number of neighbors in small-world must be less than total nodes
