# The Kuramoto Model

## 1. General Kuramoto Model

### 1.1 Basic Model

The Kuramoto model describes the dynamics of a system of coupled oscillators. For $N$ oscillators, each with phase $\theta_i$ and natural frequency $\omega_i$:

$$
\frac{d\theta_i}{dt} = \omega_i + \frac{K}{N} \sum_{j=1}^N A_{ij} \sin(\theta_j - \theta_i - \alpha), \quad i = 1,\ldots,N
$$

where:

- $\theta_i(t)$ is the phase of oscillator $i$ at time $t$
- $\omega_i$ is the natural frequency of oscillator $i$
- $K$ is the coupling strength
- $A_{ij}$ is the adjacency matrix entry (connection weight from $j$ to $i$)
- $\alpha$ is the phase lag parameter
- $N$ is the number of oscillators

### 1.2 Implementation Variants

1. **Dense Network**:
   - Full adjacency matrix $A_{ij}$
   - Direct summation over all $j$

2. **Sparse Network**:
   $$
   \frac{d\theta_i}{dt} = \omega_i + \frac{K}{N} \sum_{j \in \mathcal{N}_i} w_{ij} \sin(\theta_j - \theta_i - \alpha)
   $$
   where $\mathcal{N}_i$ is the set of neighbors of node $i$ and $w_{ij}$ are the edge weights.

3. **Parallel Implementation**:
   - System divided into chunks of size $\lceil N/n_{\text{threads}} \rceil$
   - Each thread processes a subset of oscillators independently

## 2. Special Modular Kuramoto Model

### 2.1 Model Definition

For a system divided into modules of size $M$:

$$
\frac{d\theta_i}{dt} = \omega_i + \sum_{j \neq i} K_{ij} \sin(\theta_j - \theta_i - \alpha)
$$

where:
$$
K_{ij} = \begin{cases}
\frac{K_{\text{intra}}}{N} & \text{if } i,j \text{ in same module} \\
\frac{K_{\text{inter}}}{N} & \text{if } i,j \text{ in different modules}
\end{cases}
$$

Parameters:

- $K_{\text{intra}}$: Coupling strength within modules
- $K_{\text{inter}}$: Coupling strength between modules
- Module assignment: $m_i = \lfloor i/M \rfloor$

### 2.2 Order Parameters

#### 2.2.1 Global Order Parameter

The global synchronization is measured by:

$$
r(t)e^{i\psi(t)} = \frac{1}{N}\sum_{j=1}^N e^{i\theta_j(t)}
$$

where:

- $r(t)$ is the magnitude of synchronization
- $\psi(t)$ is the average phase
- $0 \leq r \leq 1$

#### 2.2.2 Local Order Parameters

For each module $m$:

$$
r_m(t)e^{i\psi_m(t)} = \frac{1}{M}\sum_{j \in \mathcal{M}_m} e^{i\theta_j(t)}
$$

where $\mathcal{M}_m$ is the set of oscillators in module $m$.

## 3. Numerical Considerations

### 3.1 Stability and Accuracy

1. **Time Step Selection**:
   - Should be smaller than fastest oscillation period
   - Typical choice: $\Delta t < \frac{2\pi}{\max_i|\omega_i|}$

2. **Phase Wrapping**:
   - Phases should be kept in $[-\pi, \pi]$ or $[0, 2\pi]$
   - Implement periodic boundary conditions

### 3.2 Performance Optimization

1. **Parallel Implementation**:
   - Thread count: $\min(N, n_{\text{cores}})$
   - Load balancing through even chunk distribution

2. **Sparse Networks**:
   - Only compute contributions from connected nodes
   - Use sparse matrix format when density $< 0.5$

3. **Vectorization**:
   - Align memory for SIMD operations
   - Group similar computations for efficient vectorization

## 4. Implementation Notes

### 4.1 Key Functions

1. **General Kuramoto**:

   ```cpp
   dVec kuramoto_general(
       double time,
       const dVec& theta,
       const dVec& omega,
       double K,
       const Matrix& adj,
       double alpha
   );
   ```

2. **Sparse Kuramoto**:

   ```cpp
   dVec kuramoto_sparse(
       double time,
       const dVec& theta,
       const dVec& omega,
       double K,
       const SparseMatrix& sparse_adj,
       double alpha
   );
   ```

3. **Modular Kuramoto**:

   ```cpp
   dVec kuramoto_special_modular(
       double time,
       const dVec& theta,
       const dVec& omega,
       double intra_K,
       double inter_K,
       double alpha,
       const size_t module_size
   );
   ```

### 4.2 Error Handling

Critical validation checks:

1. $N > 0$: System size must be positive
2. Dimensions must match: $|\theta| = |\omega| = N$
3. Adjacency matrix must be $N \times N$
4. Module size must divide $N$ evenly in modular case
