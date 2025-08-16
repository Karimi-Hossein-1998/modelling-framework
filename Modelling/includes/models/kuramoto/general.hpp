#pragma once
#include "../../typedefs/header.hpp"

// General Kuramoto model with phase-lag and flexible Matrix
inline dVec  kuramoto_general(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              K,
    const Matrix&       adj,
    double              alpha
)
{
    const size_t N = theta.size(); // N is derived from theta, so N > 0 here.

    // Original logic
    dVec  dtheta_dt = dVec( N, 0.0 );
    auto  k_norm    = K / static_cast<double>( N ); // Renamed k to k_norm for clarity
    for ( size_t i = 0; i < N; ++i )
    {
        double sum = 0.0;
        for ( size_t j = 0; j < N; ++j )
        {
            sum += adj[i][j] * std::sin( theta[j] - theta[i] - alpha );
        }
        dtheta_dt[i] = omega[i] + k_norm * sum;
    }
    return dtheta_dt;
}

// General Kuramoto model with phase-lag and flexible Matrix (parallel)
inline dVec  kuramoto_general_parallel(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              K,
    const Matrix&       adj,
    double              alpha
)
{
    const size_t N = theta.size(); // N is derived from theta, so N > 0 here.

    std::vector<std::jthread> threads;
    size_t num_threads = std::min(N, static_cast<size_t>(std::max(1u, std::thread::hardware_concurrency())));
    if ( num_threads == 0 ) num_threads = 1;
    size_t chunk_size  = N / num_threads;
    auto   k_norm      = K / static_cast<double>( N ); // Renamed k = K/N to k_norm for clarity
    dVec   dtheta_dt   = dVec( N, 0.0 );
    
    for (size_t t = 0; t < num_threads; ++t) 
    {
        threads.emplace_back([&, t]() 
        {
            size_t start = t * chunk_size;
            size_t end = (t == num_threads - 1) ? N : (t + 1) * chunk_size;
            
            for (size_t i = start; i < end; ++i) 
            {
                double sum = 0.0;
                for (size_t j = 0; j < N; ++j) 
                {
                    sum += adj[i][j] * std::sin(theta[j] - theta[i] - alpha);
                }
                dtheta_dt[i] = omega[i] + k_norm * sum;
            }
        });
    }
    // No need to join, jthread automatically joins in destructor
    return dtheta_dt;
} 