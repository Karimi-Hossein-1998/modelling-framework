#pragma once
#include "../../typedefs/header.hpp"

// Special modular Kuramoto model (ordinary)
inline dVec  kuramoto_special_modular(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              intra_K,
    double              inter_K,
    double              alpha,
    const size_t        module_size
)
{
    const size_t N = theta.size(); // N is derived from theta, so N > 0 here.

    std::vector<std::jthread> threads;
    auto   intra_k     = intra_K / static_cast<double>( N ); // Renamed k = K/N to k_norm for clarity
    auto   inter_k     = inter_K / static_cast<double>( N ); // Renamed k = K/N to k_norm for clarity
    dVec    dtheta_dt   = dVec( N, 0.0 );

    for (size_t i = 0; i < N; ++i) 
    {
        size_t module_i = i / module_size;
        double sum = 0.0;
        for (size_t j = 0; j < N; ++j) 
        {
            size_t module_j = j / module_size;
            if ( i != j )
            {
                if ( module_i == module_j )
                {
                    sum += intra_k * std::sin(theta[j] - theta[i] - alpha);
                }
                else
                {
                    sum += inter_k * std::sin(theta[j] - theta[i] - alpha);
                }
            }
        }
        dtheta_dt[i] = omega[i] + sum;
    }
    return dtheta_dt;
}

// Special modular Kuramoto model (parallel)
inline dVec  kuramoto_special_modular_parallel(
    double              time,
    const dVec&         theta,
    const dVec&         omega,
    double              intra_K,
    double              inter_K,
    double              alpha,
    const size_t        module_size
)
{
    const size_t N = theta.size(); // N is derived from theta, so N > 0 here.

    std::vector<std::jthread> threads;
    size_t num_threads = std::min(N, static_cast<size_t>(std::max(1u, std::thread::hardware_concurrency())));
    if ( num_threads == 0 ) num_threads = 1;
    auto   intra_k     = intra_K / static_cast<double>( N ); // Renamed k = K/N to k_norm for clarity
    auto   inter_k     = inter_K / static_cast<double>( N ); // Renamed k = K/N to k_norm for clarity
    dVec    dtheta_dt   = dVec( N, 0.0 );
    size_t chunk_size  = N / num_threads;

    for (size_t thread_index = 0; thread_index < num_threads; ++thread_index) 
    {
        threads.emplace_back([&, thread_index]() 
        {
            size_t start = thread_index * chunk_size;
            size_t end   = (thread_index == num_threads - 1) ? N : (thread_index + 1) * chunk_size;
            
            for (size_t i = start; i < end; ++i) 
            {
                size_t module_i = i / module_size;
                double sum = 0.0;
                for (size_t j = 0; j < N; ++j) 
                {
                    size_t module_j = j / module_size;
                    if ( i != j )
                    {
                        if ( module_i == module_j )
                        {
                            sum += intra_k * std::sin(theta[j] - theta[i] - alpha);
                        }
                        else
                        {
                            sum += inter_k * std::sin(theta[j] - theta[i] - alpha);
                        }
                    }
                }
                dtheta_dt[i] = omega[i] + sum;
            }
        });
    }
    // No need to join, jthread automatically joins in destructor
    return dtheta_dt;
} 