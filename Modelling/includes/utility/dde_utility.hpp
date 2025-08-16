#pragma once
#include "../typedefs/header.hpp"

// Delay interpolation utility for DDE solvers
// Interpolates the solution at time t - tau using linear interpolation
inline dVec  interpolate_delay(
    double        t_delayed,             // Time to interpolate at (t - tau)
    const dVec&   time_points,           // Available time points
    const Matrix& solution,              // Solution at those time points
    double        t0,                    // Initial time
    const Func&   history_function       // Pre-history function values
)
{
    const size_t N = solution[0].size();
    dVec  y_delayed(N, 0.0);

    // If delayed time is before t0, use history function
    if (t_delayed < t0)
    {
        // Simple constant extrapolation from history
        // In a more sophisticated implementation, you might want to interpolate
        // between multiple history points or use a more complex history function
        return history_function(t_delayed);
    }
    
    // Find the two time points that bracket t_delayed
    size_t idx = 0;
    for (size_t i = 0; i < time_points.size() - 1; ++i)
    {
        if (time_points[i] <= t_delayed && t_delayed <= time_points[i + 1])
        {
            idx = i;
            break;
        }
    }
    
    // Linear interpolation
    const double t1    = time_points[idx];
    const double t2    = time_points[idx + 1];
    const double alpha = (t_delayed - t1) / (t2 - t1);
    
    for (size_t i = 0; i < N; ++i)
    {
        y_delayed[i] = (1.0 - alpha) * solution[idx][i] + alpha * solution[idx + 1][i];
    }
    
    return y_delayed;
}

// Multiple delay interpolation utility for general DDE solvers
// Interpolates the solution at multiple delayed times t - τ_i
// Returns a matrix where each row i contains the state at time t - τ_i
inline Matrix interpolate_multiple_delays(
    const dVec&   delay_times,           // Vector of delay times τ_i
    double        t,                     // Current time
    const dVec&   time_points,           // Available time points
    const Matrix& solution,              // Solution at those time points
    double        t0,                    // Initial time
    const Func&   history_function       // Pre-history function values
)
{
    const size_t num_delays = delay_times.size();
    const size_t N = solution[0].size();
    Matrix y_delayed_states(num_delays, dVec(N, 0.0));
    
    for (size_t delay_idx = 0; delay_idx < num_delays; ++delay_idx)
    {
        double t_delayed = t - delay_times[delay_idx];
        y_delayed_states[delay_idx] = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
    }
    
    return y_delayed_states;
}