#pragma once
#include "../../../../utility/dde_utility.hpp"
#include "../../../../typedefs/header.hpp"

// General DDE version of Euler's method (RK1) for multiple delays
// Solves delay differential equations of the form:
// dy/dt = f(t, y(t), y(t - τ₁), y(t - τ₂), ..., y(t - τₙ))
// where τ_i are the delay times
inline SolverResults rk1_dde_general_solver(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f  = Params.delay_derivative_multi;  // DDE function f(t, y, y_delayed_states)
    const auto&  y0 = Params.initial_conditions;
    const double t0 = Params.t0;
    const double t1 = Params.t1;
    const double dt = Params.dt;
    const size_t N  = y0.size();
    
    // DDE-specific parameters
    const auto& delay_times      = Params.delay_times;        // Vector of delay times
    const auto& history_function = Params.history_function; // Pre-history values
    
    // Validate DDE parameters
    if (!f)
    {
        throw std::invalid_argument("[rk1_dde_general_solver] DDE derivative function is required");
    }
    if (delay_times.empty())
    {
        throw std::invalid_argument("[rk1_dde_general_solver] At least one delay time is required for DDE");
    }
    
    // Find maximum delay for pre-allocation
    double max_delay = *std::max_element(delay_times.begin(), delay_times.end());
    
    // Initialize solution storage
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = Matrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1);
    
    solution[0]    = y0;
    time_points[0] = t0;
    
    auto y = y0;
    dVec  k1(N, 0.0);
    Matrix y_delayed_states(delay_times.size(), dVec(N, 0.0));
    
    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];
        y_delayed_states = interpolate_multiple_delays(delay_times, t, time_points, solution, t0, history_function);
        k1 = f(t, y, y_delayed_states);
        for (size_t j = 0; j < N; ++j)
            y[j] += dt * k1[j];
        solution[i + 1]    = y;
        time_points[i + 1] = t + dt;
    }
    
    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

// Basic interface wrapper for rk1_dde_general_solver
inline Matrix rk1_dde_general_solver(
    DelayFuncMulti deriv,
    const dVec&    y0,
    double         t0,
    double         t1,
    double         dt,
    const dVec&    delay_times      = {1.0},
    const Func&    history_function = {}
)
{
    auto params                    = SolverParameters{};
    params.delay_derivative_multi  = deriv;
    params.initial_conditions      = y0;
    params.t0                      = t0;
    params.t1                      = t1;
    params.dt                      = dt;
    params.delay_times             = delay_times;
    params.history_function        = history_function;
    params.is_dde                  = true;
    return rk1_dde_general_solver(params).solution;
} 