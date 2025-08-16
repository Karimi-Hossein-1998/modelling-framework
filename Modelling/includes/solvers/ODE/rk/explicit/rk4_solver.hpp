#pragma once
#include "../../../../typedefs/header.hpp"

// Classical fourth-order Runge-Kutta method (RK4)
// A widely used, well-balanced method for numerical integration, known for its
// accuracy and stability. It requires four derivative evaluations per step.
inline SolverResults rk4_solver(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f    = Params.derivative;
    const auto&  y0   = Params.initial_conditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();

    // Initialize solution storage
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = Matrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1);

    solution[0]    = y0;
    time_points[0] = t0;

    auto y        = y0;
    auto y_temp   = dVec(N, 0.0);
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    auto dt_half  = dt * 0.5;
    auto dt_third = dt / 3.0;
    auto dt_sixth = dt_third * 0.5;

    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];

        k1 = f(t, y);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_half * k1[j];
        k2 = f(t + dt_half, y_temp);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_half * k2[j];
        k3 = f(t + dt_half, y_temp);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt * k3[j];
        k4 = f(t + dt, y_temp);
        for (size_t j = 0; j < N; ++j)
            y[j] += dt_sixth * (k1[j] + k4[j]) + dt_third * (k2[j] + k3[j]);
        
        solution[i + 1]    = y;
        time_points[i + 1] = t + dt;
    }

    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

// Basic interface wrapper for rk4_solver
inline Matrix rk4_solver(
    MyFunc       deriv,
    const dVec& y0,
    double       t0,
    double       t1,
    double       dt
)
{
    auto params               = SolverParameters{};
    params.derivative         = deriv;
    params.initial_conditions = y0;
    params.t0                 = t0;
    params.t1                 = t1;
    params.dt                 = dt;
    return rk4_solver(params).solution;
} 