#pragma once
#include "../../../../typedefs/header.hpp"
#include "../../../../utility/dde_utility.hpp"

// RK4 DDE coefficient constants for better code aesthetics
// 3/8 Rule RK4 for General DDE
// A variant of the classical RK4 method using a different set of coefficients
// for the final weighted average, adapted for general delay differential equations.
// Solves delay differential equations of the form:
// dy/dt = f(t, y(t), y(t - tau1), y(t - tau2), ..., y(t - taun))
// where tau1, tau2, ..., taun are multiple delay times
// Uses 3/8 rule: y_{n+1} = y_n + dt/8 * (k1 + 3*k2 + 3*k3 + k4)
inline SolverResults rk4_38_dde_general_solver(const SolverParameters& Params)
{
    const auto&  f    = Params.delay_derivative_multi;
    const auto&  y0   = Params.initial_conditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();
    const auto& delay_times      = Params.delay_times;
    const auto& history_function = Params.history_function;
    if (!f)
        throw std::invalid_argument("[rk4_38_dde_general_solver] DDE derivative function is required");
    if (delay_times.empty())
        throw std::invalid_argument("[rk4_38_dde_general_solver] At least one delay time is required for DDE");
    double max_delay = *std::max_element(delay_times.begin(), delay_times.end());
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = Matrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1);
    solution[0]    = y0;
    time_points[0] = t0;
    auto y      = y0;
    auto y_temp       = dVec(N, 0.0);
    Matrix y_delayed_states(delay_times.size(), dVec(N, 0.0));
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    auto dt_third     = dt / 3.0;
    auto dt_2_3       = dt_third * 2.0;
    auto dt_8         = dt / 8.0;
    auto dt_3_8       = dt_8 * 3.0;
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];
        y_delayed_states = interpolate_multiple_delays(delay_times, t, time_points, solution, t0, history_function);
        k1 = f(t, y, y_delayed_states);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_third * k1[j];
        y_delayed_states = interpolate_multiple_delays(delay_times, t + dt_third, time_points, solution, t0, history_function);
        k2 = f(t + dt_third, y_temp, y_delayed_states);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] - dt_third * k1[j] + dt * k2[j];
        y_delayed_states = interpolate_multiple_delays(delay_times, t + dt_2_3, time_points, solution, t0, history_function);
        k3 = f(t + dt_2_3, y_temp, y_delayed_states);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt * (k1[j] - k2[j] + k3[j]);
        y_delayed_states = interpolate_multiple_delays(delay_times, t + dt, time_points, solution, t0, history_function);
        k4 = f(t + dt, y_temp, y_delayed_states);
        for (size_t j = 0; j < N; ++j)
            y[j] += dt_8 * (k1[j] + k4[j]) + dt_3_8 * (k2[j] + k3[j]);
        solution[i + 1]    = y;
        time_points[i + 1] = t + dt;
    }
    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

inline Matrix rk4_38_dde_general_solver(
    DelayFuncMulti deriv,
    const dVec&   y0,
    double         t0,
    double         t1,
    double         dt,
    const dVec&   delay_times      = {1.0},
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
    return rk4_38_dde_general_solver(params).solution;
}

// General DDE version of Ralston's RK4 for multiple delays
// An RK4 variant optimized for a smaller error bound with a specific set of
// coefficients for the stage evaluations and the final update, adapted for delay differential equations.
// Solves delay differential equations of the form:
// dy/dt = f(t, y(t), y(t - τ₁), y(t - τ₂), ..., y(t - τₙ))
// where τ_i are the delay times
inline SolverResults rk4_ralston_dde_general_solver(const SolverParameters& Params)
{
    const auto&  f    = Params.delay_derivative_multi;
    const auto&  y0   = Params.initial_conditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();
    const auto& delay_times      = Params.delay_times;
    const auto& history_function = Params.history_function;
    if (!f)
        throw std::invalid_argument("[rk4_ralston_dde_general_solver] DDE derivative function is required");
    if (delay_times.empty())
        throw std::invalid_argument("[rk4_ralston_dde_general_solver] At least one delay time is required for DDE");
    double       max_delay   = *std::max_element(delay_times.begin(), delay_times.end());
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = Matrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1);
    solution[0]    = y0;
    time_points[0] = t0;
    auto y         = y0;
    auto y_temp    = dVec(N, 0.0);
    Matrix y_delayed_states(delay_times.size(), dVec(N, 0.0));
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    auto dt_2_5  = dt * 2.0 / 5.0;
    auto sqrt_5  = std::sqrt(5.0);
    auto dt_1_4  = dt * 0.25;
    auto dt_1_40 = dt_1_4 * 0.1;
    auto k3_t    = dt_1_4 * ( sqrt_5 - 1.0 );
    auto k3_1    = dt_1_40 * ( 11.0 - sqrt_5 );
    auto k3_2    = dt_1_40 * ( 25.0 - 13.0 * sqrt_5 );
    auto k4_1    = dt_1_40 * ( 10.0 - 7.0 * sqrt_5 );
    auto k4_3    = dt_1_4 * ( 5.0 + sqrt_5 );
    auto dt_1_20 = dt_1_4 * 0.2;
    auto dy_1    = dt_1_20 * ( 5.0 - sqrt_5 );
    auto dy_3    = dt_1_20 * ( 5.0 + 3.0 * sqrt_5 );
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];
        y_delayed_states = interpolate_multiple_delays(delay_times, t, time_points, solution, t0, history_function);
        k1 = f(t, y, y_delayed_states);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + dt_2_5 * k1[j];
        y_delayed_states = interpolate_multiple_delays(delay_times, t + dt_2_5, time_points, solution, t0, history_function);
        k2 = f(t + dt_2_5, y_temp, y_delayed_states);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + k3_1 * k1[j] + k3_2 * k2[j];
        y_delayed_states = interpolate_multiple_delays(delay_times, t + k3_t, time_points, solution, t0, history_function);
        k3 = f(t + k3_t, y_temp, y_delayed_states);
        for (size_t j = 0; j < N; ++j)
            y_temp[j] = y[j] + k4_1 * (k1[j] + k2[j]) + k4_3 * k3[j];
        y_delayed_states = interpolate_multiple_delays(delay_times, t + dt, time_points, solution, t0, history_function);
        k4 = f(t + dt, y_temp, y_delayed_states);
        for (size_t j = 0; j < N; ++j)
            y[j] += dy_1 * (k1[j] + k2[j] + k4[j]) + dy_3 * k3[j];
        solution[i + 1]    = y;
        time_points[i + 1] = t + dt;
    }
    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

inline Matrix rk4_ralston_dde_general_solver(
    DelayFuncMulti deriv,
    const dVec&   y0,
    double         t0,
    double         t1,
    double         dt,
    const dVec&   delay_times      = {1.0},
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
    return rk4_ralston_dde_general_solver(params).solution;
}

// General DDE version of Gill's RK4 for multiple delays
// A memory-efficient RK4 variant that reduces the number of temporary
// variables required by reusing storage during the stage calculations,
// adapted for delay differential equations.
// Solves delay differential equations of the form:
// dy/dt = f(t, y(t), y(t - τ₁), y(t - τ₂), ..., y(t - τₙ))
// where τ_i are the delay times
inline SolverResults rk4_gill_dde_general_solver(const SolverParameters& Params)
{
    const auto&  f    = Params.delay_derivative_multi;
    const auto&  y0   = Params.initial_conditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();
    const auto& delay_times      = Params.delay_times;
    const auto& history_function = Params.history_function;
    if (!f)
        throw std::invalid_argument("[rk4_gill_dde_general_solver] DDE derivative function is required");
    if (delay_times.empty())
        throw std::invalid_argument("[rk4_gill_dde_general_solver] At least one delay time is required for DDE");
    double max_delay = *std::max_element(delay_times.begin(), delay_times.end());
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = Matrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1);
    solution[0]    = y0;
    time_points[0] = t0;
    auto y            = y0;
    auto dt_half      = dt * 0.5;
    auto dt_sixth     = dt / 6.0;
    auto sqrt2        = std::sqrt(2.0);
    auto sqrt2inverse = sqrt2 / 2.0;
    auto dt_s2inverse = dt * sqrt2inverse;
    auto a            = 1.0 - sqrt2inverse;
    auto dt_a_half    = dt_half * a;
    auto b            = 1.0 + sqrt2inverse;
    auto dt_b_half    = dt_half * b;
    auto dt_b         = dt * b;
    auto dy_2         = dt_sixth * a * 2.0;
    auto dy_3         = dt_sixth * b * 2.0;
    auto y_temp       = dVec(N, 0.0);
    Matrix y_delayed_states(delay_times.size(), dVec(N, 0.0));
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];
        y_delayed_states = interpolate_multiple_delays(delay_times, t, time_points, solution, t0, history_function);
        k1 = f(t, y, y_delayed_states);
        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_half * k1[j];
        y_delayed_states = interpolate_multiple_delays(delay_times, t + dt_half, time_points, solution, t0, history_function);
        k2 = f(t + dt_half, y_temp, y_delayed_states);
        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_a_half * k1[j] + dt_b_half * k2[j];
        y_delayed_states = interpolate_multiple_delays(delay_times, t + dt_half, time_points, solution, t0, history_function);
        k3 = f(t + dt_half, y_temp, y_delayed_states);
        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_b * k3[j] - dt_s2inverse * k2[j];
        y_delayed_states = interpolate_multiple_delays(delay_times, t + dt, time_points, solution, t0, history_function);
        k4 = f(t + dt, y_temp, y_delayed_states);
        for (size_t j = 0; j < N; ++j)
            y[j] += dt_sixth * (k1[j] + k4[j]) + dy_2 * k2[j] + dy_3 * k3[j];
        solution[i + 1]    = y;
        time_points[i + 1] = t + dt;
    }
    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

inline Matrix rk4_gill_dde_general_solver(
    DelayFuncMulti deriv,
    const dVec&   y0,
    double         t0,
    double         t1,
    double         dt,
    const dVec&   delay_times      = {1.0},
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
    return rk4_gill_dde_general_solver(params).solution;
} 