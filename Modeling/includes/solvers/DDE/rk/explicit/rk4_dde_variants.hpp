#pragma once
#include "../../../../typedefs/header.hpp"
#include "../../../../utility/dde_utility.hpp"

// 3/8 Rule RK4 for DDE
// A variant of the classical RK4 method using a different set of coefficients
// for the final weighted average, adapted for delay differential equations.
// Solves delay differential equations of the form:
//   dy/dt = f(t, y(t), y(t - tau))
// where tau is the delay time.
// Uses 3/8 rule: y_{n+1} = y_n + dt/8 * (k1 + 3*k2 + 3*k3 + k4)
inline SolverResults rk4_38_dde_solver(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f    = Params.delay_derivative;
    const auto&  y0   = Params.initial_conditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();
    const auto& delay_times      = Params.delay_times;
    const auto& history_function = Params.history_function;
    if (!f)
        throw std::invalid_argument("[rk4_38_dde_solver] DDE derivative function is required");
    if (delay_times.empty())
        throw std::invalid_argument("[rk4_38_dde_solver] At least one delay time is required for DDE");
    double max_delay = *std::max_element(delay_times.begin(), delay_times.end());
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = Matrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1);
    solution[0]    = y0;
    time_points[0] = t0;
    auto y      = y0;
    auto y_temp       = dVec(N, 0.0);
    auto y_delayed    = dVec(N, 0.0);
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    double t_delayed  = 0.0;
    auto dt_third     = dt / 3.0;
    auto dt_2_3       = dt * 2.0 / 3.0;
    auto dt_8         = dt / 8.0;
    auto dt_3_8       = dt * 3.0 / 8.0;
    
    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];
        
        // For single delay, interpolate the delayed solution at current time
        t_delayed = t - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // RK4 3/8 rule step 1: k1 = f(t, y, y_delayed)
        k1 = f(t, y, y_delayed);
        
        // RK4 3/8 rule step 2: Compute first intermediate values y_temp = y + dt/3 * k1
        for (size_t j = 0; j < N; ++j)
        {
            y_temp[j] = y[j] + dt_third * k1[j];
        }
        
        // Interpolate delayed solution at first intermediate time
        t_delayed = t + dt_third - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // RK4 3/8 rule step 3: k2 = f(t + dt/3, y_temp, y_delayed_1)
        k2 = f(t + dt_third, y_temp, y_delayed);
        
        // RK4 3/8 rule step 4: Compute second intermediate values y_temp = y - dt/3 * k1 + dt * k2
        for (size_t j = 0; j < N; ++j)
        {
            y_temp[j] = y[j] - dt_third * k1[j] + dt * k2[j];
        }
        
        // Interpolate delayed solution at second intermediate time
        t_delayed = t + dt_2_3 - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // RK4 3/8 rule step 5: k3 = f(t + 2*dt/3, y_temp, y_delayed_2)
        k3 = f(t + dt_2_3, y_temp, y_delayed);
        
        // RK4 3/8 rule step 6: Compute final intermediate values y_temp = y + dt * (k1 - k2 + k3)
        for (size_t j = 0; j < N; ++j)
        {
            y_temp[j] = y[j] + dt * (k1[j] - k2[j] + k3[j]);
        }
        
        // Interpolate delayed solution at full step time
        t_delayed = t + dt - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // RK4 3/8 rule step 7: k4 = f(t + dt, y_temp, y_delayed_full)
        k4 = f(t + dt, y_temp, y_delayed);
        
        // RK4 3/8 rule step 8: Update solution using 3/8 rule: y_{n+1} = y_n + dt/8 * (k1 + 3*k2 + 3*k3 + k4)
        for (size_t j = 0; j < N; ++j)
        {
            y[j] += dt_8 * (k1[j] + k4[j])  + dt_3_8 * (k2[j] + k3[j]);
        }
        
        solution[i + 1]    = y;
        time_points[i + 1] = t + dt;
    }
    
    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

// Basic interface wrapper for 3/8 Rule RK4 DDE solver
inline Matrix rk4_38_dde_solver(
    DelayFunc    deriv,
    const dVec& y0,
    double       t0,
    double       t1,
    double       dt,
    const dVec& delay_times      = {1.0},
    const Func&  history_function = {}
)
{
    auto params               = SolverParameters{};
    params.delay_derivative   = deriv;
    params.initial_conditions = y0;
    params.t0                 = t0;
    params.t1                 = t1;
    params.dt                 = dt;
    params.delay_times        = delay_times;
    params.history_function   = history_function;
    params.is_dde            = true;
    return rk4_38_dde_solver(params).solution;
}

// Ralston's RK4 for DDE
// An RK4 variant optimized for a smaller error bound with a specific set of
// coefficients for the stage evaluations and the final update, adapted for delay differential equations.
// Solves delay differential equations of the form:
//   dy/dt = f(t, y(t), y(t - tau))
// where tau is the delay time.
inline SolverResults rk4_ralston_dde_solver(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f    = Params.delay_derivative;
    const auto&  y0   = Params.initial_conditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();
    const auto& delay_times      = Params.delay_times;
    const auto& history_function = Params.history_function;
    if (!f)
        throw std::invalid_argument("[rk4_ralston_dde_solver] DDE derivative function is required");
    if (delay_times.empty())
        throw std::invalid_argument("[rk4_ralston_dde_solver] At least one delay time is required for DDE");
    double max_delay = *std::max_element(delay_times.begin(), delay_times.end());
    const size_t num_steps   = static_cast<size_t>((t1 - t0) / dt);
    auto         solution    = Matrix(num_steps + 1, dVec(N));
    auto         time_points = dVec(num_steps + 1);
    solution[0]    = y0;
    time_points[0] = t0;
    auto y      = y0;
    auto y_temp       = dVec(N, 0.0);
    auto y_delayed    = dVec(N, 0.0);
    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);
    double t_delayed  = 0.0;
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
    
    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];
        
        // For single delay, interpolate the delayed solution at current time
        t_delayed = t - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // Ralston RK4 step 1: k1 = f(t, y, y_delayed)
        k1 = f(t, y, y_delayed);
        
        // Ralston RK4 step 2: Compute first intermediate values y_temp = y + a21*dt*k1
        for (size_t j = 0; j < N; ++j)
        {
            y_temp[j] = y[j] + dt_2_5 * k1[j];
        }
        
        // Interpolate delayed solution at first intermediate time
        t_delayed = t + dt_2_5 - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // Ralston RK4 step 3: k2 = f(t + c2*dt, y_temp, y_delayed_1)
        k2 = f(t + dt_2_5, y_temp, y_delayed);
        
        // Ralston RK4 step 4: Compute second intermediate values y_temp = y + a31*dt*k1 + a32*dt*k2
        for (size_t j = 0; j < N; ++j)
        {
            y_temp[j] = y[j] + k3_1 * k1[j] + k3_2 * k2[j];
        }
        
        // Interpolate delayed solution at second intermediate time
        t_delayed = t + k3_t - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // Ralston RK4 step 5: k3 = f(t + c3*dt, y_temp, y_delayed_2)
        k3 = f(t + k3_t, y_temp, y_delayed);
        
        // Ralston RK4 step 6: Compute final intermediate values y_temp = y + a41*dt*k1 + a42*dt*k2 + a43*dt*k3
        for (size_t j = 0; j < N; ++j)
        {
            y_temp[j] = y[j] + k4_1 * (k1[j] + k2[j]) + k4_3 * k3[j];
        }
        
        // Interpolate delayed solution at full step time
        t_delayed = t + dt - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // Ralston RK4 step 7: k4 = f(t + dt, y_temp, y_delayed_full)
        k4 = f(t + dt, y_temp, y_delayed);
        
        // Ralston RK4 step 8: Update solution using Ralston coefficients
        for (size_t j = 0; j < N; ++j)
        {
            y[j] += dy_1 * (k1[j] + k2[j] + k4[j]) + dy_3 * k3[j];
        }
        
        solution[i + 1]    = y;
        time_points[i + 1] = t + dt;
    }
    
    auto results        = SolverResults{};
    results.solution    = solution;
    results.time_points = time_points;
    return results;
}

// Basic interface wrapper for Ralston's RK4 DDE solver
inline Matrix rk4_ralston_dde_solver(
    DelayFunc    deriv,
    const dVec& y0,
    double       t0,
    double       t1,
    double       dt,
    const dVec& delay_times      = {1.0},
    const Func&  history_function = {}
)
{
    auto params               = SolverParameters{};
    params.delay_derivative   = deriv;
    params.initial_conditions = y0;
    params.t0                 = t0;
    params.t1                 = t1;
    params.dt                 = dt;
    params.delay_times        = delay_times;
    params.history_function   = history_function;
    params.is_dde            = true;
    return rk4_ralston_dde_solver(params).solution;
}

// Gill's RK4 for DDE
// A memory-efficient RK4 variant that reduces the number of temporary
// variables required by reusing storage during the stage calculations,
// adapted for delay differential equations.
// Solves delay differential equations of the form:
//   dy/dt = f(t, y(t), y(t - tau))
// where tau is the delay time.
inline SolverResults rk4_gill_dde_solver(const SolverParameters& Params)
{
    // Extract parameters for clarity
    const auto&  f    = Params.delay_derivative;
    const auto&  y0   = Params.initial_conditions;
    const double t0   = Params.t0;
    const double t1   = Params.t1;
    const double dt   = Params.dt;
    const size_t N    = y0.size();
    const auto& delay_times      = Params.delay_times;
    const auto& history_function = Params.history_function;
    if (!f)
        throw std::invalid_argument("[rk4_gill_dde_solver] DDE derivative function is required");
    if (delay_times.empty())
        throw std::invalid_argument("[rk4_gill_dde_solver] At least one delay time is required for DDE");
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
    auto y_delayed    = dVec(N, 0.0);

    dVec  k1(N, 0.0), k2(N, 0.0), k3(N, 0.0), k4(N, 0.0);

    double t_delayed = 0.0;
    // Main integration loop
    for (size_t i = 0; i < num_steps; ++i)
    {
        const double t = time_points[i];
        
        // For single delay, interpolate the delayed solution at current time
        t_delayed = t - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // Gill RK4 step 1: k1 = dt/2 * f(t, y, y_delayed)
        k1     = f(t, y, y_delayed);
        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_half * k1[j];
        
        // Interpolate delayed solution at first intermediate time
        t_delayed = t + dt_half - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // Gill RK4 step 2: k2 = (1 - 1/sqrt(2)) * (dt*f2 - k1)
        k2 = f(t + dt_half, y_temp, y_delayed);
        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_a_half * k1[j] + dt_b_half * k2[j];
        
        // Interpolate delayed solution at second intermediate time
        t_delayed = t + dt_half - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // Gill RK4 step 3: k3 = (1 + 1/sqrt(2)) * (dt*f3 - k2)
        k3 = f(t + dt_half, y_temp, y_delayed);
        y_temp = y;
        for(size_t j = 0; j < N; ++j)
            y_temp[j] += dt_b * k3[j] - dt_s2inverse * k2[j];
        
        // Interpolate delayed solution at full step time
        t_delayed = t + dt - delay_times[0];
        y_delayed = interpolate_delay(t_delayed, time_points, solution, t0, history_function);
        
        // Gill RK4 step 4: f4 = f(t + dt, y_temp, y_delayed_full)
        k4 = f(t + dt, y_temp, y_delayed);
        
        // Gill RK4 step 5: Update solution using Gill coefficients
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

// Basic interface wrapper for Gill's RK4 DDE solver
inline Matrix rk4_gill_dde_solver(
    DelayFunc    deriv,
    const dVec& y0,
    double       t0,
    double       t1,
    double       dt,
    const dVec& delay_times      = {1.0},
    const Func&  history_function = {}
)
{
    auto params               = SolverParameters{};
    params.delay_derivative   = deriv;
    params.initial_conditions = y0;
    params.t0                 = t0;
    params.t1                 = t1;
    params.dt                 = dt;
    params.delay_times        = delay_times;
    params.history_function   = history_function;
    params.is_dde             = true;
    return rk4_gill_dde_solver(params).solution;
} 