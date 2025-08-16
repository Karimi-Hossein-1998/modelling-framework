#pragma once
#include "../../../typedefs/header.hpp"
#include "../../../utility/adaptive_solvers_utility.hpp"
#include "../../../utility/dde_utility.hpp"

// Heun-Euler RK1(2) General DDE solver for multiple delays
SolverResults heun_euler_dde_general_solver(const SolverParameters& Params)
{
    auto deriv = Params.delay_derivative_multi;
    auto y     = Params.initial_conditions;
    auto t     = Params.t0;
    auto t1    = Params.t1;
    auto dt    = Params.dt;

    Matrix result;
    dVec   Times;
    dVec   Errors;
    dVec   Steps;
    result.push_back(y);
    Times.push_back(t);

    Matrix      Errors_History;
    Matrix      Steps_History;
    wVec Failed_Trials;
    dVec        Error_Attempt(Params.max_trial, 0.0);
    dVec        Step_Attempt(Params.max_trial, 0.0);

    auto y_temp = y;
    auto y_heun = y;
    auto y_old  = y;
    auto k1     = y;
    auto k2     = y;
    Matrix y_delayed_states(Params.delay_times.size(), dVec(y.size(), 0.0));

    // Case 1: No error estimation - simple Euler
    if (!Params.error_estimate)
    {
        while (t < t1)
        {
            y_delayed_states = interpolate_multiple_delays(Params.delay_times, t, Times, result, Params.t0, Params.history_function);
            k1 = deriv(t, y, y_delayed_states);
            for (size_t i = 0; i < y.size(); ++i)
                y[i] += dt * k1[i];
            
            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 2: Error estimation without variable steps
    else if (!Params.variable_steps)
    {
        y_old        = y;
        y_heun       = y_old;
        auto dt_half = dt * 0.5;
        while (t < t1)
        {
            // Compute Heun solution (2nd order)
            y_delayed_states = interpolate_multiple_delays(Params.delay_times, t, Times, result, Params.t0, Params.history_function);
            k1 = deriv(t, y, y_delayed_states);
            
            // Second stage
            y_heun = y;
            for (size_t i = 0; i < y_heun.size(); ++i)
                y_heun[i] += dt * k1[i];
            y_delayed_states = interpolate_multiple_delays(Params.delay_times, t + dt, Times, result, Params.t0, Params.history_function);
            k2 = deriv(t + dt, y_heun, y_delayed_states);

            // Update using Heun
            for (size_t i = 0; i < y_heun.size(); ++i)
                y_heun[i] = y[i] + dt_half * (k1[i] + k2[i]);

            // Compute Euler solution (1st order)
            for (size_t i = 0; i < y.size(); ++i)
                y[i] = y[i] + dt * k1[i];

            // Calculate error
            auto error = estimate_error(y, y_heun, Params);
            Errors.push_back(error);
            Steps.push_back(dt);

            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 3: Full adaptive step size with error estimation
    else
    {
        while (t < t1)
        {
            bool   step_accepted = false;
            y_old                = y;
            double current_dt    = dt;
            size_t trials        = 0;

            while (!step_accepted)
            {
                y            = y_old;  // Reset to previous state
                y_heun       = y_old;
                auto dt_half = current_dt * 0.5;
                
                // Compute Heun solution (2nd order)
                y_delayed_states = interpolate_multiple_delays(Params.delay_times, t, Times, result, Params.t0, Params.history_function);
                k1 = deriv(t, y_old, y_delayed_states);
                
                // Second stage
                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += current_dt * k1[i];
                y_delayed_states = interpolate_multiple_delays(Params.delay_times, t + current_dt, Times, result, Params.t0, Params.history_function);
                k2 = deriv(t + current_dt, y_temp, y_delayed_states);

                // Update using Heun
                for (size_t i = 0; i < y_heun.size(); ++i)
                    y_heun[i] = y_old[i] + dt_half * (k1[i] + k2[i]);

                // Compute Euler solution (1st order)
                for (size_t i = 0; i < y.size(); ++i)
                    y[i] = y_old[i] + current_dt * k1[i];

                // Calculate error and adjust step size
                auto error = estimate_error(y, y_heun, Params);
                Error_Attempt[trials] = error;
                Step_Attempt[trials] = current_dt;
                trials++;

                auto new_dt = adjust_step_size(current_dt, error, Params, 1);
                if (new_dt != current_dt)
                {
                    current_dt = new_dt;
                    continue;
                }

                step_accepted = true;
                dt = new_dt;
            }

            if (Params.attempts_history)
            {
                Errors_History.push_back(Error_Attempt);
                Steps_History.push_back(Step_Attempt);
                Failed_Trials.push_back(trials);
            }

            Errors.push_back(Error_Attempt[trials - 1]);
            Steps.push_back(current_dt);
            t += current_dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }

    SolverResults Results;
    Results.solution = result;
    Results.time_points = Times;
    if (Params.error_estimate)
    {
        Results.errors = Errors;
        if (Params.variable_steps)
        {
            Results.step_sizes = Steps;
            if (Params.attempts_history)
            {
                Results.errors_history = Errors_History;
                Results.steps_history = Steps_History;
                Results.failed_trials = Failed_Trials;
            }
        }
    }
    return Results;
}

// Bogacki-Shampine RK2(3) General DDE solver for multiple delays
SolverResults bogacki_shampine_dde_general_solver(const SolverParameters& Params)
{
    auto deriv = Params.delay_derivative_multi;
    auto y     = Params.initial_conditions;
    auto t     = Params.t0;
    auto t1    = Params.t1;
    auto dt    = Params.dt;

    Matrix result;
    dVec   Times;
    dVec   Errors;
    dVec   Steps;
    result.push_back(y);
    Times.push_back(t);

    Matrix      Errors_History;
    Matrix      Steps_History;
    wVec Failed_Trials;
    dVec        Error_Attempt(Params.max_trial, 0.0);
    dVec        Step_Attempt(Params.max_trial, 0.0);

    auto y_temp = y;
    auto y_rk3  = y;
    auto y_old  = y;
    auto k1     = y;
    auto k2     = y;
    auto k3     = y;
    Matrix y_delayed_states(Params.delay_times.size(), dVec(y.size(), 0.0));

    // Case 1: No error estimation - simple RK2
    if (!Params.error_estimate)
    {
        auto dt_half = dt * 0.5;
        while (t < t1)
        {
            y_delayed_states = interpolate_multiple_delays(Params.delay_times, t, Times, result, Params.t0, Params.history_function);
            k1 = deriv(t, y, y_delayed_states);
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_half * k1[i];
            y_delayed_states = interpolate_multiple_delays(Params.delay_times, t + dt_half, Times, result, Params.t0, Params.history_function);
            k2 = deriv(t + dt_half, y_temp, y_delayed_states);
            
            for (size_t i = 0; i < y.size(); ++i)
                y[i] += dt * k2[i];
            
            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 2: Error estimation without variable steps
    else if (!Params.variable_steps)
    {
        y_rk3         = y;
        auto dt_half  = dt * 0.5; 
        auto dt_34    = dt * 0.75;
        auto dt_14    = dt * 0.25;
        auto dt_third = dt / 3.0;
        auto dt_29    = dt * 2.0 / 9.0;
        auto dt_49    = dt * 4.0 / 9.0;
        while (t < t1)
        {
            // Compute RK3 solution
            y_delayed_states = interpolate_multiple_delays(Params.delay_times, t, Times, result, Params.t0, Params.history_function);
            k1 = deriv(t, y, y_delayed_states);
            
            // Second stage
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_half * k1[i];
            y_delayed_states = interpolate_multiple_delays(Params.delay_times, t + dt_half, Times, result, Params.t0, Params.history_function);
            k2 = deriv(t + dt_half, y_temp, y_delayed_states);

            // Third stage
            y_temp = y;
            for (size_t i = 0; i < y_temp.size(); ++i)
                y_temp[i] += dt_34 * k2[i];
            y_delayed_states = interpolate_multiple_delays(Params.delay_times, t + dt_34, Times, result, Params.t0, Params.history_function);
            k3 = deriv(t + dt_34, y_temp, y_delayed_states);

            // Update using RK3
            for (size_t i = 0; i < y_rk3.size(); ++i)
                y_rk3[i] = y[i] + dt_29 * k1[i] + dt_third * k2[i] + dt_49 * k3[i];

            // Compute RK2 solution
            for (size_t i = 0; i < y.size(); ++i)
                y[i] = y[i] + dt_14 * k1[i] + dt_34 * k3[i];

            // Calculate error
            auto error = estimate_error(y, y_rk3, Params);
            Errors.push_back(error);
            Steps.push_back(dt);

            t += dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }
    // Case 3: Full adaptive step size with error estimation
    else
    {
        while (t < t1)
        {
            bool   step_accepted = false;
            y_old                = y;
            double current_dt    = dt;
            size_t trials        = 0;

            while (!step_accepted)
            {
                y = y_old;  // Reset to previous state
                auto dt_half  = current_dt * 0.5; 
                auto dt_34    = current_dt * 0.75;
                auto dt_14    = current_dt * 0.25;
                auto dt_third = current_dt / 3.0;
                auto dt_29    = current_dt * 2.0 / 9.0;
                auto dt_49    = current_dt * 4.0 / 9.0;

                // Compute RK3 solution
                y_rk3 = y_old;
                y_delayed_states = interpolate_multiple_delays(Params.delay_times, t, Times, result, Params.t0, Params.history_function);
                k1 = deriv(t, y_old, y_delayed_states);
                
                // Second stage
                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += dt_half * k1[i];
                y_delayed_states = interpolate_multiple_delays(Params.delay_times, t + dt_half, Times, result, Params.t0, Params.history_function);
                k2 = deriv(t + dt_half, y_temp, y_delayed_states);

                // Third stage
                y_temp = y_old;
                for (size_t i = 0; i < y_temp.size(); ++i)
                    y_temp[i] += dt_34 * k2[i];
                y_delayed_states = interpolate_multiple_delays(Params.delay_times, t + dt_34, Times, result, Params.t0, Params.history_function);
                k3 = deriv(t + dt_34, y_temp, y_delayed_states);

                // Update using RK3
                for (size_t i = 0; i < y_rk3.size(); ++i)
                    y_rk3[i] = y_old[i] + dt_29 * k1[i] + dt_third * k2[i] + dt_49 * k3[i];

                // Compute RK2 solution
                for (size_t i = 0; i < y.size(); ++i)
                    y[i] = y_old[i] + dt_14 * k1[i] + dt_34 * k3[i];

                // Calculate error and adjust step size
                auto error = estimate_error(y, y_rk3, Params);
                Error_Attempt[trials] = error;
                Step_Attempt[trials] = current_dt;
                trials++;

                auto new_dt = adjust_step_size(current_dt, error, Params, 2);
                if (new_dt != current_dt)
                {
                    current_dt = new_dt;
                    continue;
                }

                step_accepted = true;
                dt = new_dt;
            }

            if (Params.attempts_history)
            {
                Errors_History.push_back(Error_Attempt);
                Steps_History.push_back(Step_Attempt);
                Failed_Trials.push_back(trials);
            }

            Errors.push_back(Error_Attempt[trials - 1]);
            Steps.push_back(current_dt);
            t += current_dt;
            result.push_back(y);
            Times.push_back(t);
        }
    }

    SolverResults Results;
    Results.solution = result;
    Results.time_points = Times;
    if (Params.error_estimate)
    {
        Results.errors = Errors;
        if (Params.variable_steps)
        {
            Results.step_sizes = Steps;
            if (Params.attempts_history)
            {
                Results.errors_history = Errors_History;
                Results.steps_history = Steps_History;
                Results.failed_trials = Failed_Trials;
            }
        }
    }
    return Results;
} 