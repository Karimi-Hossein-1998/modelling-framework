#pragma once
#include "../../../typedefs/header.hpp"
#include "../../../utility/adaptive_solvers_utility.hpp"
#include "../../../utility/dde_utility.hpp"

// RK12 General DDE solver with restructured error estimation
SolverResults rk12_dde_general_custom_solver(const SolverParameters& Params)
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

    // Case 1: No error estimation - simple RK1
    if (!Params.error_estimate)
    {
        while (t < t1)
        {
            auto y_delayed_states = interpolate_multiple_delays(Params.delay_times, t, Times, result, Params.t0, Params.history_function);
            auto k1 = deriv(t, y, y_delayed_states);
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
        auto dt_half = dt * 0.5;
        auto y_old   = y;
        auto y_rk2   = y_old;
        while (t < t1)
        {
            auto y_delayed_states = interpolate_multiple_delays(Params.delay_times, t, Times, result, Params.t0, Params.history_function);
            auto k1 = deriv(t, y_old, y_delayed_states);
            y_rk2 = y;
            y_old = y;
            for (size_t i = 0; i < y.size(); ++i)
                y_old[i] += dt_half * k1[i];
            y_delayed_states = interpolate_multiple_delays(Params.delay_times, t + dt_half, Times, result, Params.t0, Params.history_function);
            auto k2 = deriv(t + dt_half, y_old, y_delayed_states);
            for (size_t i = 0; i < y_rk2.size(); ++i)
                y_rk2[i] += dt * k2[i];
            for (size_t i = 0; i < y.size(); ++i)
                y[i] += dt * k1[i];
            auto error = estimate_error(y, y_rk2, Params);
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
            auto   y_old         = y;
            double current_dt    = dt;
            size_t trials        = 0;
            while (!step_accepted)
            {
                y = y_old;
                auto y_rk2   = y_old;
                auto dt_half = current_dt * 0.5;
                auto y_delayed_states = interpolate_multiple_delays(Params.delay_times, t, Times, result, Params.t0, Params.history_function);
                auto k1 = deriv(t, y_old, y_delayed_states);
                auto y_mid = y_old;
                for (size_t i = 0; i < y_mid.size(); ++i)
                    y_mid[i] += dt_half * k1[i];
                y_delayed_states = interpolate_multiple_delays(Params.delay_times, t + dt_half, Times, result, Params.t0, Params.history_function);
                auto k2 = deriv(t + dt_half, y_mid, y_delayed_states);
                for (size_t i = 0; i < y_rk2.size(); ++i)
                    y_rk2[i] += current_dt * k2[i];
                for (size_t i = 0; i < y.size(); ++i)
                    y[i] += current_dt * k1[i];
                auto error = estimate_error(y, y_rk2, Params);
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

// RK23 General DDE solver with restructured error estimation
// ... (implement following the same pattern as above)
// RK34 General DDE solver with restructured error estimation
// ... (implement following the same pattern as above)
// RK45 General DDE solver with restructured error estimation
// ... (implement following the same pattern as above) 