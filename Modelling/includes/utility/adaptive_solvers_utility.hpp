#pragma once
#include "../typedefs/header.hpp"

// Helper function for error estimation
double estimate_error(const dVec& y, const dVec& y_ref, const SolverParameters& Params)
{
    double error = 0.0;
    if (Params.weighted_error)
    {
        // Weighted component-wise error
        for (size_t i = 0; i < y.size(); ++i)
        {
            double diff = std::abs(y[i] - y_ref[i]);
            error = std::max(error, diff / ( ( Params.absolute_tol + 1.0 ) * ( std::abs(y_ref[i]) + std::abs(y[i]) ) ));
        }
    }
    else if (Params.norm_error)
    {
        // L2 norm-based error
        double norm_diff = 0.0;
        double norm_sol  = 0.0;
        for (size_t i = 0; i < y.size(); ++i)
        {
            double diff = y[i] - y_ref[i];
            norm_diff += diff * diff;
            norm_sol += y_ref[i] * y_ref[i];
        }
        error = std::sqrt(norm_diff) / ( ( Params.absolute_tol + 1.0 ) * std::sqrt(norm_sol));
    }
    else
    {
        for (size_t i = 0; i < y.size(); ++i)
        {
            double diff = std::abs(y[i] - y_ref[i]);
            error = std::max(error, diff / ( ( Params.absolute_tol + 1.0 ) * std::abs(y_ref[i])));
        }
    }
    return error;
}

// Helper function for step size adjustment
double adjust_step_size(double current_dt, double error, const SolverParameters& Params, size_t p)
{
    if (error > Params.local_tol && current_dt > Params.min_dt)
    {
        return std::max(std::min(current_dt * Params.decrease_factor * std::pow(Params.local_tol / error, 1.0 / (p + 1.0)), 
                                Params.max_dt), 
                       Params.min_dt);
    }
    else if (error < Params.local_tol * Params.local_tol_error_ratio || current_dt <= Params.min_dt)
    {
        return std::max(std::min(current_dt * Params.increase_factor * std::pow(Params.local_tol / error, 1.0 / (p + 1.0)), 
                                Params.max_dt), 
                       Params.min_dt);
    }
    return current_dt;
}