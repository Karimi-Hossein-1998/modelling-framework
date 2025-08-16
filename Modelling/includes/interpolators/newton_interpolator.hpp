#pragma once
#include "../typedefs/header.hpp"

// Newton interpolator: returns a function that interpolates at any x
inline func newton_interpolator(const dVec& x_points, const dVec& y_points)
{
    // Compute divided differences
    size_t n = x_points.size();
    dVec coef = y_points;
    for (size_t j = 1; j < n; ++j)
    {
        for (size_t i = n - 1; i >= j; --i)
        {
            coef[i] = (coef[i] - coef[i - 1]) / (x_points[i] - x_points[i - j]);
        }
    }
    return [x_points, coef](double x) -> double {
        double result = coef.back();
        for (int i = static_cast<int>(coef.size()) - 2; i >= 0; --i)
        {
            result = result * (x - x_points[i]) + coef[i];
        }
        return result;
    };
}

// Vector-valued Newton interpolator: returns a function that interpolates at any x
inline Func newton_interpolator(const dVec& x_points, const Matrix& y_points)
{
    size_t n = x_points.size();
    size_t N = y_points[0].size();
    // Compute divided differences for each component
    Matrix coef = y_points;
    for (size_t j = 1; j < n; ++j)
    {
        for (size_t i = n - 1; i >= j; --i)
        {
            for (size_t k = 0; k < N; ++k)
            {
                coef[i][k] = (coef[i][k] - coef[i - 1][k]) / (x_points[i] - x_points[i - j]);
            }
        }
    }
    return [x_points, coef, n, N](double x) -> dVec {
        dVec result(N, 0.0);
        for (size_t k = 0; k < N; ++k)
        {
            double value = coef[n - 1][k];
            for (int i = static_cast<int>(n) - 2; i >= 0; --i)
            {
                value = value * (x - x_points[i]) + coef[i][k];
            }
            result[k] = value;
        }
        return result;
    };
} 