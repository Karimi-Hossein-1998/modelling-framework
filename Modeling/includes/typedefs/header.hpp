#pragma once    
#include <vector>
#include <thread>
#include <array>
#include <functional>
#include <cmath>
#include <iostream>
#include <random>
#include <utility>
#include <algorithm>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <ostream>
// -----------------------------------------------------------------------------
// Constant expression for pi
constexpr double PI        = 3.14159265358979323846;
constexpr int    max_order = 10;
constexpr double myexp     = 2.7;

// -----------------------------------------------------------------------------
// Type aliases for clarity and convenience
template<typename T>                                                              // A template for std::vector
using Vec            = std::vector<T>;                                            // that is basically a renaming (Vec = std::vector)
using dVec           = Vec<double>;                                               // 1D vector of doubles
using wVec           = Vec<size_t>;
using triple         = std::array<double, 3>;
using Matrix         = Vec<dVec >;                                                // 2D matrix (vector of vectors)
using func           = std::function<double(double)>;
using Func           = std::function<dVec(double)>;
using MyFunc         = std::function<dVec(double, const dVec&)>;                  // ODE function type
using DelayFunc      = std::function<dVec(double, const dVec&, const dVec&)>;     // DDE function type (single delay)
using DelayFuncMulti = std::function<dVec(double, const dVec&, const Matrix&)>;   // DDE function type (multiple delays)

using interPolater   = std::function<func(const dVec&, const dVec&)>;
using InterPolater   = std::function<Func(const dVec&, const Matrix&)>;
// Type definition for solvers
// -------------------------------------------------------------------------------------------------------------
struct SolverParameters 
{
    MyFunc derivative;                           // ODE function
    dVec    initial_conditions;
    double t0                    = 0.0;         // Start time
    double t1                    = 10.0;        // End time
    double dt                    = 1.0 / 16.0;  // Initial step size
    int    order                 = 4;           // Order for AB/ABM
    int    iterations            = 2;           // Iterations for ABM
    double min_dt                = 1.0 / 128.0; // Minimum step size (for variable step size control)
    double max_dt                = 0.5;         // Maximum step size (for variable step size control)
    double decrease_factor       = 0.9;
    double increase_factor       = 2.0;
    double local_tol             = 1.0e-7;      // Local error tolerance
    double absolute_tol          = 1.0e-10;     // Takes care of division by zero (NaN and/or inf values)
    double local_tol_error_ratio = 0.1;  
    size_t max_trial             = 25;
    bool   error_estimate        = false;      // Wether to use error estimation (or not)
    bool   variable_steps        = false;      // Wether to use variable step size (or not)
    bool   attempts_history      = false;      // Controls wether the errors and step sizes are stored for all attempts or not.
    bool   weighted_error        = true;       // Calculate error based on weighted error formula
    bool   norm_error            = false;      // Calculate error based on norm error formula (norm or the solution is used instead of maximum of the absolute value of the elements of the solution)
    // DDE-specific parameters
    DelayFunc delay_derivative   = nullptr;     // DDE function f(t, y, y_delayed)
    DelayFuncMulti delay_derivative_multi = nullptr; // DDE function f(t, y, y_delayed_states) for multiple delays
    dVec      delay_times;                      // Vector of delay times τ_i
    Func      history_function;                 // Initial history function y(t) for t < t0
    bool      is_dde            = false;        // Flag to indicate DDE vs ODE
    size_t    max_delay_order   = 10;           // Maximum order for delay terms
    double    interpolation_tol = 1.0e-6;       // Tolerance for delay interpolation
};

struct SolverResults 
{
    Matrix      solution;       // The solution matrix.
    dVec         time_points;    // Time points where solution is computed.
    dVec         errors;         // Error estimates.
    dVec         step_sizes;     // Actual step sizes used.
    Matrix      steps_history;  // History of step sizes used (contains failed steps attempted).
    Matrix      errors_history; // History of errors (contains failed steps attempted).
    wVec failed_trials;  // Number of failed trials in each step.
};

// -----------------------------------------------------------------------------
// Sparse vector and matrix structures for efficient storage of sparse data
struct SparsedVec  
{
    std::vector<std::pair<size_t, double>> values; // (index, value) pairs
    SparsedVec(size_t N) : values(N) {}
    void add(size_t i, double value) { values.emplace_back(i, value); }
};

struct SparseMatrix {
    std::vector<std::vector<std::pair<size_t, double>>> rows; // Each row: (col, value) pairs
    SparseMatrix(size_t N) : rows(N) {}
    void add(size_t i, size_t j, double weight) { rows[i].emplace_back(j, weight); }
};
