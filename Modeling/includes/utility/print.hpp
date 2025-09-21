#pragma once
#include "../typedefs/header.hpp"

// -----------------------------------------------------------------------------
// Utility functions for printing matrices/vectors
//
// Note: These are useful for debugging and data export. For large matrices/vectors,
// consider writing to file instead of printing to console to avoid overwhelming output.

// Print a matrix to std::cout (space-separated, row by row)
template <typename t>
void print_matrix(const Vec<t>& matrix)
{
    for (const auto& row : matrix)
    {
        for (const auto& val : row)
        {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }
}       

// Print a vector to std::cout (space-separated)
template <typename t>
void print_vector(const Vec<t>& vector)
{
    for (const auto& val : vector)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

// Call this function to print the progress bar
template <typename T>
void print_progress_bar_eta(size_t current, size_t total, T eta, size_t bar_width = 50)
{
    double progress = static_cast<double>(current) / total;
    size_t pos = static_cast<size_t>(bar_width * progress);
    std::cout << "\rProgress: [";
    for (size_t i = 0; i < bar_width; ++i)
    {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << ".";
    }
    std::cout << "] ";
    auto percentage = static_cast<size_t>(progress * 100.0);
    if ( percentage != 100)
        std::cout << std::setw(3) << percentage << "%";
    else
        std::cout << "Done!";
    std::cout << " : eta = " << eta;
    std::cout.flush();
}

void print_progress_bar(size_t current, size_t total, size_t bar_width = 50)
{
    double progress = static_cast<double>(current) / total;
    size_t pos = static_cast<size_t>(bar_width * progress);
    std::cout << "\rProgress: [";
    for (size_t i = 0; i < bar_width; ++i)
    {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << ".";
    }
    std::cout << "] ";
    auto percentage = static_cast<size_t>(progress * 100.0);
    if ( percentage != 100)
        std::cout << std::setw(3) << percentage << "%";
    else
        std::cout << "Done!";
    std::cout.flush();
}