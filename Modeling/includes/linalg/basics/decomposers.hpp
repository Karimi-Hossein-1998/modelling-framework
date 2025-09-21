#pragma once
#include "matrix_operations.hpp"
#include "eigen.hpp"

// Contains functions for performing various matrix decompositions.

// Returns LU decomposition with partial pivoting (PA = LU)
// Returns a tuple of {L, U, P} where P is a permutation vector.
inline std::tuple<Matrix, Matrix, wVec> lu_decomposition(Matrix A)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n)
        throw std::invalid_argument("Matrix must be square for LU decomposition.");

    wVec P(n);
    for (size_t i = 0; i < n; ++i) P[i] = i;

    for (size_t k = 0; k < n; ++k)
    {
        size_t max_row = k;
        for (size_t i = k + 1; i < n; ++i)
            if (std::abs(A[i][k]) > std::abs(A[max_row][k]))
                max_row = i;
        if (std::abs(A[max_row][k]) < 1e-12)
            throw std::runtime_error("Matrix is singular or nearly singular.");

        std::swap(A[k], A[max_row]);
        std::swap(P[k], P[max_row]);

        for (size_t i = k + 1; i < n; ++i)
        {
            A[i][k] /= A[k][k];
            for (size_t j = k + 1; j < n; ++j)
                A[i][j] -= A[i][k] * A[k][j];
        }
    }

    Matrix L(n, dVec(n, 0.0));
    Matrix U(n, dVec(n, 0.0));
    for (size_t i = 0; i < n; ++i)
    {
        L[i][i] = 1.0;
        for (size_t j = 0; j < n; ++j)
        {
            if (i > j)
                L[i][j] = A[i][j];
            else
                U[i][j] = A[i][j];
        }
    }
    return {L, U, P};
}

// Returns QR decomposition using modified Gram-Schmidt for a m x n matrix A.
// Returns a pair of {Q, R}.
inline std::pair<Matrix, Matrix> qr_decomposition(const Matrix& A)
{
    const size_t m = A.size();
    if (m == 0) return {{}, {}};
    const size_t n = A[0].size();

    Matrix Q(m, dVec(n, 0.0));
    Matrix R(n, dVec(n, 0.0));
    Matrix V = A; // Make a mutable copy of A

    for (size_t j = 0; j < n; ++j)
    {
        double norm_v = 0.0;
        for (size_t i = 0; i < m; ++i)
            norm_v += V[i][j] * V[i][j];
        R[j][j] = std::sqrt(norm_v);

        if (R[j][j] < 1e-12)
            throw std::runtime_error("Matrix is rank deficient during QR decomposition.");

        for (size_t i = 0; i < m; ++i)
            Q[i][j] = V[i][j] / R[j][j];

        for (size_t k = j + 1; k < n; ++k)
        {
            R[j][k] = 0.0;
            for (size_t i = 0; i < m; ++i)
                R[j][k] += Q[i][j] * V[i][k];

            for (size_t i = 0; i < m; ++i)
                V[i][k] -= R[j][k] * Q[i][j];
        }
    }
    return {Q, R};
}

// Returns Cholesky decomposition L for a symmetric positive definite matrix A (A = LL^T).
inline Matrix cholesky_decomposition(const Matrix& A)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n)
        throw std::invalid_argument("Matrix must be square for Cholesky decomposition.");

    Matrix L(n, dVec(n, 0.0));
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            double sum = A[i][j];
            for (size_t k = 0; k < j; ++k)
                sum -= L[i][k] * L[j][k];
            if (i == j)
            {
                if (sum <= 0.0)
                    throw std::runtime_error("Matrix is not positive definite.");
                L[i][j] = std::sqrt(sum);
            }
            else
            {
                if (std::abs(L[j][j]) < 1e-12)
                    throw std::runtime_error("Zero pivot in Cholesky decomposition.");
                L[i][j] = sum / L[j][j];
            }
        }
    }
    return L;
}

// Returns LDLT decomposition for a symmetric matrix A (A = LDL^T).
// Returns a pair of {L, D} where D is a diagonal vector.
inline std::pair<Matrix, dVec> ldlt_decomposition(const Matrix& A)
{
    const size_t n = A.size();
    if (A.empty() || A[0].size() != n)
        throw std::invalid_argument("Matrix must be square for LDLT decomposition.");

    Matrix L(n, dVec(n, 0.0));
    dVec D(n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            double sum = A[i][j];
            for (size_t k = 0; k < j; ++k)
                sum -= L[i][k] * D[k] * L[j][k];
            if (std::abs(D[j]) < 1e-12)
                throw std::runtime_error("Zero pivot in LDLT decomposition.");
            L[i][j] = sum / D[j];
        }
        double sum = A[i][i];
        for (size_t k = 0; k < i; ++k)
            sum -= L[i][k] * L[i][k] * D[k];
        D[i] = sum;
        L[i][i] = 1.0;
    }
    return {L, D};
}

// Returns Singular Value Decomposition (U, S, V) for a square matrix A (A = USV^T).
// Returns a tuple of {U, S, V}, where S is a vector of singular values.
inline std::tuple<Matrix, dVec, Matrix> svd(const Matrix& A)
{
    if (A.empty() || A.size() != A[0].size())
        throw std::invalid_argument("SVD currently supports square matrices only.");

    Matrix AtA = matrix_multiply(transpose(A), A);
    auto eigen_pair = eigen_qr(AtA);
    dVec s_values = eigen_pair.first;
    Matrix V = eigen_pair.second;

    for(double& val : s_values)
        val = std::sqrt(std::max(0.0, val));
    
    Matrix S_inv(s_values.size(), dVec(s_values.size(), 0.0));
    for(size_t i = 0; i < s_values.size(); ++i)
        if (s_values[i] > 1e-12)
            S_inv[i][i] = 1.0 / s_values[i];
    
    Matrix U = matrix_multiply(A, matrix_multiply(V, S_inv));
    
    return {U, s_values, V};
}

// Returns the Eigen decomposition (eigenvalues, eigenvectors) for a square matrix A.
inline std::pair<dVec, Matrix> eigen_decomposition(const Matrix& A)
{
    if (A.empty() || A.size() != A[0].size())
        throw std::invalid_argument("Matrix must be square for Eigen decomposition.");
    
    return eigen_qr(A);
} 