#ifndef MATRIX_VECTOR_MUL_HPP
#define MATRIX_VECTOR_MUL_HPP

#include <vector>
#include <mpi.h>
#include <iostream>
/**
 * @brief Perform matrix-vector multiplication for a compressed sparse row (CSR) matrix.
 *
 * This function computes the product of a matrix, stored in compressed sparse row (CSR) format,
 * and a vector, and stores the result in the provided output vector.
 *
 * @param sub_A_values Non-zero values of the matrix.
 * @param sub_A_columns Column indices corresponding to the values in sub_A_values.
 * @param sub_A_rowptrs Row pointer array indicating the start of each row in sub_A_values.
 * @param v Input vector to be multiplied with the matrix.
 * @param sub_result Output vector to store the result of the multiplication.
 */
void mat_times_vec(const std::vector<double> &sub_A_values, const std::vector<int> &sub_A_columns, const std::vector<int> &sub_A_rowptrs, const std::vector<double> &v, std::vector<double> &sub_result);
void vec_lin_combination_result_1st_vec(double a, std::vector<double> &x, double b, const std::vector<double> &p);
double mpi_dot_product(const std::vector<double> &sub_u, const std::vector<double> &sub_v);
double dot_product(const std::vector<double> &u, const std::vector<double> &v);
double mpi_dot_product_debug(const std::vector<double> &sub_u, const std::vector<double> &sub_v);

void forward_substitution_mpi(  const std::vector<std::vector<int>> &which_row_which_rank,
                                // const std::vector<int> &L_rowptrs,
                                const std::vector<double> &sub_L_values,
                                const std::vector<int> &sub_L_columns,
                                const std::vector<int> &sub_L_rowptrs,
                                const std::vector<double> &sub_r,
                                std::vector<double> &sub_s,
                                MPI_Request *request);
void backward_substitution_mpi(  const std::vector<std::vector<int>> &which_row_which_rank,
                                // const std::vector<int> &L_rowptrs,
                                const std::vector<double> &sub_L_values,
                                const std::vector<int> &sub_L_columns,
                                const std::vector<int> &sub_L_rowptrs,
                                const std::vector<double> &sub_r,
                                std::vector<double> &sub_s,
                                MPI_Request *request);
#endif // MATRIX_VECTOR_MUL_HPP
