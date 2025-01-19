#ifndef READ_MATRIX_HPP
#define READ_MATRIX_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept> // For std::invalid_argument and std::out_of_range

// Function to read values, columns, and row pointers from text files
void read_val_col_rowptrs_from_txts(std::vector<double> &values, std::vector<int> &columns, std::vector<int> &row_ptrs, std::string data_name);
void read_b_from_txt(std::vector<double> &b, std::string data_name);

#endif // READ_MATRIX_HPP
