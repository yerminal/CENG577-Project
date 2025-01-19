#include <mpi.h>
#include "matrix_vector_mul.hpp"
#include "read_matrix.hpp"
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <chrono>
#include <thread>

std::vector<int> distributeScalarInteger(int scalar, int n) {
    std::vector<int> result(n, scalar / n); // Divide equally
    int remainder = scalar % n;

    // Distribute the remainder
    for (int i = 0; i < remainder; ++i) {
        result[i] += 1;
    }
    
    return result;
}

int main(int argc, char **argv) { // START
    std::string current_exec_name = argv[0];
    std::vector<std::string> all_args;

    if (argc > 1) {
        all_args.assign(argv + 1, argv + argc);
    }

    MPI_Init (&argc, &argv);
    int nprocs, rank;
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    double time_start, time_end;
    

    int max_iteration = 100000;
    std::vector<double> values;
    std::vector<int> columns;
    std::vector<int> row_ptrs;

    std::vector<double> b;
    std::string data_name = all_args[0];
    // Read the matrix and right hand side vector
    read_val_col_rowptrs_from_txts(values, columns, row_ptrs, data_name);
    read_b_from_txt(b, data_name);


    
    if(nprocs > 1) // START LOAD BALANCING FOR PARALLEL COMPUTATION
    {

        time_start = MPI_Wtime();

        std::vector<double> x(b.size(), 0.0);
        std::vector<double> p(b);

        int datasize_per_process = row_ptrs.back()/nprocs;
        int control_datasize_per_process = datasize_per_process;
        int index_rank = 0, total_used_procs = 0;
        int fix_loop_start_rank = -1;

        // FINDING THE RANK WHERE THE EQUAL DISTRIBUTION STARTS TO GO WRONG
        for(int i=0; i<row_ptrs.size(); i++){
            if(row_ptrs[i] < control_datasize_per_process){
                if(nprocs - total_used_procs == 0 && i >= row_ptrs.size()-1){
                    fix_loop_start_rank = 0;
                    break;
                }
                continue;
            }
            else{
                index_rank++;
                total_used_procs = index_rank;                
                if( ( (row_ptrs.back() - row_ptrs[i]) / datasize_per_process < nprocs - total_used_procs ) || (row_ptrs.size()-1-i < nprocs - total_used_procs) || nprocs - total_used_procs == 0){
                    fix_loop_start_rank = index_rank-1;
                    break;
                }
                if(row_ptrs[i] != row_ptrs.back()){
                    control_datasize_per_process = row_ptrs[i] + datasize_per_process;
                }
            }
        }

        index_rank = 0;
        control_datasize_per_process = datasize_per_process;

        int fix_loop_start_row_ptr = 0;
        
        std::vector<int> sub_row_ptrs;

        // Declarations of sub vectors
        std::vector<double> sub_r;
        std::vector<double> sub_p;
        std::vector<double> sub_x;
        
        // A vector to track which row belongs to which rank
        std::vector<std::vector<int>> which_row_which_rank (nprocs, std::vector<int>(2));

        // Main Loop for distribution of matrix (Load Balancing)
        if(fix_loop_start_rank != 0){
            // Max size can be datasize_per_process+1
            for(int i=0; i < row_ptrs.size(); i++){
                if(row_ptrs[i] < control_datasize_per_process){
                    if(rank == index_rank){
                        sub_row_ptrs.push_back(row_ptrs[i]);
                        sub_r.push_back(b[i]);
                        sub_p.push_back(b[i]);
                        sub_x.push_back(0.0);
                    }
                }
                else{
                    if(rank == index_rank){
                        sub_row_ptrs.push_back(row_ptrs[i]);
                    }
                    which_row_which_rank[index_rank][1] = i-1;
                    index_rank++;
                    if(index_rank >= nprocs){
                        break;
                    }

                    if(index_rank == fix_loop_start_rank){
                        fix_loop_start_row_ptr = i;
                        break;
                    }

                    which_row_which_rank[index_rank][0] = i;

                    control_datasize_per_process = row_ptrs[i] + datasize_per_process;
                    if(rank == index_rank){
                        sub_row_ptrs.push_back(row_ptrs[i]);
                        sub_r.push_back(b[i]);
                        sub_p.push_back(b[i]);
                        sub_x.push_back(0.0);
                    }
                }
            }
        }

        // IF THERE IS A RANK WHERE THE EQUAL DISTRIBUTION STARTS TO GO WRONG
        // AFTER DISTRIBUTING NORMALLY, SPREAD THE REMAINING ROWS TO THE REMAINING RANKS
        if(fix_loop_start_rank != -1){
            int index = fix_loop_start_row_ptr;
            std::vector<int> remaining_rows_distribution( distributeScalarInteger(row_ptrs.size()-1-fix_loop_start_row_ptr, nprocs-fix_loop_start_rank) );
            for(int remain_row_per_rank : remaining_rows_distribution){
                which_row_which_rank[index_rank][0] = index;
                for(; index<fix_loop_start_row_ptr+remain_row_per_rank+1; index++){
                    // std::cout << "XXXXXXXXXXXXXXXXXXXXX index_rank: " << index_rank << " //////// index: " << index << "////// fix_loop_start_row_ptr+remain_row_per_rank+1: "<< fix_loop_start_row_ptr+remain_row_per_rank+1 << std::endl;
                    if(rank == index_rank){
                        sub_row_ptrs.push_back(row_ptrs[index]);
                        if(index < fix_loop_start_row_ptr+remain_row_per_rank){
                            sub_r.push_back(b[index]);
                            sub_p.push_back(b[index]);
                            sub_x.push_back(0.0);
                        }
                    }
                }
                index--;
                which_row_which_rank[index_rank][1] = index-1;
                index_rank++;
                fix_loop_start_row_ptr = index;
            }
        }
        
        auto start = values.begin() + sub_row_ptrs.front();
        auto end = values.begin() + sub_row_ptrs.back();

        std::vector<double> sub_values(start, end);

        auto start2 = columns.begin() + sub_row_ptrs.front();
        auto end2 = columns.begin() + sub_row_ptrs.back();
        
        std::vector<int> sub_columns(start2, end2);

        std::vector<int> row_cnts_per_rank (nprocs);
        std::vector<int> row_starts_per_rank (nprocs);

        for(int i = 0; i < nprocs; i++){
            row_cnts_per_rank[i] = which_row_which_rank[i][1] - which_row_which_rank[i][0] + 1;
            row_starts_per_rank[i] = which_row_which_rank[i][0];
        }

        double dot_r = mpi_dot_product(sub_r, sub_r);
        double dot_r_old = dot_r;
        
        std::vector<double> sub_Ap (sub_row_ptrs.size()-1, 0.0);
        double dot_p_by_Ap;
        double alpha, beta;
        int conv_iter = -1;
        double algo_error = -1;

        // PARALLEL CG ALGORITHM STARTS after distribution of matrices and vectors
        for(int iter = 0; iter < max_iteration; iter++)
        {
            mat_times_vec(sub_values, sub_columns, sub_row_ptrs, p, sub_Ap);
            
            dot_p_by_Ap = mpi_dot_product(sub_p, sub_Ap);
            
            alpha = dot_r/dot_p_by_Ap;
            
            vec_lin_combination_result_1st_vec(1.0, sub_x, alpha, sub_p);
            vec_lin_combination_result_1st_vec(1.0, sub_r, -alpha, sub_Ap);
            dot_r = mpi_dot_product(sub_r, sub_r);
            if (sqrt(dot_r) < 0.0000000001) {
                conv_iter = iter;
                algo_error = sqrt(dot_r);
                break;
            }
            beta = dot_r/dot_r_old;
            dot_r_old = dot_r;

            vec_lin_combination_result_1st_vec(beta, sub_p, 1.0, sub_r);

            MPI_Allgatherv(&sub_p.front(), row_cnts_per_rank[rank], MPI_DOUBLE, &p.front(), &row_cnts_per_rank.front(), &row_starts_per_rank.front(), MPI_DOUBLE, MPI_COMM_WORLD);
        }

        MPI_Gatherv(&sub_x.front(), row_cnts_per_rank[rank], MPI_DOUBLE, &x.front(), &row_cnts_per_rank.front(), &row_starts_per_rank.front(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // PARALLEL CG ALGORITHM ENDS

        time_end = MPI_Wtime();

        if(rank == 0){
            
            // std::cout << "************************** PARALLEL SUMMARY **************************" << std::endl;
            // std::cout << std::endl;
            // std::cout << "PARALLEL Wall-Clock time: " << time_end-time_start << "s." << std::endl;
            // std::cout << "PARALLEL Converged at iter = " << conv_iter << std::endl;
            // std::cout << "PARALLEL Algorithm Error = " << algo_error << std::endl;
            if(nprocs < 10){
                std::cout << "0" << nprocs << ": " << time_end-time_start << std::endl;
            }
            else {
                std::cout << nprocs << ": " << time_end-time_start << std::endl;
            }

            //------------------------- Verification Test ----------------------------------------------------------------------
            // we will compare the A*x result to the right hand side
            std::vector<double> A_times_x(x.size());
            mat_times_vec(values, columns, row_ptrs, x, A_times_x);

            std::vector<double> error(x.size());
            for (int i = 0; i < x.size(); i++) {
                error[i] = abs(A_times_x[i] - b[i]);
            }
            if (*std::max_element(error.begin(), error.end()) > 0.0000000001)
                std::cout << "PARALLEL Error in solution is larger than " << 0.0000000001 << std::endl;
            // std::cout << "INFO - PARALLEL Error in solution is " << *std::max_element(error.begin(), error.end()) << std::endl;
        }    
    }
    else {
        // SEQUENTIAL CG ALGORITHM STARTS
        time_start = MPI_Wtime();
        std::vector<double> x(b.size(), 0.0);
        std::vector<double> p(b);
        std::vector<double> r(b);
        std::vector<double> Ap (b.size(), 0.0);
        
        double dot_r = dot_product(r, r);
        double dot_r_old = dot_r;
        
        double alpha, beta;
        int conv_iter = -1;
        double algo_error = -1;

        for(int iter = 0; iter < max_iteration; iter++)
        {
            mat_times_vec(values, columns, row_ptrs, p, Ap);
            alpha = dot_r_old/dot_product(Ap, p);

            vec_lin_combination_result_1st_vec(1.0, x, alpha, p);
            vec_lin_combination_result_1st_vec(1.0, r, -alpha, Ap);
            dot_r = dot_product(r, r);

            if (sqrt(dot_r) < 0.0000000001) {
                conv_iter = iter;
                algo_error = sqrt(dot_r);
                break;
            }
            beta = dot_r/dot_r_old;
            vec_lin_combination_result_1st_vec(beta, p, 1.0, r);
            dot_r_old = dot_r;
        }

        time_end = MPI_Wtime();
        // std::cout << "************************** SEQUENTIAL SUMMARY **************************" << std::endl;
        // std::cout << std::endl;
        // std::cout << "SEQUENTIAL Wall-Clock time: " << time_end-time_start << "s." << std::endl;
        // std::cout << "SEQUENTIAL Converged at iter = " << conv_iter << std::endl;
        // std::cout << "SEQUENTIAL Algorithm Error = " << algo_error << std::endl;
        std::cout << "0" << nprocs << ": " << time_end-time_start << std::endl;

        //------------------------- Verification Test ----------------------------------------------------------------------
        // we will compare the A*x result to the right hand side
        std::vector<double> A_times_x(x.size());
        std::vector<double> error(x.size());
        mat_times_vec(values, columns, row_ptrs, x, A_times_x);
        
        for (int i = 0; i < x.size(); i++) {
            error[i] = abs(A_times_x[i] - b[i]);
        }
        if (*std::max_element(error.begin(), error.end()) > 0.0000000001)
            std::cout << "SEQ Error in solution is larger than " << 0.0000000001 << std::endl;
        // std::cout << "INFO - SEQ Error in solution is " << *std::max_element(error.begin(), error.end()) << std::endl;
    }
    
    MPI_Finalize();
    return 0;
}