#include "matrix_vector_mul.hpp"
#include <iomanip>
// Matrix times vector
void mat_times_vec(const std::vector<double> &sub_A_values, const std::vector<int> &sub_A_columns, const std::vector<int> &sub_A_rowptrs, const std::vector<double> &v, std::vector<double> &sub_result)
{
    double dot_prod;
    int start = 0, end = 0;

    for (int i = 0; i < sub_A_rowptrs.size()-1; i++) { // Per row iteration
        dot_prod = 0;
        start = sub_A_rowptrs[i]-sub_A_rowptrs.front();
        end = sub_A_rowptrs[i+1]-sub_A_rowptrs.front();
        for(; start < end; start++){
            dot_prod += sub_A_values[start]*v[sub_A_columns[start]];            
        }
        sub_result[i] = dot_prod;
    }

}


// Linear combination of vectors
void vec_lin_combination_result_1st_vec(double a, std::vector<double> &x, double b, const std::vector<double> &p)
{
    for (int i = 0; i < x.size(); i++)
        x[i] = a * x[i] + b * p[i];
}

double mpi_dot_product(const std::vector<double> &sub_u, const std::vector<double> &sub_v)
{
    double product = 0.0;
    int length = sub_u.size();

    double sub_prod = 0.0;
    for (int i = 0; i < length; i++) {
        sub_prod += sub_u[i] * sub_v[i];
    }

    // reduction over sub_prod to get the total dot product
    MPI_Allreduce(&sub_prod, &product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return product;
}

double mpi_dot_product_debug(const std::vector<double> &sub_u, const std::vector<double> &sub_v)
{
    int nprocs, rank;
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    double product = 0.0;
    int length = sub_u.size();
    if(rank == 1){
        std::cout << "length: " << length << std::endl;
        std::cout << "sub_u.size(): " << sub_u.size() << std::endl;
        std::cout << "sub_v.size(): " << sub_v.size() << std::endl;
    }
    double sub_prod = 0.0;
    for (int i = 0; i < length; i++) {
        if(rank == 1){
            std::cout << "i: " << i << ", sub_u[i]" << std::setprecision(17) << sub_u[i] << ", sub_v[i]" << std::setprecision(17) << sub_v[i] << ", sub_prod: " << std::setprecision(17) << sub_prod << std::endl;
        }
        sub_prod += sub_u[i] * sub_v[i];
    }

    std::cout << "rank: " << rank << ", sub_prod: " << std::setprecision(17) << sub_prod << std::endl;

    // reduction over sub_prod to get the total dot product
    MPI_Allreduce(&sub_prod, &product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return product;
}

double dot_product(const std::vector<double> &u, const std::vector<double> &v)
{
    int length = u.size();

    double prod = 0.0;

    for (int i = 0; i < length; i++) {
        prod += u[i] * v[i];
    }

    return prod;
}

// Forward substitution: Solves Ls = r
void forward_substitution_mpi(  const std::vector<std::vector<int>> &which_row_which_rank,
                                // const std::vector<int> &L_rowptrs,
                                const std::vector<double> &sub_L_values,
                                const std::vector<int> &sub_L_columns,
                                const std::vector<int> &sub_L_rowptrs,
                                const std::vector<double> &sub_r,
                                std::vector<double> &sub_s,
                                MPI_Request *request) {
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    sub_s.resize(sub_r.size(), 0.0);

    int which_row = which_row_which_rank[rank][0];
    // for(int i = 0; i<L_rowptrs.size(); i++){
    //     if(L_rowptrs[i] == sub_L_rowptrs.front()){
    //         which_row = i;
    //         break;
    //     }
    // }
    // for (int i = 0; i < n; ++i) {
    //     double sum = 0.0;
    //     if (i % size == rank) { // Row assigned to this process
    //         for (int j = 0; j < i; ++j)
    //             sum += L[i][j] * x[j];

    //         x[i] = (b[i] - sum) / L[i][i];
    //     }

    //     // Broadcast the computed value of x[i] to all processes
    //     MPI_Bcast(&x[i], 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
    // }


    // MPI_Request request[(nprocs-1-rank)*(sub_L_rowptrs.size()-1)];
    int request_counter = 0;
    ////////////
    int start, end;
    int index_rank = rank + 1;
    double sub_L_value_current_diagonal = 0;
    double received_s_value = 0.0;
    int source_rank;



    // if(rank == 1){
    //     std::cout << "rank: " << rank << std::endl;
    //     std::cout << "sub_L_rowptrs:" << std::endl;
    //     for (int ptr : sub_L_rowptrs) {
    //         std::cout << ptr << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    // if(rank == 1){
    //     std::cout << "rank: " << rank << std::endl;
    //     std::cout << "sub_L_columns:" << std::endl;
    //     for (int ptr : sub_L_columns) {
    //         std::cout << ptr << std::endl;
    //     }
    //     std::cout << std::endl;
    // }    
    
    for (int i = 0; i < sub_L_rowptrs.size()-1; i++) { // Per row iteration
        start = sub_L_rowptrs[i]-sub_L_rowptrs.front();
        end = sub_L_rowptrs[i+1]-sub_L_rowptrs.front()-1;

        if(sub_L_columns[end] == which_row + i){
            sub_s[i] += sub_r[i]/sub_L_values[end];
            sub_L_value_current_diagonal = sub_L_values[end];
            // if(rank == 1){
            //     std::cout << "above RANK: " << rank << "//// sub_s[" << i << "]: " << sub_s[i] << std::endl;
            //     std::cout << "above RANK: " << rank << "//// start: " << start << std::endl; 
            //     std::cout << "above RANK: " << rank << "//// end: " << end << std::endl;  
            // }   
        }
        else{
            std::cout << "HATATAA" << std::endl;
        }

        end--;
        for(; start <= end; end--){ //per element in a row
            // if(rank == 1){
            //     std::cout << "RANK: " << rank << "//// sub_L_columns[end]: " << sub_L_columns[end] << std::endl;
            //     std::cout << "RANK: " << rank << "//// which_row: " << which_row << std::endl; 
            // }
            if(which_row <= sub_L_columns[end]){
                // if(rank == 1){
                //     std::cout << "RANK: " << rank << "//// sub_s[" << i << "]: " << sub_s[i] << std::endl; 
                // }

                sub_s[i] -= sub_L_values[end]/sub_L_value_current_diagonal*sub_s[sub_L_columns[end]-which_row];

                // if(rank == 1){
                //     std::cout << "RANK: " << rank << "//// sub_s[" << i << "]: " << sub_s[i] << std::endl; 
                // }            
            }
            else{
                source_rank = -1;
                for(int x = 0; x < rank; x++){
                    if(which_row_which_rank[x][0] <= sub_L_columns[end] && sub_L_columns[end] <= which_row_which_rank[x][1]){
                        source_rank = x;
                        break;
                    }
                }

                MPI_Recv(&received_s_value, 1, MPI_DOUBLE, source_rank, sub_L_columns[end], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(rank == 1){
                    std::cout << "RANK: " << rank << "//// received_s_value: " << received_s_value << "which_row: " << which_row << " i: " << i <<  std::endl; 
                }  
                // std::cout << "RANK: " << rank << "//// received_s_value: " << received_s_value << std::endl; 
                sub_s[i] -= sub_L_values[end]/sub_L_value_current_diagonal*received_s_value;
            }
        }

        // if(rank == 0){
        //     std::cout << "RANK: " << rank << "//// sub_s[" << i << "]: " << sub_s[i] << std::endl; 
        // }
        // empty_entire_s[which_row+i] = sub_s[i];

        while(index_rank < nprocs){
            MPI_Isend(&sub_s[i], 1, MPI_DOUBLE, index_rank, which_row + i, MPI_COMM_WORLD, &request[request_counter]);
            request_counter++;
            index_rank++;
        }
        index_rank = rank + 1;
    }

    // if(rank == nprocs-1){
    //     MPI_Request request2[nprocs-1];
    //     int isFWSBSdone2 = 1;
    //     for(int x = 0; x<nprocs-1; x++){
    //         MPI_Isend(&isFWSBSdone2, 1, MPI_INT, x, 0, MPI_COMM_WORLD, &request2[x]);
    //     }
    // }
    // MPI_Ibcast(&isFWSBSdone, 1, MPI_INT, nprocs-1, MPI_COMM_WORLD, isFWSBSdone_request);

}

// Backward substitution: Solves L^T * s = r
void backward_substitution_mpi(  const std::vector<std::vector<int>> &which_row_which_rank,
                                // const std::vector<int> &L_rowptrs,
                                const std::vector<double> &sub_L_values,
                                const std::vector<int> &sub_L_columns,
                                const std::vector<int> &sub_L_rowptrs,
                                const std::vector<double> &sub_r,
                                std::vector<double> &sub_s,
                                MPI_Request *request) {
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    sub_s.resize(sub_r.size(), 0.0);

    int which_row = which_row_which_rank[rank][0];
    // for(int i = 0; i<L_rowptrs.size(); i++){
    //     if(L_rowptrs[i] == sub_L_rowptrs.front()){
    //         which_row = i;
    //         break;
    //     }
    // }
    // for (int i = 0; i < n; ++i) {
    //     double sum = 0.0;
    //     if (i % size == rank) { // Row assigned to this process
    //         for (int j = 0; j < i; ++j)
    //             sum += L[i][j] * x[j];

    //         x[i] = (b[i] - sum) / L[i][i];
    //     }

    //     // Broadcast the computed value of x[i] to all processes
    //     MPI_Bcast(&x[i], 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
    // }


    // MPI_Request request[(nprocs-1-rank)*(sub_L_rowptrs.size()-1)];
    int request_counter = 0;
    ////////////
    int start, end;
    int index_rank = rank - 1;
    double sub_L_value_current_diagonal = 0;
    double received_s_value = 0.0;
    int source_rank;



    // if(rank == 1){
    //     std::cout << "rank: " << rank << std::endl;
    //     std::cout << "sub_L_rowptrs:" << std::endl;
    //     for (int ptr : sub_L_rowptrs) {
    //         std::cout << ptr << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

    // if(rank == 1){
    //     std::cout << "rank: " << rank << std::endl;
    //     std::cout << "sub_L_columns:" << std::endl;
    //     for (int ptr : sub_L_columns) {
    //         std::cout << ptr << std::endl;
    //     }
    //     std::cout << std::endl;
    // }    
    
    for (int i = sub_L_rowptrs.size()-2; i >= 0; i--) { // Per row iteration
        start = sub_L_rowptrs[i]-sub_L_rowptrs.front();
        end = sub_L_rowptrs[i+1]-sub_L_rowptrs.front()-1;

        if(sub_L_columns[start] == which_row + i){
            sub_s[i] += sub_r[i]/sub_L_values[start];
            sub_L_value_current_diagonal = sub_L_values[start];
            // if(rank == nprocs-1){
            //     std::cout << "above RANK: " << rank << "//// sub_s[" << i << "]: " << sub_s[i] << std::endl;
            //     std::cout << "above RANK: " << rank << "//// start: " << start << std::endl; 
            //     std::cout << "above RANK: " << rank << "//// end: " << end << std::endl;  
            // }   
        }
        else{
            std::cout << "HATATAA" << std::endl;
        }

        start++;
        for(; start <= end; start++){ //per element in a row
            // if(rank == 1){
            //     std::cout << "RANK: " << rank << "//// sub_L_columns[end]: " << sub_L_columns[end] << std::endl;
            //     std::cout << "RANK: " << rank << "//// which_row: " << which_row << std::endl; 
            // }
            if(which_row+sub_L_rowptrs.size()-2 >= sub_L_columns[start]){
                // if(rank == 1){
                //     std::cout << "RANK: " << rank << "//// sub_s[" << i << "]: " << sub_s[i] << std::endl; 
                // }

                sub_s[i] -= sub_L_values[start]/sub_L_value_current_diagonal*sub_s[sub_L_columns[start]-which_row];

                // if(rank == 1){
                //     std::cout << "RANK: " << rank << "//// sub_s[" << i << "]: " << sub_s[i] << std::endl; 
                // }            
            }
            else{
                source_rank = -1;
                for(int x = 0; x < rank; x++){
                    if(which_row_which_rank[x][0] <= sub_L_columns[start] && sub_L_columns[start] <= which_row_which_rank[x][1]){
                        source_rank = x;
                        break;
                    }
                }

                MPI_Recv(&received_s_value, 1, MPI_DOUBLE, source_rank, sub_L_columns[start], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(rank == 1){
                    std::cout << "RANK: " << rank << "//// received_s_value: " << received_s_value << "which_row: " << which_row << " i: " << i <<  std::endl; 
                }  
                // std::cout << "RANK: " << rank << "//// received_s_value: " << received_s_value << std::endl; 
                sub_s[i] -= sub_L_values[start]/sub_L_value_current_diagonal*received_s_value;
            }
        }

        // if(rank == 0){
        //     std::cout << "RANK: " << rank << "//// sub_s[" << i << "]: " << sub_s[i] << std::endl; 
        // }
        // empty_entire_s[which_row+i] = sub_s[i];

        while(index_rank >= 0){
            MPI_Isend(&sub_s[i], 1, MPI_DOUBLE, index_rank, which_row + i, MPI_COMM_WORLD, &request[request_counter]);
            request_counter++;
            index_rank--;
        }
        index_rank = rank - 1;
    }

    // if(rank == nprocs-1){
    //     MPI_Request request2[nprocs-1];
    //     int isFWSBSdone2 = 1;
    //     for(int x = 0; x<nprocs-1; x++){
    //         MPI_Isend(&isFWSBSdone2, 1, MPI_INT, x, 0, MPI_COMM_WORLD, &request2[x]);
    //     }
    // }
    // MPI_Ibcast(&isFWSBSdone, 1, MPI_INT, nprocs-1, MPI_COMM_WORLD, isFWSBSdone_request);

}