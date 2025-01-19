matfilename = 'Trefethen_20000.mat';
load(matfilename)


matrix = Problem.A;
rng(1,"twister");
b = -10 + (10-(-10)) .* rand(size(matrix,1),1);

[j,i,v] = find(sparse(matrix)');

row_ptrs_CSR = zeros(size(matrix, 1)+1, 1);

whichrow = 1;
temp = i(1);
row_ptrs_CSR(whichrow) = temp;
for x = 2:length(i)
    if i(x) ~= temp
        whichrow = whichrow + 1;
        temp = i(x);
        row_ptrs_CSR(whichrow) = x;
    end
end
row_ptrs_CSR(whichrow+1) = length(i)+1;

col_ind_CSR = j;
val_CSR = v;

% making the start index begin from 0 not 1 as in the MATLAB
row_ptrs_CSR = row_ptrs_CSR - 1;
col_ind_CSR = col_ind_CSR - 1;

% Specify the output file name
[~, matfilename, ~] = fileparts(matfilename);
filename_r = strcat('row_ptrs_', matfilename, '.txt');
filename_c = strcat('col_inds_', matfilename, '.txt');
filename_v = strcat('values_', matfilename, '.txt');
filename_b = strcat('b_', matfilename, '.txt');

% Open the file for writing
f_row_ptr = fopen(filename_r, 'w');
f_col_ind = fopen(filename_c, 'w');
f_val = fopen(filename_v, 'w');
f_b = fopen(filename_b, 'w');

% Check if the file is opened successfully
if f_row_ptr == -1 || f_col_ind == -1 || f_val == -1 || f_b == -1
    error('Could not open the file for writing.');
end

% Write the CSR vectors to the files
for x = 1:length(row_ptrs_CSR)
    if x ~= length(row_ptrs_CSR)
        fprintf(f_row_ptr, '%d\n', row_ptrs_CSR(x));
    else
        fprintf(f_row_ptr, '%d', row_ptrs_CSR(x));
    end
end

for x = 1:length(col_ind_CSR)
    if x ~= length(col_ind_CSR)
        fprintf(f_col_ind, '%d\n', col_ind_CSR(x));
    else
        fprintf(f_col_ind, '%d', col_ind_CSR(x));
    end
    
end

for x = 1:length(val_CSR)
    if x ~= length(val_CSR) 
        fprintf(f_val, '%.10g\n', val_CSR(x));
    else
        fprintf(f_val, '%.10g', val_CSR(x));
    end
end

for x = 1:length(b)
    if x ~= length(b) 
        fprintf(f_b, '%.10g\n', b(x));
    else
        fprintf(f_b, '%.10g', b(x));
    end
end

% Close the file
fclose(f_row_ptr);
fclose(f_col_ind);
fclose(f_val);
fclose(f_b);
