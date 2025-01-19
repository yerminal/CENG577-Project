load("LFAT5.mat");
% n = 10; % Size of the matrix 
% A = randomSPD(n);
% rand_matr = rand(n);
% A = rand_matr' * rand(n)*-10;
% A(4,:) = A(5,:);
% L2 = ichol(Problem.A);
A = full(Problem.A);
% isZero = isZeroRowExist(A)

% symt = issymmetric(A)
eigenvalues = eig(A)
% b = rand(size(A,1), 1)*10;
% b = zeros(length(A(1,:)), 1);

load("LFAT5_b.mat"); %% loading b
x = zeros(length(A(1,:)), 1);

max_iteration = 10000;
Lxx = ichol_cus(A);
[Lxx2, Dxx] = incomplete_chol_with_thresholding(A, 1e-3, 1e-6);
Lxx2 = full(Lxx2);
Dxx = full(Dxx);
opts.type = 'ict'; % Incomplete Cholesky with thresholding
opts.droptol = 1e-3; % Drop tolerance
L = full(ichol(Problem.A, opts));

truth_x = A\b;
seq_cg_x = cg(A, b, x, max_iteration);
parcg_x = parcg(A, b, x, max_iteration);

% precision_n = 5;
% test = getNdecimal(seq_cg_x, precision_n);
simil_score_seq = similarity_percent_vector(seq_cg_x, truth_x)
simil_score_par = similarity_percent_vector(parcg_x, truth_x)
%test for substitution
% L\b
% forward_substitution(L, b)
% L'\b
% backward_substitution(L', b)

% spA = sparse(A)
% Sp_ichol(spA)

function [x] = cg(A, b, x, max_iteration)
    r = b - A * x; % residual vector
    p = r; % search direction vector
    rsold = r' * r; % residual square norm old

    for iters = 1:max_iteration
        Ap = A * p;
        alpha = rsold / (p' * Ap); % step size
        x = x + alpha * p; % update solution
        r = r - alpha * Ap; % update residual
        rsnew = r' * r; % residual square norm new
        
        if sqrt(rsnew) < 1e-10
            iters
            break;
        end
        
        beta = (rsnew / rsold); % Conjugate direction coeff.
        p = r + beta * p; % new search direction
        rsold = rsnew;
    end
    error = sqrt(rsnew)

end

function a = ichol_cus(a)
	n = size(a,1);

	for k = 1:n
        % a(k,k) = abs(sqrt(a(k,k))) * abs(a(k,k)) / a(k,k);
        a(k,k) = sqrt(a(k,k));
		for i = (k+1):n
		    if (a(i,k) ~= 0)
		        a(i,k) = a(i,k)/a(k,k);
            end
        end
		for j = (k+1):n
		    for i = j:n
		        if (a(i,j) ~= 0)
		            a(i,j) = a(i,j) - a(i,k)*a(j,k);
                end
            end
        end
    end

    for i = 1:n
        for j = i+1:n
            a(i,j) = 0;
        end
    end
    

end

function x = forward_substitution(L, b)
    % Solve Lx = b using forward substitution.
    % Inputs:
    %   L - Lower triangular matrix (n x n)
    %   b - Right-hand side vector (n x 1)
    % Output:
    %   x - Solution vector (n x 1)
    
    n = length(b);
    x = zeros(n, 1); % Initialize solution vector
    
    for i = 1:n
        x(i) = (b(i) - L(i, 1:i-1) * x(1:i-1)) / L(i, i);
    end
end

function x = backward_substitution(U, b)
    % Solve Ux = b using backward substitution.
    % Inputs:
    %   U - Upper triangular matrix (n x n)
    %   b - Right-hand side vector (n x 1)
    % Output:.21412122
    %   x - Solution vector (n x 1)
    
    n = length(b);
    x = zeros(n, 1); % Initialize solution vector
    
    for i = n:-1:1
        x(i) = (b(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
    end
end

function [x] = parcg(A, b, x, max_iteration)
    % L = ichol_cus(A);
    
    % opts.type = 'ict'; % Incomplete Cholesky with thresholding
    % opts.droptol = 1e-3; % Drop tolerance
    % L = full(ichol(sparse(A), opts));
    
    [L, ~] = incomplete_chol_with_thresholding(A, 1e-3, 1e-6);

    L_T = L';
    r = b - A * x; % residual vector
    s = forward_substitution(L, r);
    
    rho_old = 1;
    p_old = zeros(length(b), 1);
    alpha = 0;

    for iters = 1:max_iteration
        rho_new = s'*s;
        w = backward_substitution(L_T, s);
        beta = rho_new/rho_old;
        
        
        p_new = w + beta * p_old;
        q = A*p_new;

        gamma = p_new' * q;
        x = x + alpha*p_old;
        alpha = rho_new/gamma;
        r = r - alpha*q;

        r_norm = r' * r;
        s = forward_substitution(L, r);

        if sqrt(r_norm) < 1e-10
            x = x + alpha*p_new;
            iters
            break;
        end

        % Updates old new
        rho_old = rho_new;
        p_old = p_new;
    end
    error = sqrt(r_norm)
end

% function A=Sp_ichol(A)
% 	n=size(A,1);
% 	ncols=A(n).col;
%     c_end=0;
%     for col=1:ncols
%         is_next_col=0;
%         c_start=c_end+1;
%         for i=c_start:n
%             if A(i).col==col % in the current column (col):
%                 if A(i).col==A(i).row 
%                     A(i).val=sqrt(A(1i).val); % take the square root of the current column's diagonal element
%                     div=A(i).val;
%                 else
%                     A(i).val=A(i).val/div; % divide the other current column's elements by the square root of the diagonal element
%                 end
%             end
%             if A(i).col>col % in the next columns (col+1 ... ncols):
%                 if is_next_col==0
%                     c_end=i-1;
%                     is_next_col=1;
%                 end
%                 v1=0;
%                 v2=0;
%                 for j=c_start:c_end
%                     if A(j).col==col
%                         if A(j).row==A(i).row % search for current column's (col) elements A(j) whose row index is the same as current element's A(i) row index
%                             v1=A(j).val;
%                         end
%                         if A(j).row==A(i).col % search for current column's (col) elements A(j) whose row index is the same as current element's A(i) column index
%                             v2=A(j).val;
%                         end
%                         if v1~=0 && v2~=0 % if these elements exist in the current column (col), recalculate the current element A(i):
%                             A(i).val=A(i).val-v1*v2;
%                             break;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

function A = randomSPD(n) 
    % Create a random n x n matrix 
    R = randn(n); 
    % Create a diagonal matrix with positive entries 
    D = diag(rand(n, 1) + n); % Ensure positive diagonal entries 
    % Generate a symmetric positive definite matrix 
    A = R' * R + D; 
end

function A = getNdecimal(A, N)
    for i=1:size(A,1)
        for j=1:size(A,2)
            A(i,j) = str2num(sprintf(strcat('%.', num2str(N), 'f'), A(i,j)));
        end
    end
end

function isIt = isZeroRowExist(A)
    for i=1:size(A,1)
        isIt = all(A(i,:) == zeros(1,size(A,2)));
    end
end

function percent = similarity_percent_vector(A, B)
    Min_values = min([abs(A); abs(B)]);
    error = (A - B)./Min_values*100/length(B);
    error(B==0 & A==0) = 0;
    error(B==0 & A~=0) = A(B==0 & A~=0)./1e-2;
    error(abs(error)>100/length(B)) = 100/length(B);
    sum_error = sum(abs(error));
    percent = getNdecimal(100-sum_error, 2);
end

function [L, D] = incomplete_chol_with_thresholding(A, droptol, epsilon)
    % Ensure the matrix is sparse
    if ~issparse(A)
        A = sparse(A);
    end

    % Get the size of the matrix
    n = size(A, 1);

    % Initialize L and D
    L = sparse(n, n);
    D = sparse(n, n);

    % Loop through rows
    for i = 1:n
        for j = 1:i
            if i == j
                % Compute diagonal element
                sum_diag = 0;
                for k = 1:j-1
                    sum_diag = sum_diag + L(j, k)^2;
                end
                val = A(j, j) - sum_diag;
                if val > 0
                    L(j, j) = sqrt(val);
                else
                    L(j, j) = sqrt(epsilon);
                    D(j, j) = D(j, j) + epsilon;
                end
            else
                % Compute off-diagonal element
                if A(i, j) ~= 0
                    sum_off_diag = 0;
                    for k = 1:j-1
                        sum_off_diag = sum_off_diag + L(i, k) * L(j, k);
                    end
                    val = (A(i, j) - sum_off_diag) / L(j, j);
                    if abs(val) > droptol
                        L(i, j) = val;
                    end
                end
            end
        end
    end
end
