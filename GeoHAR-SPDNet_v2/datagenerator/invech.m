function A = invech(v)
    % INVECH build a symmetric matrix starting from its vech representation
    % (v)
    %
    % A = INVECH(v) returns a symmetric square matrix (A) of size (n x n)
    % where n is automatically determined
    % where the lower triangular is given from vector v
    n = floor((sqrt(1 + 8 * width(v)) - 1) / 2);
    A = zeros(n, n);
    
    % creates logical matrix (indices) for lower triangular matrix
    idx = tril(true(n));
    
    % assigns v values based on indices
    A(idx) = table2array(v);
    
    % completes the matrix
    A = A + tril(A, -1)';
end