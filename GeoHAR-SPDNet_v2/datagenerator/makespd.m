function A_spd = makespd(A, epsilon)
    % This function takes a matrix A and approximates it to be SPD by:
    % 1. Symmetrizing the matrix.
    % 2. Clipping any negative eigenvalues.
    % 3. adding a small regularization term to ensure positive definiteness.
    
    if nargin < 2
        epsilon = 1e-6;  % Default small regularization value
    end

    %A = 0.5 * (A + A');  % Ensure symmetry
    
    [Q, Lambda] = eig(A);

    Lambda = diag(max(diag(Lambda), epsilon));  % Clip eigenvalues
    
    A_spd = Q * Lambda * Q';
    
end
