function is_spd = is_spd(A)
    % Check if the matrix is SPD
    if ~issymmetric(A)
        is_spd = false;
        return;
    end
    eigenvalues = eig(A);
    is_spd = all(eigenvalues > 0);
end
