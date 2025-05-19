function block_diag_matrix = diagblockhar(varargin)
    % DIAGBLOCKHAR Build a diagonal block input matrix (A) using HAR structure
    % varargin - 22 SPD matrices and an optional method for Fréchet mean computation

    % Default to Log-Euclidean if no method is specified
    if ischar(varargin{end}) || isstring(varargin{end})
        method = varargin{end}; % Extract the method
        matrices = varargin(1:end-1); % Exclude the method from matrices
    else
        method = 'log-euclidean';
        matrices = varargin;
    end

    % Compute Fréchet means using the specified method
    frecw = frechet_mean(matrices{end-4:end}, method);
    frecm = frechet_mean(matrices{:}, method);

    % Build the block diagonal matrix
    block_diag_matrix = blkdiag(matrices{22}, frecw, frecm);
end
