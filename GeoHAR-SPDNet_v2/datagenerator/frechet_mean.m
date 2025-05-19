function M = frechet_mean(varargin)
    % Input:
    % varargin - SPD matrices passed as individual arguments
    % varargin{end} can contain the computation method as a string: 'procrustes' or 'log-euclidean'
    % Output:
    % M - the Fréchet mean of the SPD matrices

    % Default method
    if ischar(varargin{end}) || isstring(varargin{end})
        method = varargin{end}; % Extract the method
        varargin = varargin(1:end-1); % Exclude the method from matrices
    else
        method = 'log-euclidean'; % Default method
    end

    n = numel(varargin); % Number of matrices passed
    d = size(varargin{1}, 1); % Dimensionality of the SPD matrices

    switch lower(method)
        case 'log-euclidean'
            % Compute mean in the Log-Euclidean metric
            sumLog = zeros(d);
            for i = 1:n
                sumLog = sumLog + logm(makespd(varargin{i}));
            end
            M = expm(sumLog / n);

        case 'procrustes'
            % Compute mean using the Procrustes method
            M = varargin{1}; % Start with the first matrix
            tol = 1e-6; % Convergence tolerance
            maxIter = 100; % Maximum iterations
            iter = 0;
            while true
                sumProcrustes = zeros(d);
                for i = 1:n
                    [U, ~, V] = svd(M^(-0.5) * makespd(varargin{i}) * M^(-0.5));
                    R = U * V';
                    sumProcrustes = sumProcrustes + R;
                end
                M_next = sumProcrustes / n;
                if norm(M_next - M, 'fro') < tol || iter >= maxIter
                    break;
                end
                M = M_next;
                iter = iter + 1;
            end

        otherwise
            error('Unsupported method. Choose "log-euclidean" or "procrustes".');
    end
end


