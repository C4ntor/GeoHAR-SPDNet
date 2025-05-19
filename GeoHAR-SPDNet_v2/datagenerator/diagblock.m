function block_diag_matrix = diagblock(varargin)
    % DIAGBLOCKS build a diagonal block input matrix (A) starting from a
    % collection of lagged RCOV matrices
    n = nargin; % Number of matrices passed
    inputs = {}
     for i = 1:n
        inputs{end+1} =makespd(varargin{i});
        
     end
    block_diag_matrix = blkdiag(inputs{:});
   
end
