function Y = vl_frobloss(X, c, dzdy)
% Frobenius Norm Metric/Loss function

batchSize = length(c);  
n = size(X{1}, 1);      

if nargin < 3
    % Loss computation (forward pass)
    Y = 0;
    for i = 1:batchSize
        frobeniusLoss = norm(X{i} - c{i}, 'fro');
        Y = Y + frobeniusLoss;
    end
    % Step 4: Average the loss over the batch
    Y = Y / batchSize;
else
    % Gradient computation (backward pass)
    Y = cell(1, batchSize);  % Initialize cell array for gradients
    for i = 1:batchSize
        grad = 2 * (X{i} - c{i}) / norm(X{i} - c{i}, 'fro');
        
        % Step 2: Multiply by dzdy to propagate through the chain rule
        Y{i} = grad * dzdy;
    end
end
end