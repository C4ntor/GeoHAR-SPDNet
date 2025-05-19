function Y = vl_logeloss(X, c, dzdy)
% Log-Euclidean Metric/Loss function

batchSize = length(c);

n = size(X{1}, 1);  % Assuming square matrices
if nargin < 3
    Y = 0;
    for i = 1:batchSize
        diff = logm(X{i}) - logm(c{i});
        logEuclideanLoss = sum(diff(:).^2) / (n^2);
        Y = Y + logEuclideanLoss;
    end
    Y = Y / batchSize;
else
    Y = cell(1, batchSize);  % Initialize cell array for gradients
    for i = 1:batchSize
        grad = 2*(logm(X{i}) - logm(c{i}))*inv(X{i}) / (n^2);
        Y{i} = grad * dzdy;
    end
end
end

