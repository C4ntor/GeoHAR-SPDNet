function Y = vl_mseloss(X,c,dzdy)
% Softmax layer
batchSize = length(c);
n=size(X{1},1);
%n = sz(1)*sz(2) ;
if nargin < 3
   Y=0;
   for i = 1:batchSize
            % Compute the squared difference between prediction and ground truth
            diff = X{i} - c{i};
            % Sum element-wise squared differences for the current pair
            mseLoss = sum(diff(:).^2) / (n^2);  % Average over elements in the matrix
            % Accumulate the total loss
            Y = Y + mseLoss;
   end
  Y = Y / batchSize;
else
  Y = cell(1, batchSize);  % Initialize cell array to store gradients
        for i = 1:batchSize
            % Compute the gradient for the current pair
            grad = (X{i} - c{i}) / (n^2);  % Average over elements in the matrix
            % Multiply by dzdy to propagate through the chain rule
            Y{i} = grad * dzdy;
        end
  %disp('loss_grad (Y) = dzdx')
  %disp(size(Y))
  %disp(Y{1})
  %Y = single(cat(3, Y{:}));
  %disp(size(ressss))
  %disp(ressss)
end

