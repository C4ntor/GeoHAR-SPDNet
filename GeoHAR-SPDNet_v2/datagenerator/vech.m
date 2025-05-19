function vec = vech(A)
% It computes the vectorized form of the symmetric input matrix (A)
% the vector contains in order elements from and below the main diagonal
% stacking by column.

[M,N] = size(A);
if (M == N)
    vec  = [];
    for ii=1:M
        vec = [vec; A(ii:end,ii)];
    end
else
end