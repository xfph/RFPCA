function T = cdas(X, label)
% convert data array to struct

T.X = X;
if nargin > 1
    T.label = label;
end
if length(size(X)) < 3
    warning('X is a 2-dimensional array');
end
[d(1), d(2), ndata] = size(X);
pd = prod(d);
T.X1V = reshape(T.X, [pd, ndata]);
T.X1 = T.X(:, :); 
T.X2 = permute(T.X, [2, 1, 3]); 
T.X2 = T.X2(:, :);
T.XV2 = sum(T.X1V.^2, 1);
