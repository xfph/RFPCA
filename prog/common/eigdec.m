function [evals, evec] = eigdec(X, N)
%EIGDEC	Sorted eigendecomposition
%
%	Description
%	 EVALS = EIGDEC(X, N computes the largest N eigenvalues of the
%	matrix X in descending order.  [EVALS, EVEC] = EIGDEC(X, N) also
%	computes the corresponding eigenvectors.
%	Copyright (c) Ian T Nabney (1996-2001)

if nargout == 1
   evals_only = logical(1);
else
   evals_only = logical(0);
end

if N ~= round(N) || N < 1 || N > size(X, 2)
   error('Number of PCs must be integer, >0, < dim');
end

% Find the eigenvalues of the data covariance matrix
if evals_only
   % Use eig function as always more efficient than eigs here
   temp_evals = eig(X);
else
   % Use eig function unless fraction of eigenvalues required is tiny
   if (N/size(X, 2)) > 0.01
      [temp_evec, temp_evals] = eig(X);
   else
      options.disp = 0;
      [temp_evec, temp_evals] = eigs(X, N, 'LM', options);
   end
   temp_evals = diag(temp_evals);
end

% Eigenvalues nearly always returned in descending order, but just
% to make sure.....
[evals,perm] = sort(-temp_evals);
evals = -evals(1:N);
if ~evals_only
   if evals == temp_evals(1:N)
      % Originals were in order
      evec = temp_evec(:, 1:N);
      return
   else
      % Need to reorder the eigenvectors
      for i=1:N
         evec(:,i) = temp_evec(:,perm(i));
      end
   end
end