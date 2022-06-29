function res = mvdigamma(x,p)
% mvdigamma: Multivariate Digamma Function
% The multivariate digamma function is the derivative of the log of the multivariate gamma function; 
% for p = 1 it is the same as the univariate digamma function.
% psi_p(a) = ∑psi(a+(1-i)/2)
% where psi is the univariate digamma function (the derivative of the log-gamma function).

y = x+(1-(1:p))/2;
res = sum(digamma(y));