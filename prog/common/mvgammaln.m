function res = mvgammaln(x,p)
% mvgammaln: the log of the multivariate gamma

y = x+(1-(1:p))/2;
res = log(pi)*p*(p-1)/4+sum(gammaln(y));
