function [bp,opts] = mt(bp,x,opts)
%   mt performs maximum likelihood estimation for multivariate t distribution
%
%   bp:      Structure of all parameters in multivariate t distribution
%             bp.M: mean matrix;
%             bp.nu: degree of freedom;
%             bp.tau: Weights;
%             bp.S: Sigma;
%  
%  opts:      Structure of related settings in multivariate t  distribution
%             opts.alg: 'PX-ECME'
%             opts.errlog:  a vector recording the log likelihood at each iteration.
%	          opts.disp_it = 1 display error values; logs error;
%	          opts.tol: threshold (relative change in log-likelihood).
%	          opts.maxit: the maximum number of iterations; default 100.
%             opts.logL: final log-likelihood.
%             opts.itnum: number of iterations
%             opts.time: CPU time

bp.model = 'tPCA';
tic
niters = 100;
if isfield(opts, 'maxit') 
    niters = opts.maxit;
end
disp_it = 0;
if isfield(opts, 'disp_it') 
    disp_it = opts.disp_it;
end
store = 0;
if (nargout > 1)
    store = 1;
    opts.errlog = zeros(1, niters);
    opts.time.it = zeros(1, niters);
end
test = 0;
if isfield(opts, 'tol') 
    test = 1;
end

nu1 = 2; nu2 = 1000; 
[d,N]=size(x);  Id=ones(d,1);
opts.time.preit = toc;

t2=0;
for n = 1:niters
    tic
    opts.itnum = n;
    [bp.lmd,bp.U] = eigdec2(bp.S,d,1e-6);
    C = bp.U*diag(bp.lmd.^(-0.5));
    xC = C'*x;  tr1 = sum(xC.^2,1);
    MTC = bp.M'*C; tr2 = MTC*xC; 
    tr3 = sum(sum(MTC.^2));
    tr = tr1 - 2*tr2 + tr3; 
    
    % update nu
    limit = [nu1,nu2];
    f1 = f(limit(1),d,N,tr);
    f2 = f(limit(2),d,N,tr);
    if f1*f2 < 0
        bp.nu = fzero(@(x) f(x,d,N,tr),limit);
    elseif f1 < 0
        bp.nu = nu1;
    else
        bp.nu = nu2;
    end
    
    if (disp_it || store || test)
        % the actual log-likelihood
        logdet_S = sum(log(bp.lmd));
        ll1 = -N/2*(logdet_S+d*log(pi*bp.nu))+N*(gammaln((d+bp.nu)/2)-gammaln(bp.nu/2));
        ll2 = sum(log(tr./bp.nu+1))*((d+bp.nu)/2);
        e = ll1 - ll2;
        t1 = toc;
        opts.time.it(n) = t1 + t2;
        if store
            opts.errlog(n) = e;
        end
        if (disp_it > 0)
            if n>1
                fprintf(1, 'Iter %4d  logL %11.6f, relative increment %e\n', n, e/N, (e-eold)/abs(eold));
            else
                fprintf(1, 'Iter %4d  logL %11.6f\n', n, e/N);
            end
            if (n > 1 && e<eold)
                fprintf('----> LogL decreased in iteration %4d\n', n);
            end
        end
        if test
            if (n > 1 && abs(e-eold) / abs(e) < opts.tol)
                opts.logL = e;
                opts.errlog = opts.errlog(1:opts.itnum);
                opts.time.it = opts.time.it(1:opts.itnum);
                return;
            else
                eold = e;
            end
        end
    end
    
    tic
    % E-step: E[tau]
    bp.tau = (bp.nu+d)./(bp.nu+tr); 
    sumtx = x*bp.tau';
    srt = bp.tau.^.5;  st = sum(bp.tau); 
    % update M
    bp.M = sumtx/st;
    
    % update sigma
    tx = x.*kron(srt,Id); 
    stxx = tx*tx'; 
    tMtM = bp.M*bp.M'*st;
    %bp.S = (stxx-tMtM)/N;
    bp.S = (stxx-tMtM)/st;
end
end

function f=f(x,d,N,delta)
f=1-digamma(x/2)+log(x/2)+digamma((d+x)/2)-log((d+x)/2)+(sum(log((x+d)./(x+delta)))-sum((x+d)./(x+delta)))/N;
end
