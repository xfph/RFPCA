function [bp,opts] = mxvt(bp,X,opts)
%   mxvt performs maximum likelihood estimation for T distribution
%
%   bp:      Structure of all parameters in T distribution
%             bp.M: mean matrix;
%             bp.nu: degree of freedom;
%             bp.W: weight matrix;
%             bp.Sc: Sigma_c;
%             bp.Sr: Sigma_r;
%  
%  opts:      Structure of related settings in T distribution
%             opts.alg: 'ECME'
%             opts.errlog:  a vector recording the log likelihood at each iteration.
%	          opts.disp_it = 1 display error values; logs error;
%	          opts.tol: threshold (relative change in log-likelihood).
%	          opts.maxit: the maximum number of iterations; default 100.
%             opts.logL: final log-likelihood.
%             opts.itnum: number of iterations
%             opts.time: CPU time
%
%  reference: [1] Thompson G Z, Maitra R, Meeker W Q, et al. Classification with the matrix-variate-t distribution[J]. 

bp.model = 'TPCA';
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

nu1=2.01; nu2=1000; eta = 1e-6;
[d(1), d(2), N] = size(X);
d = bp.d; sd = sum(d); pd =prod(d);
opts.time.preit = toc;

t2 = 0; bp.W = zeros(d(1),d(1),N); 
for n = 1:niters
    tic
    opts.itnum = n;
    Xc = X - bp.M;
    [bp.lmd{1},bp.U{1}] = eigdec2(bp.Sc,d(1),eta); Ci = bp.U{1}.*(bp.lmd{1}.^(-.5))'; Sci = Ci*Ci';
    [bp.lmd{2},bp.U{2}] = eigdec2(bp.Sr,d(2),eta); Ri = bp.U{2}.*(bp.lmd{2}.^(-.5))'; 
    
    % E-step: E[W]
    Slogdet_W = 0;  SWX = zeros(d(1),d(2));
    SXTWX = zeros(d(2),d(2)); Slogdet_juli = 0; 
    for i = 1:N
        XR = Xc(:,:,i)*Ri; 
        XSrXT = XR*XR'; ScXSrXT = Sci*XSrXT;  
        bp.W(:,:,i) = (bp.nu+sd-1)*inv(XSrXT+bp.Sc);
        WX = bp.W(:,:,i)*X(:,:,i);
        SWX = SWX+WX;  
        SXTWX = SXTWX+X(:,:,i)'*WX;
        Slogdet_W = Slogdet_W+logdet2(bp.W(:,:,i),'chol');
        Slogdet_juli = Slogdet_juli+logdet2(eye(d(1))+ScXSrXT);
    end
    SW = sum(bp.W,3); 
    mSW = SW/N; V = chol(mSW); Vi = inv(V); mSWi = Vi*Vi'; 
    logdet_SW = logdet2(SW,'chol');
    mSWX = SWX/N; mSXTWX = SXTWX/N;
    
    % update nu
    limit=[nu1,nu2];   
    g1 = g(limit(1),d(1),d(2),N,Slogdet_W,logdet_SW);
    g2 = g(limit(2),d(1),d(2),N,Slogdet_W,logdet_SW);
    if g1*g2 < 0
        bp.nu = fzero(@(x) g(x,d(1),d(2),N,Slogdet_W,logdet_SW),limit);
    elseif g1 < 0
        bp.nu = nu1;
    else
        bp.nu = nu2;
    end
    bp.nu = max(bp.nu,nu1);
    
    if (disp_it || store || test)
        % the true log-likelihood
        logdet_Sc = sum(log(bp.lmd{1})); logdet_Sr = sum(log(bp.lmd{2}));
        ll1= -N/2 * (d(2) * logdet_Sc + d(1) * logdet_Sr + pd * log(pi));
        ll2 = N * ( mvgammaln((bp.nu+sd-1)/2,d(1)) - mvgammaln((bp.nu+d(1)-1)/2,d(1)) );
        ll3 = Slogdet_juli*(bp.nu+sd-1)/2;
        e = ll1 + ll2 - ll3;
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
    % update M
    bp.M = mSWi*mSWX;
    
    % update Sigma_r
    bp.Sr = (mSXTWX - mSWX'*mSWi*mSWX)/d(1);
    bp.Sr = (bp.Sr+bp.Sr')/2;
    
    % update Sigma_c
    bp.Sc = (bp.nu+d(1)-1)*mSWi;
    t2 = toc;
end
opts.logL = e;
opts.errlog = opts.errlog(1:opts.itnum);
opts.time.it = opts.time.it(1:opts.itnum);
end

function g = g(x,dc,dr,N,Slogdet_W,logdet_SW)
g = -mvdigamma((x+dc-1)/2,dc) + mvdigamma((x+dc+dr-1)/2,dc) + Slogdet_W/N - dc*log(x+dc+dr-1) - logdet_SW + dc*log(N*(x+dc-1));
end
