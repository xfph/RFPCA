function [bp,opts]=mvt(bp,X,opts)
%   mvt performs maximum likelihood estimation for matrix-variate t distribution
%
%   bp:      Structure of all parameters in matrix-variate t distribution
%             bp.M: mean matrix;
%             bp.nu: degree of freedom;
%             bp.tau: Weights;
%             bp.Sc: Sigma_c;
%             bp.Sr: Sigma_r;
%  
%  opts:      Structure of related settings in matrix-variate t distribution
%             opts.alg: 'ECME, 'PX-ECME'
%             opts.errlog:  a vector recording the log likelihood at each iteration.
%	          opts.disp_it = 1 display error values; logs error;
%	          opts.tol: threshold (relative change in log-likelihood).
%	          opts.maxit: the maximum number of iterations; default 100.
%             opts.logL: final log-likelihood.
%             opts.itnum: number of iterations
%             opts.time: CPU time

bp.model = 'RFPCA';
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

nu1=2.01; nu2=1000; 
[d(1), d(2), N] = size(X);
d = bp.d; pd = prod(d);
Ic = ones(1,d(1)); Ir = ones(1,d(2));
data.X = X; data = cdas(data.X);
opts.time.preit = toc;

t2=0; eta=1e-6;
for n = 1:niters
    tic
    opts.itnum = n;
    [bp.lmd{1},bp.U{1}]=eigdec2(bp.Sc,d(1),eta); Ci=bp.U{1}.*(bp.lmd{1}.^(-.5))'; %Ci=bp.U{1}*diag(bp.lmd{1})^(-0.5); 
    [bp.lmd{2},bp.U{2}]=eigdec2(bp.Sr,d(2),eta); Ri=bp.U{2}.*(bp.lmd{2}.^(-.5))'; %Ri=bp.U{2}*diag(bp.lmd{2})^(-0.5);
    CX = reshape(permute(reshape(Ci'*data.X1,[d(1) d(2) N]),[1 3 2]),[d(1)*N d(2)]);
    CXR = reshape(permute(reshape(CX*Ri,[d(1) N d(2)]),[1 3 2]),[pd N]);
    tr1 = sum(CXR.^2,1);
    CMR = reshape(Ci'*bp.M*Ri,[1 pd]);
    tr2 = CMR*CXR;                                        
    tr3 = sum(sum(CMR.^2));
    tr = tr1 - 2*tr2 + tr3;
    
    % ECME¡¢PX-ECME: update nu
    limit=[nu1,nu2];
    g1 = g(limit(1),d(1),d(2),N,tr);
    g2 = g(limit(2),d(1),d(2),N,tr);
    if g1*g2 < 0
        bp.nu = fzero(@(x) g(x,d(1),d(2),N,tr),limit);
    elseif g1 < 0
        bp.nu = nu1;
    else
        bp.nu = nu2;
    end
  
    if (disp_it || store || test)
        logdet_Sc = sum(log(bp.lmd{1}));  logdet_Sr = sum(log(bp.lmd{2}));
        ll1= -N/2 * (d(2) * logdet_Sc + d(1) * logdet_Sr + pd * log(pi*bp.nu));
        ll2 = N * ( gammaln((pd+bp.nu)/2) - gammaln((bp.nu/2)) );
        tr0 = sum(log(tr./bp.nu+1))*((pd+bp.nu)/2);
        e = ll1 + ll2 - tr0;
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
    bp.tau = ((bp.nu+d(1)*d(2))./(bp.nu+tr)); 
        
    % update M
    sumtX = reshape(data.X1V*bp.tau',[d(1) d(2)]);
    srt = bp.tau.^.5;  st = sum(bp.tau); 
    bp.M = sumtX/st;
    
    % update sigma_C
    tX = data.X2 .* kron(srt,Ic); %[d(2) d(1)N]
    tXR = reshape(tX'*Ri,[d(1) N*d(2)]);
    tXRX = tXR * tXR';
    MR = bp.M * Ri;
    tMRM = MR * MR' * st;
    sumtXRX = tXRX - tMRM;
    switch opts.alg
        case 'ECME'
            bp.Sc = sumtXRX/(N*d(2));
        case 'PX-ECME'
            bp.Sc = sumtXRX/(d(2)*st);
    end
    
    % update sigma_R
    [bp.lmd{1},bp.U{1}]=eigdec2(bp.Sc,d(1),eta); Ci=bp.U{1}.*(bp.lmd{1}.^(-.5))'; %Ci=bp.U{1}*diag(bp.lmd{1})^(-0.5); 
    tXT = data.X1 .* kron(srt,Ir);
    tXC = reshape(tXT'*Ci,[d(2) N*d(1)]); tXCX = tXC * tXC';
    MC = bp.M' * Ci; tMCM = MC * MC' * st;
    sumtXCX = tXCX - tMCM;
    switch opts.alg
        case 'ECME'
            bp.Sr = sumtXCX/(N*d(1));
        case 'PX-ECME'
            bp.Sr = sumtXCX/(d(1)*st); 
    end
    t2=toc;
end
opts.logL = e;
opts.errlog = opts.errlog(1:opts.itnum);
opts.time.it = opts.time.it(1:opts.itnum);
end

function g = g(x,dc,dr,N,tr)
g = 1-digamma(x/2)+log(x/2)+digamma((dc*dr+x)/2)-log((dc*dr+x)/2)+(sum(log((x+dc*dr)./(x+tr)))-sum((x+dc*dr)./(x+tr)))/N;
end