function [bp,opts] = mvn(bp,X,opts)
%   mvn performs ML estimation for matrix normal distribution.
%
%   bp:      Structure of all parameters in matrix-variate t distribution
%            bp.M: mean matrix;
%            bp.Sc: sigma_c;
%            bp.Sr: sigma_r;
%  
%  opts:      Structure of related settings in matrix-variate normal distribution
%             opts.alg: 'CM'
%             opts.errlog:  a vector recording the log likelihood at each iteration.
%	          opts.disp_it = 1 display error values; logs error;
%	          opts.tol: threshold (relative change in log-likelihood).
%	          opts.maxit: the maximum number of iterations; default 100.
%             opts.logL: final log-likelihood.
%             opts.itnum: number of iterations
%             opts.time: CPU time

bp.model = 'FPCA';
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

[d(1), d(2), N] = size(X);
d = bp.d; bp.M = bp.M;
data.X = X; data = cdas(data.X);
eta = 1e-6; % 1e-6; 1e-2; 0.05
opts.time.preit = toc;

for n = 1:niters
    tic
    opts.itnum = n;
    % update Sigma_c
    [bp.lmd{2},bp.U{2}]=eigdec2(bp.Sr,d(2),eta); Ri=bp.U{2}*diag(bp.lmd{2})^(-0.5); 
    XR = reshape(data.X2'*Ri,[d(1),N,d(2)]); XR = XR(:,:); % [d_c d_r*N]
    sumXRXT = XR*XR';
    bp.Sc = sumXRXT/(N*d(2));
    
    % update Sigma_r
    [bp.lmd{1},bp.U{1}]=eigdec2(bp.Sc,d(1),eta); Ci=bp.U{1}*diag(bp.lmd{1})^(-0.5); 
    XTC = reshape(data.X1'*Ci,[d(2),N,d(1)]); XTC = XTC(:,:); % [d_r d_c*N]
    sumXTCX = XTC*XTC';
    bp.Sr = sumXTCX/(N*d(1));

    [bp.lmd{2},bp.U{2}]=eigdec2(bp.Sr,d(2),eta); Ri=bp.U{2}*diag(bp.lmd{2})^(-0.5); 
    SiR = Ri*Ri';
    sumtr = trace(sumXTCX*SiR);
    
    if (disp_it || store || test) 
        logdet_Sc = sum(log(bp.lmd{1}));  logdet_Sr = sum(log(bp.lmd{2}));
        e=-N/2*(d(1)*d(2)*log(2*pi)+d(2)*logdet_Sc+d(1)*logdet_Sr)-1/2*sumtr;
        opts.time.it(n) = toc;
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
            if (n > 1 && abs(e - eold)/abs(e) < opts.tol)
                opts.logL = e;
                opts.errlog = opts.errlog(1:opts.itnum);
                opts.time.it = opts.time.it(1:opts.itnum);
                return;
            else
                eold = e;
            end
        end
    end
end