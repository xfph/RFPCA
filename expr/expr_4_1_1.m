% expr4_1_1 performs the convergence of ECME and PX-ECME algorithms.
clc; clear all;
addpath('./data'); addpath('../prog'); addpath('../prog/common'); 

load simu1-data1.mat;
%load simu1-data2.mat;

[dc, dr, N] = size(X);
alg={'ECME' 'PX-ECME'};
iternum=zeros(1, length(alg)); logL=iternum;
disp_ini=1; disp_it=1;
if disp_ini
    fprintf('    data: %d points in [%d %d] dimensions\n',N,[dc dr]);
    fprintf('\n--> Fitting RFPCA model by different ML algorithms:\n');
end

if ~disp_it fprintf('\n\t\t\t\t\t\t CPU time\t\t Iterations\t\t\t logL'); end
for i=1:length(alg)
    if ~disp_it fprintf(['\n--------->%6s:'], alg{i});
    else fprintf(['\n--------->%6s:\n'], alg{i}); end
    opts = []; opts.alg = alg{i}; opts.maxit = 1000; opts.disp_it = disp_it; opts.tol = 1e-8;
    [bpt, Xc, opts] = mvt_ini(X, opts);
    [bp{i},opts] = mvt(bpt,X,opts);
    T{i} = opts.time.ini + opts.time.preit + cumsum(opts.time.it);
    T1(i) = T{i}(end);
    iternum(i) = opts.itnum;
    logL(i) = opts.logL;
    ll{i} = opts.errlog;
    if disp_it fprintf('\n\t\t CPU time\t\t Iterations\t\t\t logL\n'); end
    fprintf('%15.3f', T{i}(end)); fprintf('%16d', iternum(i)); fprintf('%22.5f', logL(i));
end
% eval(['save simu1_low.mat' ' alg T T1 iternum ll;']);
% eval(['save simu1_high.mat' ' alg T T1 iternum ll;']);
fprintf('\n'); rmpath('../prog'); rmpath('../prog/common'); 