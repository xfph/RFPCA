%  expr_4_1_5 performs the changes in CPU time for the different methods
%  with the sample size.
clc; clear all;
addpath('./data'); addpath('../prog'); addpath('../prog/common');
model = {'RFPCA','TPCA','FPCA'}; %'tPCA'

datatype='mx_g'; d = [100,100]; q = [3,3]; pd=prod(d); nu=[];
a = [100,110]; p = 0.005; N_tr = [200,500,1000,2000,4000,8000,13000];
otype = 'oc'; nj = 2;

rep = 1;
load simu5-data4-N2.mat;
T = zeros(length(N_tr),length(model)); itnum = T;
logL = zeros(1,length(model)); llh = cell(1,length(model)); T_it=llh;
for mj = 1:length(model)
    opts = []; opts.maxit = 1000; opts.tol = 1e-8; opts.disp_it = 1;
    switch model{mj}
        case 'RFPCA'
            t0=cputime;
            [bpt, Xc, opts] = mvt_ini(X_tr, opts);
            [bp{mj},opts] = mvt(bpt,Xc,opts);
        case 'TPCA'
            t0=cputime;
            [bpt, Xc, opts] = mxvt_ini(X_tr, opts);
            [bp{mj},opts] = mxvt(bpt,Xc,opts);
        case 'FPCA'
            t0=cputime;
            [bpt, Xc, opts] = mvn_ini(X_tr, opts);
            [bp{mj},opts] = mvn(bpt,Xc,opts);
        case 'tPCA'
            t0=cputime;
            [bpt, xc, opts] = mt_ini(x_tr, opts);
            [bp{mj},opts] = mt(bpt,xc,opts);
    end
    T(nj,mj) = cputime-t0; itnum(nj,mj) = opts.itnum;
    bp{mj}.model = model{mj};
    logL(mj) = opts.logL; llh{mj} = opts.errlog;
    T_it{mj} = opts.time.ini + opts.time.preit + cumsum(opts.time.it);
end
%     eval(['save simu5_t_itnum_' otype '.mat'  ' T itnum;']);
rmpath('./data');rmpath('../prog'); rmpath('../prog/common');