% expr_4_1_2 performs the accuracies of estimators
clc; clear all;
addpath('./data'); addpath('../prog'); addpath('../prog/common'); 

rng(100);
model = {'RFPCA','TPCA','FPCA','tPCA'};

load simu2-data1.mat; dataid = 1; % data from Mt distribution
% load simu2-data3.mat; dataid = 2; % data from matrix-normal distribution
dis_true={'MVT','MVN'};

[d(1),d(2),N_tr] = size(X_tr); N_ts = size(X_ts,3);
x_ts = reshape(X_ts,[prod(d),N_ts]);
x_tr = reshape(X_tr,[prod(d),N_tr]);

llh_true = llh(mo,X_ts,dis_true{dataid}); % true test log-likelihood
llh_est1 = zeros(1,length(model));
for mj=1:length(model)
    %fprintf(['\n--------->%6s:\n'],mdname{mj});
    opts = [];  opts.maxit = 1000; opts.tol=1e-8; opts.disp_it = 0;
    switch model{mj}
        case 'RFPCA'
            [bpt, Xc, opts] = mvt_ini(X_tr, opts);
            [bp{mj},opts] = mvt(bpt,Xc,opts);
            llh_est1(mj)= llh(bp{mj},X_ts,'MVT');
        case 'TPCA'
            [bpt, Xc, opts] = mxvt_ini(X_tr, opts);
            [bp{mj},opts] = mxvt(bpt,Xc,opts);
            llh_est1(mj) = llh(bp{mj},X_ts,'T');
        case 'FPCA'
            [bpt, Xc, opts] = mvn_ini(X_tr, opts);
            [bp{mj},opts] = mvn(bpt,Xc,opts);
            llh_est1(mj) = llh(bp{mj},X_ts,'MVN');
        case 'tPCA'
            [bpt, xc, opts] = mt_ini(x_tr, opts);
            [bp{mj},opts] = mt(bpt,xc,opts);
            llh_est1(mj) = llh(bp{mj},x_ts,'mt');
    end
end
fprintf('test log-likelihood for each method:\n')
array2table(llh_est1,'VariableNames',model)
% eval(['save simu2_test_llh_' int2str(dataid) '.mat' ' llh_true llh_est;']);
rmpath('./data'); rmpath('../prog'); rmpath('../prog/common'); 
