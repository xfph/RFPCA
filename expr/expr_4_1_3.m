%  expr_4_1_3 used to obtain the relative difference between the estimated
%  and true covariance matrix for different proportions of outliers.
clc; clear all;
addpath('./data'); addpath('../prog'); addpath('../prog/common');

model = {'RFPCA','TPCA','FPCA','BPCA','tPCA','PCA'};

datatype='mx_g'; ndata1=1000; nu=[];
d = [4,10]; q = [1,3];
a = [-100,100;-10000,10000;100,110;10000,11000];
caseid = {'a','b','c','d'}; cj = 3;
p = [0,0.01,0.02,0.03,0.07,0.09];
% rep = 50;
otypes = 'oc'; % 'pc' for PC outliers, 'oc' for OC outliers and 'all' for
% PC+OC outliers.

fprintf(['\n In the presence of ',otypes,' outliers: \n']);
for pj = [3,4]
    fprintf('\n Proportion of outliers: %s',string(p(pj)));
    load(['simu3-Sit-III-p',int2str(pj),'.mat']);
    [err,bp]=simu3(X,model,mo);
    fprintf('\n the difference between the true principal subspace and the estimated principal subspace: \n');
    err
    %eval(['save relED_50_' int2str(pj) '_case_' caseid{cj} '_' otypes '.mat err methods;']);
end
rmpath('./data'); rmpath('../prog'); rmpath('../prog/common');

function [err,bp]=simu3(X,model,mo)
[d(1),d(2),N] = size(X); pd = prod(d);
x = reshape(X,[pd,N]);

n_model = length(model);
err=zeros(1, n_model); bp = cell(1,n_model);
for i=1:n_model
    opts = []; opts.maxit = 1000; opts.tol = 1e-8; opts.disp_it = 0;
    switch model{i}
        case 'RFPCA'
            [bpt, ~, opts] = mvt_ini(X, opts);
            [bp{i},~] = mvt(bpt,X,opts);
        case 'TPCA'
            [bpt, ~, opts] = mxvt_ini(X, opts);
            [bp{i},~] = mxvt(bpt,X,opts);
        case 'FPCA'
            [bpt, ~, opts] = mvn_ini(X, opts);
            [bp{i},~] = mvn(bpt,X,opts);
        case 'BPCA'
            bp{i} = bpca(X);
        case 'tPCA'
            [bpt, ~, opts] = mt_ini(x, opts);
            [bp{i},~] = mt(bpt,x,opts);
        case 'PCA'
            opts.center=0; bp{i}.M = mean(x,2); xc = x-bp{i}.M;
            bpt = pca(xc',opts);
            bp{i}.lmd = bpt.lmd; bp{i}.U = bpt.U;
    end
    bp{i}.model = model{i};
end

for i=1:n_model
    switch bp{i}.model
        case {'RFPCA','FPCA'}
            bp{i}.S = kron(bp{i}.Sr,bp{i}.Sc);
        case 'BPCA'
            bp{i}.S = kron(bp{i}.Sr,bp{i}.Sc)/trace(bp{i}.Sc);
        case 'TPCA'
            bp{i}.S = kron(bp{i}.Sr,bp{i}.Sc)/bp{i}.nu;
        case 'PCA'
            bp{i}.S = bp{i}.U.*bp{i}.lmd'*bp{i}.U';
    end
    err(i) = norm(mo.S-bp{i}.S,'fro')/norm(mo.S,'fro');
end
end