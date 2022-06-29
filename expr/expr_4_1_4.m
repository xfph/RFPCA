% expr_4_1_4 performs outlier detection.
clc; clear all;
addpath('./data'); addpath('../prog'); addpath('../prog/common'); 

% case1: U(100,110); case2: U(100,102); case3: U(100000,100002)

rng(5);
model = {'RFPCA','tPCA'}; 
datatype='mx_g'; otypes='oc';
d = [4,10]; q = [1,3]; ndata1=1000; p = 0.05;
a = [100,110; 100,102; 100000,100002]; 
cj = 2;

load simu4-data3-2.mat;
bp = simu4(X,model);
%eval(['save simu4_caseoc' int2str(cj) '.mat bp']);
rmpath('../prog'); rmpath('../prog/common'); 

function bp = simu4(X,model)
[d(1),d(2),N]=size(X); pd = prod(d);
x = reshape(X,[pd,N]);

for i=1:length(model)
    opts = []; opts.maxit = 1000; opts.tol = 1e-8; opts.disp_it = 1;
    switch model{i}
        case 'RFPCA'
            [bpt, Xc, opts] = mvt_ini(X, opts);  
            [bp{i},opts] = mvt(bpt,X,opts);
        case 'tPCA'
            opts.ini = 'adhoc';
            [bpt, xc, opts] = mt_ini(x, opts); 
            [bp{i},opts] = mt(bpt,x,opts);
    end
    bp{i}.model = model{i};
end
end