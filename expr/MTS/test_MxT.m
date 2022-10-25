% This code is used to show the results of running the MxT with parsimonious covariance structures on the AUSLAN dataset.
%
% To run this code, you MUST：
% 1)Change 'Rpath' in the code to the installation path of R on your computer.
% 2)Install R packages 'R.matlab' and 'MixMatrix'. To do this, you can run 'install.packages('R.matlab')' 
% and 'install.packages('MixMatrix')' in R, and specify which CRAN mirror to use.
%
% Note: this code shows all the six cases by default, which takes about 3 minutes. If you only want to run some of these cases, 
% set 'dis_case' in the code to the corresponding case numbers. For example, if you only run cases I and IV, dis_case=[1,4].
clc;clear all;
addpath('../../prog/common');
dataname = 'AUSLAN'; dc=47; dr=22;
covstr = {'Case I: Sigma_c: AR(1) and Sigma_r: U','Case II: Sigma_c: U and Sigma_r: AR(1)',...
    'Case III: Sigma_c: AR(1) and Sigma_r: AR(1)','Case IV: Sigma_c: CS and Sigma_r: U',...
    'Case V: Sigma_c: U and Sigma_r: CS','Case VI: Sigma_c: CS and Sigma_r: CS'}; 
% U denotes the unconstrained covariance structure.

dis_case = 1:6; % Select the cases you want to test. 
Rpath = 'D:\Program Files\R\R-4.0.3\bin\x64'; % The path for the installed 'R.exe'.

load AUS_without.mat; % Load data
fprintf('Data: AUSLAN\t\t Method: MxT\n');
for cj = dis_case
    fprintf(['---> ',covstr{cj},':\n']);
    RscriptFileName = ['.\runMxT_case',int2str(cj),'.R'];
    RunRcode(RscriptFileName,Rpath);
    
    eval(['load res_case' int2str(cj) '.mat;']);
    if sum(sum(isnan(Sc))) && sum(sum(isnan(Sr)))
        fprintf('\tMxT failed to run due to numerical problems.\n')
    else
        bp.model = 'TPCA';
        [bp.lmd{1},bp.U{1}]=eigdec(Sc,dc);
        [bp.lmd{2},bp.U{2}]=eigdec(Sr,dr);
        err = nnerr2(bp,X(:,:,1:length(X_label)),Y,X_label,Y_label);
        fprintf('\tRun successfully!\n')
        fprintf('\tClassifier: 1-nearest neighbor classifier\n')
        fprintf('\t%s%3.1f\n','The error rate of MxT: ',err*100);
    end
end
rmpath('../../prog/common');

function err = nnerr2(bp,X,Y,X_label,Y_label)
N_test=size(Y,3); N_train=length(X_label);
Z_tr=compr(bp,X); Z_te=compr(bp,Y);
Dtmp=zeros(1,N_train); I=zeros(1,N_test);
for i=1:N_test
    D=Dtmp;
    for j=1:N_train
        Z=Z_te(:, :, i)-Z_tr(:, :, j);
        D(j)=sum(sum(Z.^2, 1), 2);
    end
    [~,index]=min(D);
    I(i)=index;
end
err=sum(X_label(I)~=Y_label)/N_test;
end