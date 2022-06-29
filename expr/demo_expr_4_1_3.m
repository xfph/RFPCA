% Show some of the results in Table 2 of the main text.
clc;clear all;
mdname = {'RFPCA','TPCA','FPCA','BPCA','tPCA','PCA'};

p = [0,0.01,0.02,0.03,0.07,0.09];
a = [-100,100;-10000,10000;100,110;10000,11000];
caseid = {'a','b','c','d'}; % 'a' for Sit-I: U(-100,100); 'b' for Sit-II: U(-10000,10000); 
% 'c' for Sit-III: U(100,110); 'd' for Sit-IV: U(10000,11000);
casenames = {'U(-100,100)','U(-10000,10000)','U(100,110)','U(10000,11000)'};
showid = 3; % Presenting the results for Sit-III.
showpj = [2,4]; % Show the results when the outlier proportions are 1% and 3%.
otype = 'pc'; % 'pc' for PC outliers, 'oc' for OC outliers and 'all' for
% PC+OC outliers.

addpath('./result');
otype = 'pc'; % 'pc' for PC outliers
fprintf(['\n outliers type:\t',otype,'\n']);
for cj = showid
    fprintf(['\n outliers are generated from ',casenames{cj},'\n']);
    merr = zeros(length(showpj),length(mdname)); mT = merr;
    for pj = 1:length(showpj)
        eval(['load relED_50_' int2str(showpj(pj)) '_case_' caseid{cj} '_' otype '.mat']);
        merr(pj,:) = mean(err);
    end
    merr = roundn(merr,-1);
    fprintf(' distance between the true and the estimated covariance matrix:\n')
    array2table(merr','RowNames',mdname,'VariableNames',string(p(showpj)))
end

otype = 'oc'; % 'oc' for OC outliers
fprintf(['\n outliers type:\t',otype,'\n']);
for cj = showid
    fprintf(['\n outliers are generated from ',casenames{cj},'\n']);
    merr = zeros(length(showpj),length(mdname)); mT = merr;
    for pj = 1:length(showpj)
        eval(['load relED_50_' int2str(showpj(pj)) '_case_' caseid{cj} '_' otype '.mat']);
        merr(pj,:) = mean(err);
    end
    merr = roundn(merr,-1);
    fprintf(' distance between the true and the estimated covariance matrix:\n')
    array2table(merr','RowNames',mdname,'VariableNames',string(p(showpj)))
end
rmpath('./result'); 
