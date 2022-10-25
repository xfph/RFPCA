% This code to show the lowest error rates of different methods on the AUSLAN dataset by one random splitting.
%
% To run the results on the ECG dataset, you can download ECG data from http://www.mustafabaydogan.com./files/viewcategory/20-data-sets.html, 
% and replace the AUSLAN data set in this M-file, by dividing it into training and test sets.
clc;clear all;
addpath('../../prog'); addpath('../../prog/common');
dataname = 'AUS';  dis_dataname = 'AUSLAN';
casename = {'without','caseIV'}; 
dis_casename = {'Without outliers','Case IV: p = 20%%, U(0,10)'};
mdname={'RFPCA','FPCA','BPCA','PCA'}; % TPCA,'tPCA'

fprintf(['Data set: ' dis_dataname '\n']);
for cj=1:length(casename)
    fprintf(['\n----> ' dis_casename{cj} ':\n']);
    load([dataname,'_',casename{cj}]);
    [dc,dr,N] = size(X);
    expr=struct(); expr.d = [dc,dr]; expr.mdname=mdname;
    [bp,err,T_tr,T_cls,~] = mts_1(X,Y,X_label,Y_label,expr);
    
    fprintf('The lowest error rates (mean ± std) and their corresponding dimensions by different methods:\n')
    I=zeros(length(mdname),2);
    for i=1:length(mdname)
        switch mdname{i}
            case {'RFPCA','TPCA','FPCA','BPCA'}
                [merr,i1]=min(err{i});
                [opte,i2]=min(merr);i1=i1(i2);
                mme(i)=opte; I(i, :)=[i1 i2];
                fprintf('%s\t%3.1f(%2d,%2d)\n',mdname{i},mme(i)*100,I(i,1),I(i,2));
                
            case {'tPCA','PCA'}
                [opte,i2]=min(err{i});
                mme(i)=opte; I(i, 1)=i2;
                fprintf('%s\t\t%3.1f(%2d)\n',mdname{i},mme(i)*100,I(i,1));
        end
    end
end
rmpath('../../prog'); rmpath('../../prog/common');

