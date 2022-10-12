clc;clear all;

addpath('../../prog'); addpath('../../prog/common');
dataname = 'AUS';
dis_casename = {'without','caseIV'}; % without: without outliers; caseIV: p=20%,U(0,10)
mdname={'RFPCA','FPCA','BPCA','PCA'}; % TPCA,'tPCA'

for cj=1:length(dis_casename)
    fprintf(['\n---->' dataname,'\t' dis_casename{cj} ':\n']);
    load([dataname,'_',dis_casename{cj}]);
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

