% The following codes can be used to reproduce Figure 2 in the main text. 
clc; clear all;
model = {'RFPCA','TPCA','FPCA','tPCA'};
N_tr = [60 100 200 500 1000]; nj = length(N_tr);
nu = [3,inf]; 

addpath('./result');
% plot test log-likelihood
cr={'r-p','b:*','g--o','c-.>','k:'};
for j = 1:length(nu)
    eval(['load simu2_test_llh_',int2str(j) '.mat']);
    subplot(2,2,j);
    for i = 1:length(model)+1
        switch i
            case {1,2,3,4} 
                g(i) = errorbar(1:nj,mean(llh_est{i}),std(llh_est{i}),cr{i},'linewidth',1,'MarkerSize',3);hold on %
            case 5
                g(i) = plot(repelem(llh_true,nj),cr{i},'linewidth',1);hold on
        end
    end
    legend({'RFPCA','TPCA','FPCA','{\itt}PCA','true'},'location', 'best','fontsize',8)
    xlim([-1/50,nj+1/50]);
    if nu(j)==3 ylim([-4.21e5,-2.948e5]);
    else ylim([-4.82e5,-2.448e5]);
    end
    set(gca, 'xtick', 1:nj);
    set(gca,'XTicklabel', {'60','100','200','500','1000'});
    xlabel('Sample size')
    if j==1 ylabel('Test log-likelihood'); end
end
rmpath('./result');
    
