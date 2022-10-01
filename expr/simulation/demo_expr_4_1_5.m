% The following codes can be used to reproduce Section 4.1.5 in the main
% text. 
clear all;clc;clf;
mdname = {'RFPCA','TPCA','FPCA'}; 
N_tr = [200,500,1000,2000,4000,8000,13000];
rep = 1; 
cr={'r-p','b:*','g--o'}; 

addpath('.\result');
load('simu5_t_itnum_oc.mat')
for mj = 1:length(mdname)
    g(mj) = plot(T(:,mj),cr{mj},'linewidth',0.5,'MarkerSize',6);
    hold on
end
set(gca, 'xtick', [1:length(N_tr)]);
set(gca,'XTicklabel',string(N_tr));
xlim([0.5 length(N_tr)+0.5]); ylim([-400 7500]);  
h=legend(g,mdname,'location', 'best','fontsize',8);
xlabel('Sample size','fontsize',10);
ylabel('CPU time','fontsize',10);

dis_Ni = [2,4,6,7];
fprintf('\n CPU time \n')
array2table(T(dis_Ni,1:3),'RowNames',{'N=500','N=2000','N=8000','N=13000'},'VariableNames',{'RFPCA','TPCA','FPCA'}')
fprintf('\n Iteration number \n')
array2table(itnum(dis_Ni,1:3),'RowNames',{'N=500','N=2000','N=8000','N=13000'},'VariableNames',{'RFPCA','TPCA','FPCA'}')





