% The following codes can be used to reproduce Figure 1 in the main text. 
clc; clear all; close all;
alg={'ECME' 'PX-ECME'}; 
cr={'b:','r-'};
addpath('.\result');

load simu1_low.mat;
%Iteration number
subplot(2,4,1)
for i=1:length(alg)
    g(i)=plot(ll{i},cr{i},'linewidth',1.5);hold on
end
legend(g, alg,'location','best','fontsize',8);
xlabel('Iteration number')
ylabel('Log likelihood')
xlim([-10,max(iternum)+10])
ylim([-6.6e4 -2.85e4]);
%CPU time
subplot(2,4,2)
for i=1:length(alg)
    if strcmp(alg{i},'ECME') T{i}=T{i}-(T{i}(1)-T{i+1}(1)); end
    g(i)=plot(T{i},ll{i},cr{i},'linewidth',1.5);hold on
end
legend(g, alg,'location','best','fontsize',8);
xlabel('CPU time');
xlim([-0.02,max(T1)+0.02])
ylim([-6.6e4 -2.85e4]);

load simu1_high.mat;
%Iteration number
subplot(2,4,3)
for i=1:length(alg)
    g(i)=plot(ll{i},cr{i},'linewidth',1.5);hold on
end
legend(g, alg,'location','best','fontsize',8);
xlabel('Iteration number')
ylabel('Log-likelihood')
xlim([-500,max(iternum)+500])
ylim([-4.815e6 -4.8035e6]);

%CPU time
subplot(2,4,4)
for i=1:length(alg)
    if strcmp(alg{i},'ECME') T{i}=T{i}-(T{i}(1)-T{i+1}(1)); end
    g(i)=plot(T{i}(1:iternum(i)),ll{i},cr{i},'linewidth',1.5);hold on
end
legend(g, alg,'location','best','fontsize',8);
xlabel('CPU time');
xlim([-100,max(T1)+100])
ylim([-4.815e6 -4.8035e6]);

rmpath('.\result');