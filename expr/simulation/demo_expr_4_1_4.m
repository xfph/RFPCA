% The following codes can be used to reproduce Figure 3 in the main text.
% hat_tau <c ,as an outlier.
clc; clear all; clf;
model = {'RFPCA','{\it t}PCA'};
casename = {'$U(100,110)$','$U(100,102)$','$U(100000,100002)$'};
alpha=0.05;
cr={'r','g','k','c','b','m'};  sb={'o','d','.','^'};

addpath('.\result');
for cj = 1:3
    eval(['load simu4_caseoc' int2str(cj) '.mat;']);
    N=length(bp{1}.tau); N1=1000; N2=N-N1;
    x=zeros(2,N);
    x(1,:)=bp{1}.tau; x(2,:)=bp{2}.tau;
    % critical value ----(Multivariate t nonlinear mixed....equation(24))
    k1=(1+bp{1}.d(1)*bp{1}.d(2)/bp{1}.nu)*betainv(alpha,bp{1}.nu/2,bp{1}.d(1)*bp{1}.d(2)/2);
    k2=(1+bp{2}.d/bp{2}.nu)*betainv(alpha,bp{2}.nu/2,bp{2}.d/2);
    
    subplot(2,3,cj)
    p1(1)=plot((1:N1/5),x(1,1:N1/5),[cr{1} sb{1}],'MarkerSize',3); hold on;
    % Pick something greater(less) than k2
    tmp = (N1/5+1):(N1/5+3*N2/5);
    if cj==1
        p1(2)=plot(tmp(1:3:end),x(1,(N1+1:N-4*N2/5)),[cr{2} sb{2}],'MarkerSize',3);hold on;
    elseif cj==2
        y1=[1042 1024 1046 1021];
        random_num=[1008 1030 1012 1031 1047 1033];
        y2=sort([random_num y1]);
        p1(2)=plot(tmp(1:3:end),x(1,y2),[cr{2} sb{2}],'MarkerSize',3);
        hold on;
    else cj==3
        y1=[1004 1011 1031 1030 1015];
        random_num =[1006 1047 1001 1038 1039];
        y2=sort([random_num y1]);
        p1(2)=plot(tmp(1:3:end),x(1,y2),[cr{2} sb{2}],'MarkerSize',3);
        hold on;
    end
    % critical horizontal line  (if  \tau<k, an outlier)
    xmin=ceil(0.01*N); xmax=ceil(0.225*N);
    p1(3)=plot((1:xmax-8),repmat(k1,1,xmax-8),'m');
    
    if cj==1 h=legend(p1, '$\mathbf{E}(\tau_n)$ of non-outliers', '$\mathbf{E}(\tau_n)$ of OC outliers','critical value','location', 'best');
        set(h,'Interpreter','latex','fontsize',6); end
    xlim([-xmin xmax]);
    ylim([min(x(1,:))-0.1,max(x(1,:))+0.1]);
    if cj==1 ylabel({'weight of RFPCA'},'fontsize',10); end
    title(casename{cj},'Interpreter','latex','fontsize',10);
    
    subplot(2,3,cj+3)
    p2(1)=plot((1:N1/5),x(2,1:N1/5),[cr{3} sb{3}],'MarkerSize',3);
    hold on;
    % Pick something greater(less) than k2
    tmp = (N1/5+1):(N1/5+3*N2/5);
    if cj==1
        p2(2)=plot(tmp(1:3:end),x(1,(N1+1:N-4*N2/5)),[cr{4} sb{4}],'MarkerSize',3);hold on;
    elseif cj==2
        y1=[1042 1024 1046 1021];
        random_num=[1008 1030 1012 1031 1047 1033];
        y2=sort([random_num y1]);
        p2(2)=plot(tmp(1:3:end),x(2,y2),[cr{4} sb{4}],'MarkerSize',3);
        hold on;
    else cj==3
        y1=[1004 1011 1031 1030 1015];
        random_num =[1006 1047 1001 1038 1039];
        y2=sort([random_num y1]);
        p2(2)=plot(tmp(1:3:end),x(2,y2),[cr{4} sb{4}],'MarkerSize',3);
        hold on;
    end
    xmin=ceil(0.01*N); xmax=ceil(0.225*N);
    p2(3)=plot((1:xmax-7),repmat(k2,1,xmax-7),'m');
    if cj==1
        ylim([-0.05,max(x(2,1:N1/5))+0.5]);
        set(gca,'ytick',0:0.5:max(x(2,1:N1/5))+0.5);
    else
        ylim([-0.55,max(x(2,1:N1/5))]+0.5);
        set(gca,'ytick',-0.5:0.5:max(x(2,1:N1/5))+0.5);
    end
    if cj==1 h=legend(p2, '$\mathbf{E}(\tau_n)$ of non-outliers', '$\mathbf{E}(\tau_n)$ of OC outliers','critical value','location', 'best');
        set(h,'Interpreter','latex','fontsize',6); end  %legend boxoff;
    xlim([-xmin,xmax]);
    xlabel('Sample number','fontsize',10);
    if cj==1 ylabel({'weight of {\itt}PCA'},'fontsize',10); end
end
rmpath('.\result')

