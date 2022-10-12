function [bp,err,T_tr,T_cls,itnum] = mts_1(X,Y,X_label,Y_label,expr)
d = expr.d; pd = prod(d); mdname = expr.mdname; disp_md = 0;
N_tr = size(X,3);
for mj = 1:length(mdname)
    if disp_md fprintf(['\n--> running ' mdname{mj} ':\n']); end
    if sum(strcmpi(mdname{mj},{'tPCA' 'PCA'}))
        x=reshape(X,[pd,size(X,3)]); 
        y=reshape(Y,[pd,size(Y,3)]);
    end
    
    opts=[]; opts.maxit = 1000; opts.disp_it = 0; opts.ini = 'random'; opts.tol = 1e-8; eta = 1e-6;
    switch mdname{mj}
        case 'RFPCA'
            t0 = cputime;
            [bpt, Xc, opts] = mvt_ini(X, opts);
            [bp{mj},opts] = mvt(bpt,Xc,opts);
            %elaps(mj) = opts.time.ini+opts.time.preit+sum(opts.time.it);
            itnum(mj) = opts.itnum;
            T_tr(mj) = cputime-t0;
            t1 = cputime;
            [bp{mj}.lmd{1},bp{mj}.U{1}]=eigdec2(bp{mj}.Sc,d(1),eta);
            [bp{mj}.lmd{2},bp{mj}.U{2}]=eigdec2(bp{mj}.Sr,d(2),eta);
            err{mj} = nnerr(bp{mj},X(:,:,1:length(X_label)),Y,X_label,Y_label);
            T_cls(mj) = cputime-t1;
            
        case 'TPCA'
            t0 = cputime;
            [bpt, Xc, opts] = mxvt_ini(X, opts);
            [bp{mj},opts] = mxvt(bpt,Xc,opts);
            itnum(mj) = opts.itnum;
            T_tr(mj) = cputime-t0;
            t1 = cputime;
            [bp{mj}.lmd{1},bp{mj}.U{1}]=eigdec2(bp{mj}.Sc,d(1),eta);
            [bp{mj}.lmd{2},bp{mj}.U{2}]=eigdec2(bp{mj}.Sr,d(2),eta);
            err{mj} = nnerr(bp{mj},X(:,:,1:length(X_label)),Y,X_label,Y_label);
            T_cls(mj) = cputime-t1;
            
        case 'FPCA'
            t0=cputime;
            [bpt, Xc, opts] = mvn_ini(X, opts);
            [bp{mj},opts] = mvn(bpt,Xc,opts);
            itnum(mj)= opts.itnum;
            T_tr(mj) = cputime-t0;
            t1 = cputime;
            [bp{mj}.lmd{1},bp{mj}.U{1}]=eigdec2(bp{mj}.Sc,d(1),eta);
            [bp{mj}.lmd{2},bp{mj}.U{2}]=eigdec2(bp{mj}.Sr,d(2),eta);
            err{mj} = nnerr(bp{mj},X(:,:,1:length(X_label)),Y,X_label,Y_label);
            T_cls(mj) = cputime-t1;
            
        case 'BPCA'
            t0=cputime;
            bp{mj} = bpca(X);
            T_tr(mj) = cputime-t0;
            itnum(mj) = 0;
            bp{mj}.model = 'BPCA';
            t1 = cputime;
            err{mj} = nnerr(bp{mj},X(:,:,1:length(X_label)),Y,X_label,Y_label);
            T_cls(mj) = cputime-t1;
              
        case 'tPCA'
            t0=cputime;
            [bpt, Xc, opts] = mt_ini(x, opts);
            [bp{mj},opts] = mt(bpt,Xc,opts);
            itnum(mj)=opts.itnum;
            T_tr(mj)=cputime-t0;
            t1 = cputime;
            [lmd,U]=eigdec(bp{mj}.S,pd); qtmp=min(N_tr-10,length(lmd));
            bp{mj}.lmd = lmd(1:qtmp); bp{mj}.U = U(:,1:qtmp);
            err{mj} = nnerr(bp{mj},x(:,1:length(X_label)),y,X_label,Y_label);
            T_cls(mj) = cputime-t1;
            
        case 'PCA'
            t0=cputime;
            opts.center=0;
            bp{mj}.M = mean(x,2); Xc = x-bp{mj}.M;
            bpt = pca(Xc',opts);
            T_tr(mj)=cputime-t0;
            bp{mj}.model = 'PCA';
            t1 = cputime;
            qtmp=min(N_tr-10,length(bpt.lmd));
            bp{mj}.lmd = bpt.lmd(1:qtmp); bp{mj}.U = bpt.U(:,1:qtmp); 
            err{mj} = nnerr(bp{mj},x(:,1:length(X_label)),y,X_label,Y_label);
            T_cls(mj) = cputime-t1;
    end
end
end
