function err = nnerr(bp,X,Y,X_label,Y_label)
model=bp.model;  
switch model
    case {'RFPCA', 'TPCA', 'FPCA', 'BPCA'}
        N_test = size(Y,3); N_train = length(X_label);
        Z_tr = compr(bp,X); Z_te = compr(bp,Y);
        err = zeros(bp.d(1),bp.d(2));
        Dtmp = zeros(bp.d(1), bp.d(2),N_train);
        I = zeros(bp.d(1), bp.d(2),N_test);
        for i=1:N_test
            D=Dtmp;
            for j=1:N_train
                Z=Z_te(:, :, i)-Z_tr(:, :, j);
                D(:, :, j)=cumsum(cumsum(Z.^2, 1), 2);
            end
            [~,index]=min(D,[],3);
            I(:, :, i)=index;
        end
        for i=1:bp.d(1)
            for j=1:bp.d(2)
                index=squeeze(I(i, j, :));
                err(i,j)=sum(X_label(index)~=Y_label)/N_test;
            end
        end
        
    case {'tPCA' 'PCA'}
        N_test = size(Y,2); N_train = length(X_label);
        Z_tr = compr(bp,X); Z_te = compr(bp,Y);
        d = size(Z_tr,1);
        err = zeros(d,1);  Dtmp = zeros(d,N_train); I = zeros(d,N_test);
        for i=1:N_test
            D=Dtmp;
            for j=1:N_train
                Z=Z_te(:,i)-Z_tr(:,j);
                D(:,j)=cumsum(Z.^2, 1);
            end
            [~,index]=min(D,[],2);
            I(:,i)=index;
        end
        for i=1:d
            index=I(i,:);
            err(i)=sum(X_label(index)~=Y_label)/N_test;
        end
            
end

