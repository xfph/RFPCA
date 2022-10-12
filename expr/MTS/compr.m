function [Z,bp] = compr(bp,X)

switch bp.model
    case {'RFPCA', 'TPCA', 'FPCA'}
        [dc,dr,N] = size(X); Z = zeros(dc,dr,N);
        CCi = bp.lmd{1}.^(-.5).*bp.U{1}';
        RRi = bp.U{2}.*(bp.lmd{2}.^(-.5))';
        for i = 1:N
            Z(:, :, i) = CCi*X(:,:,i)*RRi;
        end
    case 'BPCA'
        [dc,dr,N] = size(X); Z = zeros(dc,dr,N);
        for i = 1:N
            Z(:, :, i) = bp.U{1}'*X(:,:,i)*bp.U{2};
        end   
    case {'tPCA', 'PCA'}
        CCi = bp.lmd.^(-.5).*bp.U';
        Z=CCi*X;
end



