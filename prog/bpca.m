function bp = bpca(X)
[d(1),d(2),N] = size(X); bp.d = d;
M=mean(X,3); X=X-M;
X1 = X(:, :); %[X1,X2...XN]  [c r*N];
X2 = permute(X, [2, 1, 3]); X2 = X2(:, :); %[r c*N];

St=cell(2,1); bp.lmd=St; bp.U=St; 
St{1}=zeros(d(1));St{2}=zeros(d(2));
St{1}=X1*X1'/N; St{2}=X2*X2'/N;
for l=1:2
    [evals, vec] = eigdec(St{l}, d(l));
    bp.lmd{l} = evals;
    bp.U{l}=vec;
end
bp.Sc = St{1}; bp.Sr = St{2}; 