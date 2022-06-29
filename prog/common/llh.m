function L = llh(bp,X,md)
%  bp: bp.M,bp.Sc,bp.Sr,bp.nu
%  md: 'MVT','MVN','MN','mt'   
%  MVT: matrix variate t; MVN: matrix variate normal; T: matrix variate T 
%  MN: Multivariate normal; mt: Multivariate t
eta = 1e-6;
switch md
    case {'MVT','MVN','T'}
        [d(1),d(2),N] = size(X);  pd =prod(d);  X = X-bp.M;
        [lc,Uc]=eigdec2(bp.Sc,d(1),eta); Ci=Uc.*(lc.^(-.5))';
        [lr,Ur]=eigdec2(bp.Sr,d(2),eta); Ri=Ur.*(lr.^(-.5))';
        logdet_Sc = sum(log(lc));  logdet_Sr = sum(log(lr));
        switch md
            case {'MVT','MVN'}
                CX = reshape(permute(reshape(Ci'*X(:,:),[d(1) d(2) N]),[1 3 2]),[d(1)*N d(2)]);
                CXR = reshape(permute(reshape(CX*Ri,[d(1) N d(2)]),[1 3 2]),[pd N]);
                tr = sum(CXR.^2,1);
            case 'T'
                Slogdet_juli = 0; 
                for i = 1:N
                    Sci = Ci*Ci'; XR = X(:,:,i)*Ri; 
                    XSrXT = XR*XR'; ScXSrXT = Sci*XSrXT;
                    Slogdet_juli = Slogdet_juli+logdet2(eye(d(1))+ScXSrXT);
                end
        end
        switch md
            case 'MVN'
                L = -N/2*(d(1)*d(2)*log(2*pi)+d(2)*logdet_Sc+d(1)*logdet_Sr)-1/2*sum(tr);
            case 'MVT'
                ll1= -N/2 * (d(2) * logdet_Sc + d(1) * logdet_Sr + pd * log(pi*bp.nu));
                ll2 = N * ( gammaln((pd+bp.nu)/2) - gammaln((bp.nu/2)) );
                tr0 = sum(log(tr./bp.nu+1))*((pd+bp.nu)/2);
                L = ll1 + ll2 - tr0;
            case 'T'
                ll1= -N/2 * (d(2) * logdet_Sc + d(1) * logdet_Sr + pd * log(pi));
                ll2 = N * ( mvgammaln((bp.nu+sum(d)-1)/2,d(1)) - mvgammaln((bp.nu+d(1)-1)/2,d(1)) );
                ll3 = Slogdet_juli*(bp.nu+sum(d)-1)/2;
                L = ll1 + ll2 - ll3;
        end
    case {'MN','mt'}
        [d,N] = size(X); X = X-bp.M(:);
        [lmd,U] = eigdec2(bp.S,d,eta); logdet_S = sum(log(lmd));
        C = U.*(lmd.^(-.5))';  XTC = X'*C;  
        switch md
            case 'MN'
                tr = sum(sum(XTC.^2));
                L = -N/2*(d*log(2*pi)+logdet_S)-1/2*tr;
            case 'mt'
                tr = sum(XTC.^2,2);
                ll1 = -N/2*(logdet_S+d*log(pi*bp.nu))+N*(gammaln((d+bp.nu)/2)-gammaln(bp.nu/2));
                ll2 = sum(log(tr./bp.nu+1))*((d+bp.nu)/2);
                L = ll1 - ll2;
        end 
end

    

