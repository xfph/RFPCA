function [vals,vec] = eigdec2(S,d,eta)
[l, U] = eigdec(S,d); 
l(l<eta)=eta;
vals=l;
vec=U;