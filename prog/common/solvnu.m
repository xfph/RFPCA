function root=solvnu(nu, d2, dim)
ndata=length(d2);

%Caculate the value of omega 
omega=(nu+dim)./(d2+nu);
root=diag_ln((nu+dim)/2)-diag_ln(nu/2)+1/ndata*sum(log(omega)-omega)+1;


%Caculate the difference between a digamma and a log function with the same
%parameter
function s=diag_ln(x)
s=digamma(x)-log(x);


