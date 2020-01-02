function dP_n = dLegendre_fwd(n,x)
%DLEGENDRE_FWD Summary of this function goes here
%   Detailed explanation goes here

if iscolumn(x), x = x.'; end

P_n_all = legendre(n,x);
P_n_plus1_all = legendre(n+1,x);

A = ones(n+1,1) * (1./(1-x.^2));
B = ones(n+1,1) * x;
C = diag(n-(0:n)+1);

dP_n = A .* ((n+1)*B.*P_n_all - C*P_n_plus1_all(1:end-1,:));

end
