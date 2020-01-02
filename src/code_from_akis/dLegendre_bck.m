function dP_n = dLegendre_bck(n,x)
%DLEGENDRE_BCK Summary of this function goes here
%   Detailed explanation goes here

if iscolumn(x), x = x.'; end

P_n_all = legendre(n,x);
P_n_minus1_all = legendre(n-1,x);

A = ones(n,1) * (1./(x.^2-1));
B = ones(n,1) * x;
C = diag(n+(0:n-1));

dP_n = A .* (n*B.*P_n_all(1:end-1,:) - C*P_n_minus1_all(1:end,:));

end
