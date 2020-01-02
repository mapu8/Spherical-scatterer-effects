function dRSH_phi = get_dRSH_phi(N, dirs)
%DSH_PHI Summary of this function goes here
%   Detailed explanation goes here

Y_N = getSH(N,dirs,'real').';
dRSH_phi = zeros(size(Y_N));
for n=0:N
    idx = n^2+(1:2*n+1);
    Y_n = Y_N(idx,:);
    M = -diag(-n:n);
    dRSH_phi(idx,:) = M * Y_n(end:-1:1,:);
end
dRSH_phi = dRSH_phi.';

end
