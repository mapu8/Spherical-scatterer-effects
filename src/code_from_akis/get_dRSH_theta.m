function dRSH_theta = get_dRSH_theta(N, dirs)
%DSH_PHI Summary of this function goes here
%   Detailed explanation goes here

Y_N_plus1 = getSH(N+1,dirs,'real').';
nSH = (N+1)^2;
dRSH_theta = zeros(nSH, size(Y_N_plus1,2));

cot_theta = cot(dirs(:,2));
csc_theta = csc(dirs(:,2));

for n=0:N
    m = (-n:n).';
    idx_n = n^2+n+1+m;
    idx_n_plus1 = (n+1)^2+n+2+m;
    Y_n = Y_N_plus1(idx_n,:);
    Y_n_plus1 = Y_N_plus1(idx_n_plus1,:);
    
    A = -(n+1)*ones(2*n+1,1)*cot_theta.';
    B = sqrt( (n-abs(m)+1).*(n+abs(m)+1)*(2*n+1)/(2*n+3) )*csc_theta.';

    dRSH_theta(idx_n,:) = A.*Y_n + B.*Y_n_plus1;
end
dRSH_theta = dRSH_theta.';

end
