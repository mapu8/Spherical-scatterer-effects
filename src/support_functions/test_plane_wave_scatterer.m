function [p_cmplx] = test_plane_wave_scatterer(t,f_vec,m0,phi0_vec,dop_vec,X,Y,N,D,r_scat,indexes)
%TEST_PLANE_WAVE_SCATTERER Compute the pressure field with plane wave
%                              expansion with a scatterer in the field.

% ARGUMENTS:
% t - time point
% f - frequency
% m0 - initial magnitude
% phi0 - initial phase
% dop - direction of propagation [azimuth inclination] in radians
% X and Y - meshgrid for X and Y
% N - maximum order of the expansion
% D - playback radius in meters
% r_scat - radius of the scatterer
% indexes - for calculating only single points in the grid, OPTIONAL
%
% RETURNS:
% p_cmplx - complex pressure field


c = 343; % speed of sound
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X))); % spherical coords
% Calculations for all the points in the grid
if nargin < 11
    r_stack = R(:);
    azi_stack = Azi(:);
    elev_stack = Elev(:);
% Calculations for only the selected points
else
    Azi = Azi(indexes);
    Elev = Elev(indexes);
    R = R(indexes);
    r_stack = R(:);
    azi_stack = Azi(:);
    elev_stack = Elev(:);
end

% Plane wave expansion
% plane wave exp(jwt-jkr)
w = 2*pi.*f_vec; % angular frequency
k = w./c; % wavenumber
n_dop = [cos(dop_vec(1,:)').*sin(abs(dop_vec(2,:)')) sin(dop_vec(1,:)').*sin(abs(dop_vec(2,:)')) cos(abs(dop_vec(2,:)'))]; % unit vector of propagation
a0 = m0.*exp(1i*phi0_vec); % initial amplitude

n_doa = -n_dop; % direction of arrival
[doa(:,1), doa(:,2)] = cart2sph(n_doa(:,1),n_doa(:,2),n_doa(:,3));

    size(n_dop)
    size(doa)
    size(sin(doa(:,2)))
    size(sin(elev_stack))
    size(sin(doa(:,2)).*sin(elev_stack))
    size(cos(doa(:,2)).*cos(elev_stack(:)).*cos(doa(:,1).*ones(size(R(:)))))
    size(azi_stack)
    
cos_alpha = sin(doa(:,2)).*sin(elev_stack) + cos(doa(:,2)).*cos(elev_stack(:)).*cos(doa(:,1).*ones(size(R(:)))-azi_stack);
radial_terms = zeros(length(R(:)), N+1, length(f_vec));
kr = k.*r_stack;
kr_scat = k.*r_scat;
for order=0:N
    jn = sph_besselj(order, kr);
    jnprime = dsph_besselj(order, kr_scat);
    hn = sph_hankel2(order, kr);
    hnprime = dsph_hankel2(order, kr_scat);
    radial_terms(:, order+1, :) = 4*pi*1i^order * (jn-(jnprime./hnprime).*hn);
end
angular_terms = zeros(length(R(:)), N+1, length(f_vec));
for degree=0:N
    temp = legendre(degree,cos_alpha);
    angular_terms(:, degree+1 ,:) = temp(1,:);
end
temporal_terms = exp(1i*w*t);

cp_stack = sum(a0.*(ones(size(R(:))).*(2*(0:N)+1)).*radial_terms.*angular_terms.*temporal_terms,2)./(4*pi);
% Set values to infinity inside the scatterer
cp_stack(r_stack<r_scat) = NaN;
p_cmplx = reshape(cp_stack, size(R));
% Set values to infinity outside the reproduction area
p_cmplx(R>D) = NaN;

end

