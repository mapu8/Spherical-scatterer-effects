function [p_cmplx] = plane_wave_scatterer_pressure(t,f,m0,phi0,doa,X,Y,N,D,r_scat,indexes)
%PLANE_WAVE_SCATTERER_PRESSURE Compute the pressure field with plane wave
%                              expansion with a scatterer in the field.

% ARGUMENTS:
% t - time point
% f - frequency
% m0 - initial magnitude
% phi0 - initial phase
% doa - direction of arrival [azimuth elevation] in radians
% X and Y - meshgrid for X and Y
% N - maximum order of the expansion
% D - playback radius in meters
% r_scat - radius of the scatterer
% indexes - for calculating only single points in the grid, OPTIONAL
%
% RETURNS:
% p_cmplx - complex pressure field

% Number of arguments when 'indexes' is used
args_when_indexes_used = 11;

c = 343; % speed of sound
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X))); % spherical coords
% Calculations for all the points in the grid
if nargin < args_when_indexes_used
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
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
%n_dop = [cos(dop(1))*sin(abs(dop(2))) sin(dop(1))*sin(abs(dop(2))) cos(abs(dop(2)))]; % unit vector of propagation
a0 = m0.*exp(1i*phi0); % initial amplitude

%n_doa = -n_dop; % direction of arrival
%[doa(1), doa(2)] = cart2sph(n_doa(1),n_doa(2),n_doa(3));

cos_alpha = sin(doa(2))*sin(elev_stack) + cos(doa(2)).*cos(elev_stack(:)).*cos(doa(1)*ones(size(R(:)))-azi_stack);
radial_terms = zeros(length(R(:)), N+1);
kr = k*r_stack;
kr_scat = k*r_scat;
for order=0:N
    jn = sph_besselj(order, kr);
    jnprime = dsph_besselj(order, kr_scat);
    hn = sph_hankel2(order, kr);
    hnprime = dsph_hankel2(order, kr_scat);
    radial_terms(:, order+1) = 4*pi*1i^order * (jn-(jnprime./hnprime).*hn);
end
angular_terms = zeros(length(R(:)), N+1);
for degree=0:N
    temp = legendre(degree,cos_alpha);
    angular_terms(:,degree+1) = temp(1,:);
end
temporal_terms = exp(1i*w*t);

cp_stack = sum(a0*(ones(size(R(:)))*(2*(0:N)+1)).*radial_terms.*angular_terms*temporal_terms,2)/(4*pi);
% Set values to infinity inside the scatterer
cp_stack(round(r_stack,4)<round(r_scat,4)) = NaN;
p_cmplx = reshape(cp_stack, size(R));
% Set values to infinity outside the reproduction area
p_cmplx(R>D) = NaN;

end

