function [v_cmplx] = plane_wave_expansion_velocity(t,f,m0,phi0,doa,X,Y,N,D,indexes)
%PLANE_WAVE_EXPANSION_VELOCITY Compute the particle velocity field created
%                              by a number of plane waves.

% ARGUMENTS:
% t - time point
% f - frequency
% m0 - initial magnitude
% phi0 - initial phase
% doa - direction of arrival [azimuth elevation] in radians
% X and Y - meshgrid for X and Y
% N - maximum order of the expansion
% D - playback radius in meters
% indexes - for calculating only single points in the grid, OPTIONAL
%
% RETURNS:
% v_cmplx - complex particle velocity field

% Number of arguments when 'indexes' is used
args_when_indexes_used = 10;

c = 343; % speed of sound
rho_0 = 1.2;
Z_0 = c*rho_0;
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X))); % spherical coords
% Calculations for all the points in the grid
if nargin < args_when_indexes_used
    r_stack = R(:);
    azi_stack = Azi(:);
    elev_stack = Elev(:);
    cv_x = zeros(size(R,1));
    cv_y = zeros(size(R,1));
    cv_z = zeros(size(R,1));
% Calculations for only the selected points
else
    Azi = Azi(indexes);
    Elev = Elev(indexes);
    R = R(indexes);
    r_stack = R(:);
    azi_stack = Azi(:);
    elev_stack = Elev(:);
    cv_x = zeros(size(R,1),1);
    cv_y = zeros(size(R,1),1);
    cv_z = zeros(size(R,1),1);
end

% Plane wave expansion
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
%n_dop = [cos(dop(1))*sin(abs(dop(2))) sin(dop(1))*sin(abs(dop(2))) cos(abs(dop(2)))]; % unit vector of propagation
a0 = m0.*exp(1i*phi0); % initial amplitude

%n_doa = -n_dop; % direction of arrival
%[doa(1), doa(2)] = cart2sph(n_doa(1),n_doa(2),n_doa(3));

cos_alpha = sin(doa(2))*sin(elev_stack) + cos(doa(2)).*cos(elev_stack(:)).*cos(doa(1)*ones(size(R(:)))-azi_stack);
radial_terms = zeros(length(R(:)), N+1);
kr = k*r_stack;

% radial component
for order=0:N
    djn_r = dSphBesselj_fwd(order, kr);
    radial_terms(:, order+1) = 4*pi*1i^order*k*djn_r;
end
angular_terms = zeros(length(R(:)), N+1);
for degree=0:N
    temp = legendre(degree,cos_alpha);
    angular_terms(:,degree+1) = temp(1,:);
end
cv_r_stack = -1/(1i*k*Z_0)*sum((ones(size(R(:)))*(2*(0:N)+1)).*radial_terms.*angular_terms,2)/(4*pi);

% azimuthal component
radial_terms = zeros(length(R(:)), N+1);
for order=0:N
    jn_r = sph_besselj(order, kr);
    radial_terms(:, order+1) = 4*pi*1i^order*jn_r;
end
radial_terms = replicatePerOrder(radial_terms,2);
nDirs = length(azi_stack);
angular_terms = (ones(nDirs,1)*getSH(N,[doa(1) pi/2-doa(2)],'real')) .* ...
                get_dRSH_phi(N,[azi_stack pi/2-elev_stack]);
cv_phi_stack = -1/(1i*k*Z_0)* (1./(r_stack.*sin(pi/2-elev_stack))) .* sum(radial_terms.*angular_terms,2);

% polar component
radial_terms = zeros(length(R(:)), N+1);
for order=0:N
    jn_r = sph_besselj(order, kr);
    radial_terms(:, order+1) = 4*pi*1i^order*jn_r;
end
radial_terms = replicatePerOrder(radial_terms,2);
nDirs = length(azi_stack);
angular_terms = (ones(nDirs,1)*getSH(N,[doa(1) pi/2-doa(2)],'real')) .* ...
                get_dRSH_theta(N,[azi_stack pi/2-elev_stack]);
cv_theta_stack = -1/(1i*k*Z_0)* (1./r_stack) .* sum(radial_terms.*angular_terms,2);

% convert to cartesian coordinates
cv_x_stack = sin(pi/2-elev_stack).*cos(azi_stack).*cv_r_stack + ...
             cos(pi/2-elev_stack).*cos(azi_stack).*cv_theta_stack - ...
             sin(azi_stack).*cv_phi_stack;
cv_y_stack = sin(pi/2-elev_stack).*sin(azi_stack).*cv_r_stack + ...
             cos(pi/2-elev_stack).*sin(azi_stack).*cv_theta_stack + ...
             cos(azi_stack).*cv_phi_stack;
cv_z_stack = cos(pi/2-elev_stack).*cv_r_stack - ...
             sin(pi/2-elev_stack).*cv_theta_stack;

% temporal term
temporal_term = exp(-1i*w*t);
cv_x = cv_x + reshape(cv_x_stack*temporal_term, size(R));
cv_y = cv_y + reshape(cv_y_stack*temporal_term, size(R));
cv_z = cv_z + reshape(cv_z_stack*temporal_term, size(R));        

% Set values to zero outside the reproduction area
cv_x(R>D) = 0;
cv_y(R>D) = 0;
cv_z(R>D) = 0;

v_cmplx = cat(3, cv_x, cv_y, cv_z);

end

