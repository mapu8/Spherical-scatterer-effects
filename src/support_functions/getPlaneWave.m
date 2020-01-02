function [p_cmplx] = getPlaneWave(t,f,m0,phi0,dop,X,Y,D)
%GETPLANEWAVE Get the complex pressure field of a plane wave.
%
% ARGUMENTS:
% t - time point
% f - frequency
% m0 - initial magnitude
% phi0 - initial phase
% dop - direction of propagation [azimuth inclination]
% X, Y and Z - meshgrid for X, Y and Z
% D - playback radius in meters
%
% RETURNS:
% p_cmplx - complex pressure field

% Plane wave exp(jwt-jkr)
c = 343;
rho_0 = 1.2;
Z_0 = c*rho_0;
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
a0 = m0.*exp(1i*phi0); % initial amplitude
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X))); % spherical coords

%n_dop = [cos(dop(1))*cos(dop(2)) sin(dop(1))*cos(dop(2)) sin(dop(2))] % unit vector of propagation
n_dop = [cos(dop(1))*sin(abs(dop(2))) sin(dop(1))*sin(abs(dop(2))) cos(abs(dop(2)))]; % unit vector of propagation
p_cmplx = a0*exp(1i* (w*t*ones(size(X)) - k* (n_dop(1)*X + n_dop(2)*Y))); % complex pressure field
p_cmplx(R>D) = NaN;
end

