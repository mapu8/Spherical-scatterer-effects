function [p_cmplx, position_vectors] = point_sources_scatterer_pressure(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D,r_scat,indexes)
%POINT_SOURCES_SCATTERER_PRESSURE Compute the pressure field created by a number of 
%                                 point sources with a spherical scatterer in the
%                                 middle of the field.

% ARGUMENTS:
% t - time point
% frequencies - vector of frequencies of the point sources
% source_positions - Nx3 matrix of point source positions in the form of
%                    (azimuth,elevation,distance) in radians and meters
% magnitudes - vector of magnitudes of the point sources
% initial_phases - vector of initial phases of the point sources
% X and Y - meshgrid matrixes
% N - maximum order of the expansion
% D - playback radius in meters
% r_scat - radius of the scatterer in meters
% indexes - for calculating only single points in the grid, OPTIONAL
%
% RETURNS:
% p_cmplx - complex pressure field
% position_vectors - position matrix of the sources in cartesian coordinates

% Number of arguments when 'indexes' is used
args_when_indexes_used = 11;

% Point source expansion
% point source exp(jwt-jk||r_s-r||)/||r_s-r||

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
p_cmplx = zeros(1,size(R,2)); % pressures
position_vectors = zeros(length(frequencies),3); % 3D position vectors in cartesian coordinates

% Calculate the pressure field created by the point sources
for i=1:length(frequencies)
    f = frequencies(i); % frequency
    w = 2*pi*f; % angular frequency
    k = w/c; % wavenumber
    rs = source_positions(i,3); % distance of source
    doa = source_positions(i,1:2); % source direction, azimuth/elevation
    position_vectors(i,:) = rs*[cos(doa(1))*cos(doa(2)) sin(doa(1))*cos(doa(2)) sin(doa(2))]; % source position vector
    m0 = magnitudes(i); % magnitude
    phi0 = initial_phases(i); % initial phase
    a0 = m0.*exp(1i*phi0); % initial amplitude
    
    cos_alpha = sin(doa(2))*sin(elev_stack) + cos(doa(2)).*cos(elev_stack(:)).*cos(doa(1)*ones(size(R(:)))-azi_stack);
    radial_terms = zeros(length(R(:)), N+1);
    angular_terms = zeros(length(R(:)), N+1);
    kr = k*r_stack;
    kr_scat = k*r_scat;
    
    for order=0:N
        jn = sph_besselj(order, kr);
        jnprime = dsph_besselj(order, kr_scat);
        hn = sph_hankel2(order, kr);
        hnprime = dsph_hankel2(order, kr_scat);
        hns = sph_hankel2(order, k*rs);
        radial_terms(:, order+1) = (jn-(jnprime./hnprime).*hn) * hns;
    end
    for degree=0:N
        temp = legendre(degree,cos_alpha);
        angular_terms(:,degree+1) = temp(1,:);
    end
    temporal_terms = exp(1i*w*t);
    cp_stack = (-1i)*k*sum(a0*(ones(size(R(:)))*(2*(0:N)+1)).*radial_terms.*angular_terms.*temporal_terms,2)/(4*pi);
    % Set values to infinity inside the scatterer
    cp_stack(round(r_stack,4)<round(r_scat,4)) = NaN;
    p_temp = reshape(cp_stack, size(R));
    p_cmplx = p_cmplx + p_temp;
    % Set values to infinity outside the reproduction area
    p_cmplx(R>D) = NaN;
end
end

