function [v_cmplx, position_vectors] = point_sources_scatterer_velocity(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D,r_scat,indexes)
%POINT_SOURCES_SCATTERER_VELOCITY Compute the velocity field created by a number of 
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
% v_cmplx - complex velocity field
% position_vectors - position matrix of the sources in cartesian coordinates

% Point source expansion
% point source exp(jwt-jk||r_s-r||)/||r_s-r||

% Number of arguments when 'indexes' is used
args_when_indexes_used = 11;

c = 343; % speed of sound
rho_0 = 1.2;
Z_0 = c*rho_0; % impedance of air
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
    kr = k*r_stack;
    kr_scat = k*r_scat;
    
    % radial component
    radial_terms = zeros(length(R(:)), N+1);
    for order=0:N
        djn_r = dsph_besselj(order, kr);
        djn_scat = dsph_besselj(order, kr_scat);
        dhn_r = dsph_hankel2(order, kr);
        dhn_scat = dsph_hankel2(order, kr_scat);
        hns = sph_hankel2(order, k*rs);
        radial_terms(:, order+1) = (-1i)*k^2.*(djn_r-(djn_scat./dhn_scat).*dhn_r) * hns;
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
        djn_scat = dsph_besselj(order, kr_scat);
        hn_r = sph_hankel2(order, kr);
        dhn_scat = dsph_hankel2(order, kr_scat);
        hns = sph_hankel2(order, k*rs);
        radial_terms(:, order+1) = (-1i)*k.*(jn_r-(djn_scat./dhn_scat).*hn_r) * hns;
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
        djn_scat = dsph_besselj(order, kr_scat);
        hn_r = sph_hankel2(order, kr);
        dhn_scat = dsph_hankel2(order, kr_scat);
        hns = sph_hankel2(order, k*rs);
        radial_terms(:, order+1) = (-1i)*k.*(jn_r-(djn_scat./dhn_scat).*hn_r) * hns;
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
    temporal_term = exp(1i*w*t);
    cv_x = cv_x + reshape(a0*cv_x_stack*temporal_term, size(R));
    cv_y = cv_y + reshape(a0*cv_y_stack*temporal_term, size(R));
    cv_z = cv_z + reshape(a0*cv_z_stack*temporal_term, size(R));
    % Set values to zero inside the scatterer and outside the reproduction
    % area
    cv_x(round(r_stack,4)<round(r_scat,4)) = 0;
    cv_y(round(r_stack,4)<round(r_scat,4)) = 0;
    cv_z(round(r_stack,4)<round(r_scat,4)) = 0;
    cv_x(R>D) = 0;
    cv_y(R>D) = 0;
    cv_z(R>D) = 0;
end
v_cmplx = cat(3, cv_x, cv_y, cv_z);

end
