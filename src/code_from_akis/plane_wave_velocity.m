
% acoustics
c = 343;
rho_0 = 1.2;
Z_0 = c*rho_0;

% playback diameter for computations, in meters
D = 6;
% grid resolution, in meters
d = 0.2;
x = -D:d:D;
y = -D:d:D;
[X,Y] = meshgrid(x,y);
Z = zeros(size(X));
% spherical coords
[Azi,Elev,R] = cart2sph(X,Y,Z);

x_stack = X(:);
y_stack = Y(:);
r_stack = R(:);
azi_stack = Azi(:);
elev_stack = Elev(:);

%% signal params

t = 0; % time point
f = 100; % frequency
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
m0 = 1; % magnitude
phi0 = 0; % initial phase
a0 = m0.*exp(1i*phi0); % initial amplitude

dop = [pi/2 0]; % direction of propagation, azimuth/elevation
n_dop = [cos(dop(1))*cos(dop(2)) sin(dop(1))*cos(dop(2)) sin(dop(2))]; % unit vector of propagation
n_doa = -n_dop; % direction of arrival
[doa(1), doa(2)] = cart2sph(n_doa(1),n_doa(2),n_doa(3));

kr = k*r_stack;

%% direct modeling - plane wave exp(jwt-jkr)

p_cmplx = a0*exp(1i* (w*t*ones(size(X)) - k* (n_dop(1)*X + n_dop(2)*Y))); % complex pressure field
v_cmplx = (1/Z_0)*repmat(p_cmplx,[1 1 3]).*cat(3, ones(size(X))*n_dop(1), ones(size(X))*n_dop(2), ones(size(X))*n_dop(3));
v_real = real(v_cmplx);
v_mag = sqrt(sum(v_real.^2,3));

% draw
v_mag(R>D/2) = Inf;
figure
subplot(1,2,1)
surf(X,Y,v_mag); colorbar
shading interp;
arrow1 = -1.2*D/2*n_dop; arrow2 = -D/2*n_dop; darrow = arrow2-arrow1;
hold on
quiver(arrow1(1),arrow1(2),darrow(1),darrow(2),0, 'linewidth',5,'color','k','maxheadsize',1);
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

subplot(1,2,2)
vx_real = v_real(:,:,1); vx_real(R>D/2) = 0;
vy_real = v_real(:,:,2); vy_real(R>D/2) = 0;
vz_real = v_real(:,:,3); vz_real(R>D/2) = 0;
quiver(X,Y, vx_real, vy_real);
arrow1 = -1.2*D/2*n_dop; arrow2 = -D/2*n_dop; darrow = arrow2-arrow1;
hold on
quiver(arrow1(1),arrow1(2),darrow(1),darrow(2),0, 'linewidth',5,'color','k','maxheadsize',1);
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

%% plane wave expansion

cos_alpha = sin(doa(2))*sin(elev_stack) + cos(doa(2)).*cos(elev_stack(:)).*cos(doa(1)*ones(size(X(:)))-azi_stack);
order_max = ceil(exp(1)*k*(D/2)/2);
radial_terms = zeros(length(X(:)), order_max+1);
% radial component
for order=0:order_max
    djn_r = dSphBesselj_fwd(order, kr);
    radial_terms(:, order+1) = 4*pi*1i^order*k*djn_r;
end
angular_terms = zeros(length(X(:)), order_max+1);
for degree=0:order_max
    temp = legendre(degree,cos_alpha);
    angular_terms(:,degree+1) = temp(1,:);
end
cv_r_stack = -1/(1i*k*Z_0)*sum((ones(size(X(:)))*(2*(0:order_max)+1)).*radial_terms.*angular_terms,2)/(4*pi);

% azimuthal component
radial_terms = zeros(length(X(:)), order_max+1);
for order=0:order_max
    jn_r = sph_besselj(order, kr);
    radial_terms(:, order+1) = 4*pi*1i^order*jn_r;
end
radial_terms = replicatePerOrder(radial_terms,2);
nDirs = length(azi_stack);
angular_terms = (ones(nDirs,1)*getSH(order_max,[doa(1) pi/2-doa(2)],'real')) .* ...
                get_dRSH_phi(order_max,[azi_stack pi/2-elev_stack]);
cv_phi_stack = -1/(1i*k*Z_0)* (1./(r_stack.*sin(pi/2-elev_stack))) .* sum(radial_terms.*angular_terms,2);

% polar component
radial_terms = zeros(length(X(:)), order_max+1);
for order=0:order_max
    jn_r = sph_besselj(order, kr);
    radial_terms(:, order+1) = 4*pi*1i^order*jn_r;
end
radial_terms = replicatePerOrder(radial_terms,2);
nDirs = length(azi_stack);
angular_terms = (ones(nDirs,1)*getSH(order_max,[doa(1) pi/2-doa(2)],'real')) .* ...
                get_dRSH_theta(order_max,[azi_stack pi/2-elev_stack]);
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

cv_x = reshape(cv_x_stack, size(X));
cv_y = reshape(cv_y_stack, size(X));
cv_z = reshape(cv_z_stack, size(X));        
         
v_cmplx = cat(3, cv_x, cv_y, cv_z);
v_real = real(v_cmplx);
v_mag = sqrt(sum(v_real.^2,3));

% draw
v_mag(R>D/2) = Inf;
figure
subplot(1,2,1)
surf(X,Y,v_mag); colorbar
shading interp;
arrow1 = -1.2*D/2*n_dop; arrow2 = -D/2*n_dop; darrow = arrow2-arrow1;
hold on
quiver(arrow1(1),arrow1(2),darrow(1),darrow(2),0, 'linewidth',5,'color','k','maxheadsize',1);
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

subplot(1,2,2)
vx_real = v_real(:,:,1); vx_real(R>D/2) = 0;
vy_real = v_real(:,:,2); vy_real(R>D/2) = 0;
vz_real = v_real(:,:,3); vz_real(R>D/2) = 0;
quiver(X,Y, vx_real, vy_real);
arrow1 = -1.2*D/2*n_dop; arrow2 = -D/2*n_dop; darrow = arrow2-arrow1;
hold on
quiver(arrow1(1),arrow1(2),darrow(1),darrow(2),0, 'linewidth',5,'color','k','maxheadsize',1);
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

%% spherical scatterer - plane wave

% scatterer radius
r_scat = 1;
kr_scat = k*r_scat;

radial_terms = zeros(length(X(:)), order_max+1);
% radial component
for order=0:order_max
    djn_r = dsph_besselj(order, kr);
    djn_scat = dsph_besselj(order, kr_scat);
    dhn_r = dsph_hankel2(order, kr);
    dhn_scat = dsph_hankel2(order, kr_scat);
    radial_terms(:, order+1) = 4*pi*1i^order*k* (djn_r-(djn_scat./dhn_scat).*dhn_r);
end
angular_terms = zeros(length(X(:)), order_max+1);
for degree=0:order_max
    temp = legendre(degree,cos_alpha);
    angular_terms(:,degree+1) = temp(1,:);
end
cv_r_stack = -1/(1i*k*Z_0)*sum((ones(size(X(:)))*(2*(0:order_max)+1)).*radial_terms.*angular_terms,2)/(4*pi);

% azimuthal component
radial_terms = zeros(length(X(:)), order_max+1);
for order=0:order_max
    jn_r = sph_besselj(order, kr);
    djn_scat = dsph_besselj(order, kr_scat);
    hn_r = sph_hankel2(order, kr);
    dhn_scat = dsph_hankel2(order, kr_scat);
    radial_terms(:, order+1) = 4*pi*1i^order * (jn_r-(djn_scat./dhn_scat).*hn_r);
end
radial_terms = replicatePerOrder(radial_terms,2);
nDirs = length(azi_stack);
angular_terms = (ones(nDirs,1)*getSH(order_max,[doa(1) pi/2-doa(2)],'real')) .* ...
                get_dRSH_phi(order_max,[azi_stack pi/2-elev_stack]);
cv_phi_stack = -1/(1i*k*Z_0)* (1./(r_stack.*sin(pi/2-elev_stack))) .* sum(radial_terms.*angular_terms,2);

% polar component
radial_terms = zeros(length(X(:)), order_max+1);
for order=0:order_max
    jn_r = sph_besselj(order, kr);
    djn_scat = dsph_besselj(order, kr_scat);
    hn_r = sph_hankel2(order, kr);
    dhn_scat = dsph_hankel2(order, kr_scat);
    radial_terms(:, order+1) = 4*pi*1i^order * (jn_r-(djn_scat./dhn_scat).*hn_r);
end
radial_terms = replicatePerOrder(radial_terms,2);
nDirs = length(azi_stack);
angular_terms = (ones(nDirs,1)*getSH(order_max,[doa(1) pi/2-doa(2)],'real')) .* ...
                get_dRSH_theta(order_max,[azi_stack pi/2-elev_stack]);
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

cv_x = reshape(cv_x_stack, size(X));
cv_y = reshape(cv_y_stack, size(X));
cv_z = reshape(cv_z_stack, size(X));        
         
v_cmplx = cat(3, cv_x, cv_y, cv_z);
v_real = real(v_cmplx);
v_mag = sqrt(sum(v_real.^2,3));

% draw
v_mag(R>D/2) = Inf;
v_mag(R<r_scat) = Inf;
figure
subplot(1,2,1)
surf(X,Y,v_mag); colorbar
shading interp;
arrow1 = -1.2*D/2*n_dop; arrow2 = -D/2*n_dop; darrow = arrow2-arrow1;
hold on
quiver(arrow1(1),arrow1(2),darrow(1),darrow(2),0, 'linewidth',5,'color','k','maxheadsize',1);
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])

subplot(1,2,2)
vx_real = v_real(:,:,1); vx_real(R>D/2) = 0; vx_real(R<r_scat) = 0;
vy_real = v_real(:,:,2); vy_real(R>D/2) = 0; vy_real(R<r_scat) = 0;
vz_real = v_real(:,:,3); vz_real(R>D/2) = 0; vz_real(R<r_scat) = 0;
quiver(X,Y, vx_real, vy_real);
arrow1 = -1.2*D/2*n_dop; arrow2 = -D/2*n_dop; darrow = arrow2-arrow1;
hold on
quiver(arrow1(1),arrow1(2),darrow(1),darrow(2),0, 'linewidth',5,'color','k','maxheadsize',1);
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])
