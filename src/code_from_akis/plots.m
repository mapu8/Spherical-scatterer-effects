
% acoustics
c = 343;
rho_0 = 1.2;
Z_0 = c*rho_0;

% playback diameter for computations, in meters
D = 10;
% grid resolution, in meters
d = 0.2;
x = -D:d:D;
y = -D:d:D;
[X,Y] = meshgrid(x,y);
% spherical coords
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X)));

%%

% plane wave exp(jwt-jkr)
t = 0; % time point
f = 200; % frequency
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
theta_n = pi/2; % direction of propagation
n = [cos(theta_n) sin(theta_n)]; % unit vector of propagation
m0 = 1; % magnitude
phi0 = 0; % initial phase
a0 = m0.*exp(1i*phi0); % initial amplitude

p_cmplx = a0*exp(1i* (w*t*ones(size(X)) - k* (n(1)*X + n(2)*Y))); % complex pressure field
p_real = real(p_cmplx); % real pressure field
p_phase = angle(p_cmplx); % phase of pressure

% draw
p_real(R>D/2) = Inf;
p_phase(R>D/2) = Inf;
figure
subplot(1,2,1)
surf(X,Y,p_real); colorbar
shading interp;
arrow1 = -1.2*D/2*n; arrow2 = -D/2*n; darrow = arrow2-arrow1;
hold on
quiver(arrow1(1),arrow1(2),darrow(1),darrow(2),0, 'linewidth',5,'color','k','maxheadsize',1);
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

subplot(1,2,2)
surf(X,Y,p_phase); colorbar
shading interp;
arrow1 = -1.2*D/2*n; arrow2 = -D/2*n; darrow = arrow2-arrow1;
hold on
quiver(arrow1(1),arrow1(2),darrow(1),darrow(2),0, 'linewidth',5,'color','k','maxheadsize',1);
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

%%

% point source exp(jwt-jk||r_s-r||)/||r_s-r||
t = 0; % time point
f = 200; % frequency 
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
theta_n = pi/2; % direction of source
rs = 1.2*D/2; % distance of source
ns = rs*[cos(theta_n) sin(theta_n)]; % source position vector
m0 = 1; % magnitude
phi0 = 0; % initial phase
a0 = m0.*exp(1i*phi0); % initial amplitude

p_cmplx = a0*exp(1i* (w*t*ones(size(X)) - k* sqrt((X-ns(1)).^2 + (Y-ns(2)).^2) ) )./(4*pi*sqrt((X-ns(1)).^2 + (Y-ns(2)).^2)); % complex pressure field
p_real = real(p_cmplx); % real pressure field
p_phase = angle(p_cmplx); % phase of pressure

% draw
p_real(R>D/2) = Inf;
p_phase(R>D/2) = Inf;
figure
subplot(1,2,1)
surf(X,Y,p_real); colorbar
shading interp;
line(ns(1),ns(2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

subplot(1,2,2)
surf(X,Y,p_phase); colorbar
shading interp;
line(ns(1),ns(2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

%% plane wave expansion

% plane wave exp(jwt-jkr)
t = 0; % time point
f = 200; % frequency
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
dop = [pi/2 0]; % direction of propagation, azimuth/elevation
n_dop = [cos(dop(1))*cos(dop(2)) sin(dop(1))*cos(dop(2)) sin(dop(2))]; % unit vector of propagation
m0 = 1; % magnitude
phi0 = 0; % initial phase
a0 = m0.*exp(1i*phi0); % initial amplitude

n_doa = -n_dop; % direction of arrival
[doa(1), doa(2)] = cart2sph(n_doa(1),n_doa(2),n_doa(3));

x_stack = X(:);
y_stack = Y(:);
r_stack = R(:);
azi_stack = Azi(:);
elev_stack = Elev(:);
cos_alpha = sin(doa(2))*sin(elev_stack) + cos(doa(2)).*cos(elev_stack(:)).*cos(doa(1)*ones(size(X(:)))-azi_stack);
order_max = ceil(exp(1)*k*(D/2)/2);
radial_terms = zeros(length(X(:)), order_max+1);
kr = k*r_stack;
for order=0:order_max
    jn = sph_besselj(order, kr);
    radial_terms(:, order+1) = 4*pi*1i^order*jn;
end
for degree=0:order_max
    temp = legendre(degree,cos_alpha);
    angular_terms(:,degree+1) = temp(1,:);
end
cp_stack = sum((ones(size(X(:)))*(2*(0:order_max)+1)).*radial_terms.*angular_terms,2)/(4*pi);
p_cmplx = reshape(cp_stack, size(X));
p_real = real(p_cmplx);
p_phase = angle(p_cmplx); % phase of pressure

% draw
p_real(R>D/2) = Inf;
p_phase(R>D/2) = Inf;
figure
subplot(1,2,1)
surf(X,Y,p_real); colorbar
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
surf(X,Y,p_phase); colorbar
shading interp;
arrow1 = -1.2*D/2*n; arrow2 = -D/2*n; darrow = arrow2-arrow1;
hold on
quiver(arrow1(1),arrow1(2),darrow(1),darrow(2),0, 'linewidth',5,'color','k','maxheadsize',1);
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

%% point source expansion

% point source exp(jwt-jk||r_s-r||)/||r_s-r||
t = 0; % time point
f = 200; % frequency
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
rs = 1.2*D/2; % distance of source
doa = [pi/2 0]; % source direction, azimuth/elevation
ns = rs*[cos(doa(1))*cos(doa(2)) sin(doa(1))*cos(doa(2)) sin(doa(2))]; % source position vector
m0 = 1; % magnitude
phi0 = 0; % initial phase
a0 = m0.*exp(1i*phi0); % initial amplitude

x_stack = X(:);
y_stack = Y(:);
r_stack = R(:);
azi_stack = Azi(:);
elev_stack = Elev(:);
cos_alpha = sin(doa(2))*sin(elev_stack) + cos(doa(2)).*cos(elev_stack(:)).*cos(doa(1)*ones(size(X(:)))-azi_stack);
order_max = ceil(exp(1)*k*(D/2)/2);
radial_terms = zeros(length(X(:)), order_max+1);
kr = k*r_stack;
for order=0:order_max
    jn = sph_besselj(order, kr);
    hns = sph_hankel2(order, k*rs);
    radial_terms(:, order+1) = (-1i)*k.*jn.*hns;
end
for degree=0:order_max
    temp = legendre(degree,cos_alpha);
    angular_terms(:,degree+1) = temp(1,:);
end
cp_stack = sum((ones(size(X(:)))*(2*(0:order_max)+1)).*radial_terms.*angular_terms,2)/(4*pi);
p_cmplx = reshape(cp_stack, size(X));
p_real = real(p_cmplx);
p_phase = angle(p_cmplx); % phase of pressure

% draw
p_real(R>D/2) = Inf;
p_phase(R>D/2) = Inf;
figure
subplot(1,2,1)
surf(X,Y,p_real); colorbar
shading interp;
line(ns(1),ns(2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

subplot(1,2,2)
surf(X,Y,p_phase); colorbar
shading interp;
line(ns(1),ns(2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

%% spherical scatterer - plane wave

% plane wave exp(jwt-jkr)
t = 0; % time point
f = 200; % frequency
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
dop = [pi/2 0]; % direction of propagation, azimuth/elevation
n_dop = [cos(dop(1))*cos(dop(2)) sin(dop(1))*cos(dop(2)) sin(dop(2))]; % unit vector of propagation
m0 = 1; % magnitude
phi0 = 0; % initial phase
a0 = m0.*exp(1i*phi0); % initial amplitude

n_doa = -n_dop; % direction of arrival
[doa(1), doa(2)] = cart2sph(n_doa(1),n_doa(2),n_doa(3));

% scatterer radius
r_scat = 2;

x_stack = X(:);
y_stack = Y(:);
r_stack = R(:);
azi_stack = Azi(:);
elev_stack = Elev(:);
cos_alpha = sin(doa(2))*sin(elev_stack) + cos(doa(2)).*cos(elev_stack(:)).*cos(doa(1)*ones(size(X(:)))-azi_stack);
order_max = ceil(exp(1)*k*(D/2)/2);

radial_terms = zeros(length(X(:)), order_max+1);
kr = k*r_stack;
kr_scat = k*r_scat;
for order=0:order_max
    jn = sph_besselj(order, kr);
    jnprime = dsph_besselj(order, kr_scat);
    hn = sph_hankel2(order, kr);
    hnprime = dsph_hankel2(order, kr_scat);
    radial_terms(:, order+1) = 4*pi*1i^order * (jn-(jnprime./hnprime).*hn);
end
for degree=0:order_max
    temp = legendre(degree,cos_alpha);
    angular_terms(:,degree+1) = temp(1,:);
end
cp_stack = sum((ones(size(X(:)))*(2*(0:order_max)+1)).*radial_terms.*angular_terms,2)/(4*pi);
cp_stack(r_stack<r_scat) = Inf;
p_cmplx = reshape(cp_stack, size(X));
p_real = real(p_cmplx);
p_phase = angle(p_cmplx); % phase of pressure

% draw
p_real(R>D/2) = Inf;
p_phase(R>D/2) = Inf;
figure
subplot(1,2,1)
surf(X,Y,p_real); colorbar
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
surf(X,Y,p_phase); colorbar
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


%% spherical scatterer - points source

% plane wave exp(jwt-jkr)
t = 0; % time point
f = 1000; % frequency
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
rs = 1.2*D/2; % distance of source
doa = [pi/2 0]; % source direction, azimuth/elevation
ns = rs*[cos(doa(1))*cos(doa(2)) sin(doa(1))*cos(doa(2)) sin(doa(2))]; % source position vector
m0 = 1; % magnitude
phi0 = 0; % initial phase
a0 = m0.*exp(1i*phi0); % initial amplitude

% scatterer radius
r_scat = 1;

x_stack = X(:);
y_stack = Y(:);
r_stack = R(:);
azi_stack = Azi(:);
elev_stack = Elev(:);
cos_alpha = sin(doa(2))*sin(elev_stack) + cos(doa(2)).*cos(elev_stack(:)).*cos(doa(1)*ones(size(X(:)))-azi_stack);
order_max = ceil(exp(1)*k*(D/2)/2);

radial_terms = zeros(length(X(:)), order_max+1);
kr = k*r_stack;
kr_scat = k*r_scat;
for order=0:order_max
    jn = sph_besselj(order, kr);
    jnprime = dsph_besselj(order, kr_scat);
    hn = sph_hankel2(order, kr);
    hnprime = dsph_hankel2(order, kr_scat);
    hns = sph_hankel2(order, k*rs);
    radial_terms(:, order+1) = (jn-(jnprime./hnprime).*hn) * hns;
end
for degree=0:order_max
    temp = legendre(degree,cos_alpha);
    angular_terms(:,degree+1) = temp(1,:);
end
cp_stack = (-1i)*k*sum((ones(size(X(:)))*(2*(0:order_max)+1)).*radial_terms.*angular_terms,2)/(4*pi);
cp_stack(r_stack<r_scat) = Inf;
p_cmplx = reshape(cp_stack, size(X));
p_real = real(p_cmplx);
p_phase = angle(p_cmplx); % phase of pressure

% draw
p_real(R>D/2) = Inf;
p_phase(R>D/2) = Inf;
figure
subplot(1,2,1)
surf(X,Y,p_real); colorbar
shading interp;
line(ns(1),ns(2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])

subplot(1,2,2)
surf(X,Y,p_phase); colorbar
shading interp;
line(ns(1),ns(2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])