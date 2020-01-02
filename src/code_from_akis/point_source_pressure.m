
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
% spherical coords
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X)));

x_stack = X(:);
y_stack = Y(:);
r_stack = R(:);
azi_stack = Azi(:);
elev_stack = Elev(:);

%% signal parameters

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

kr = k*r_stack;

%% direct modeling

% point source exp(jwt-jk||r_s-r||)/||r_s-r||
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

%% point source expansion

% point source exp(jwt-jk||r_s-r||)/||r_s-r||
cos_alpha = sin(doa(2))*sin(elev_stack) + cos(doa(2)).*cos(elev_stack(:)).*cos(doa(1)*ones(size(X(:)))-azi_stack);
order_max = ceil(exp(1)*k*(D/2)/2);
radial_terms = zeros(length(X(:)), order_max+1);
for order=0:order_max
    jn = sph_besselj(order, kr);
    hns = sph_hankel2(order, k*rs);
    radial_terms(:, order+1) = (-1i)*k.*jn.*hns;
end
angular_terms = zeros(length(X(:)), order_max+1);
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

%% spherical scatterer - points source

% scatterer radius
r_scat = 1;
kr_scat = k*r_scat;

cos_alpha = sin(doa(2))*sin(elev_stack) + cos(doa(2)).*cos(elev_stack(:)).*cos(doa(1)*ones(size(X(:)))-azi_stack);
order_max = ceil(exp(1)*k*(D/2)/2);

radial_terms = zeros(length(X(:)), order_max+1);
for order=0:order_max
    jn = sph_besselj(order, kr);
    jnprime = dsph_besselj(order, kr_scat);
    hn = sph_hankel2(order, kr);
    hnprime = dsph_hankel2(order, kr_scat);
    hns = sph_hankel2(order, k*rs);
    radial_terms(:, order+1) = (jn-(jnprime./hnprime).*hn) * hns;
end
angular_terms = zeros(length(X(:)), order_max+1);
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