%% Analyze pressure, particle velocity and intensity of recreated plane waves.

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Lauros Pajunen, 2018
%   Department of Signal Processing and Acoustics, Aalto University, Finland
%   lauros.pajunen@aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc;clear all;%close all;

%% Initializations
f = 200;
c = 343;
rho0 = 1.225;
z0 = rho0 * c;
k = 2*pi*f/c;
t1 = 1;
t2 = t1 + 100;
x = linspace(0,10,100);
[X,Y]=meshgrid(x,x);
A = 1;
B = 1;
p1 = A*cos(-k*X);
p2 = B*cos(-k*Y);
%p1 = A*exp(1i*(-k*X));
%p2 = A*exp(1i*(-k*Y));

%% Velocity and intesity by dividing with impedance z0

% Velocities
u = p1/z0;
v = p2/z0;

% Complex intensities
C1x = (1/(2*z0))*(A^2 + A*B*exp(1i*k*(X)));
C2x = (1/(2*z0))*(A^2 + A*B*exp(1i*k*(X-Y)));
C2y = (1/(2*z0))*(B^2 + A*B*exp(1i*k*(Y-X)));
% Active intensities (p. 62 of Impedance book by Fahy)
I1x = real(C1x);
I2x = real(C2x);
I2y = real(C2y);

% Test
p11 = A*exp(1i*(-k*X));
p22 = B*exp(1i*(-k*Y));
p33 = p11 + p22;
u11 = p11/z0;
u22 = p22/z0;
C11 = (1/(2))*p33.*conj(u11);
C22 = (1/(2))*p33.*conj(u22);
C33 = (1/(2*z0))*p33*conj(p33);
I11 = real(C11);
I22 = real(C22);
%[I11x, I11y] = pol2cart(deg2rad(0),real(C2x));
%[I22x, I22y] = pol2cart(deg2rad(90),real(C2y));
%I11 = I11x + I22x;
%I22 = I11y + I22y;
%[I11, I22] = pol2cart(deg2rad(45),real(C33));



%% Plot
figure(1);
surf(X,Y,p1);
view(2)
colorbar
title("1 plane wave propagating +x direction, 'pressure'");
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(2);
surf(X,Y,p1+p2);
view(2)
colorbar
title("2 plane waves, +x and +y directions, 'pressure'");
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(3);
quiver(u,zeros(max(size(X))))
title("1 plane wave propagating +x direction, particle velocity");
axis([0 100 0 100])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(4);
quiver(u,v)
title("2 plane waves, +x and +y directions, particle velocity");
axis([0 100 0 100])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(5);
quiver(u*p1,zeros(max(size(X))))
title("1 plane wave propagating +x direction, 'intensity'");
axis([0 100 0 100])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(6);
quiver(u*p1,v*p2)
title("2 plane waves, +x and +y directions, 'intensity'");
axis([0 100 0 100])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(7);
quiver(I1x,zeros(max(size(X))))
title("1 plane wave propagating +x direction, 'active intensity'");
axis([0 100 0 100])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(8);
quiver(I2x,I2y)
title("2 plane waves, +x and +y directions, 'active intensity'");
axis([0 100 0 100])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(30);
quiver(I11,I22)
title("2 plane waves, +x and +y directions I11, 'active intensity'");
axis([0 100 0 100])
%set(gcf, 'Position', [2000, -200, 1000, 800])

%% Trying to plot pressure, impedance, velocity and energy
% (p.52 in Intensity book by Fahy)
pressure = cos(-k*x + 2*pi*f*t1);
impedance = rho0 * c;
velocity = pressure / impedance;
energy = 0.5*rho0*(velocity.^2) + (pressure.^2)/(2*rho0*(c^2));
intensity = pressure .* velocity;

figure(9);
hold on;
h1 = plot(pressure,'b');
h2 = plot(velocity,'r');
h3 = plot(intensity,'black');
h4 = plot(energy,'c');
legend([h1,h2,h3,h4], "pressure","velocity","intensity","energy")
axis([20 40 -3e-3 3e-3])
%set(gcf, 'Position', [2000, -200, 1000, 800])
hold off;


%% Recreation of a plane wave by two plane waves

% Meshgrid
xplane = linspace(0,10,300);
[Xplane,Yplane]=meshgrid(xplane,xplane);
% Frequency
fplane = 200;
% Amplitudes
C = 2.32;

% Calculate the tilted plane waves
frequencies = ones(1,2)*fplane;
amplitudes = ones(1,2)*C;
directions = [30 -30];
plane_waves = tilted_plane_waves(Xplane,Yplane,frequencies,amplitudes,directions,true);
p3c = plane_waves(:,:,1);
p4c = plane_waves(:,:,2);

% Real parts
p3 = real(p3c);
p4 = real(p4c);

% Combine waves
p5c = p3c + p4c;
p5 = p3 + p4;

% Velocities
u1_planec = p1_planec/z0;
v1_planec = zeros(max(size(u1_planec)));
u1_plane = real(u1_planec);
v1_plane = ones(size(u1_plane))*1e-6;
[u3c,v3c] = pol2cart(deg2rad(dir1),p3c/z0);
[u4c,v4c] = pol2cart(deg2rad(dir2),p4c/z0);
u3 = real(u3c);
v3 = real(v3c);
u4 = real(u4c);
v4 = real(v4c);
u5c = u3c + u4c;
v5c = v3c + v4c;
u5 = u3 + u4;
v5 = v3 + v4;
repr_velocity = u5 + v5;

% Complex intensities
C1x_plane = (1/2)*p1_planec.*conj(u1_planec);
C1y_plane = (1/2)*p1_planec.*conj(v1_plane);
C5x = (1/2)*p5c.*conj(u5c);
C5y = (1/2)*p5c.*conj(v5c);

% Active intensities
I1x_plane = real(C1x_plane);
I1y_plane = ones(size(I1x_plane))*1e-6;
I5x_plane = real(C5x);
I5y_plane = real(C5y);

% Measure reproduction error
% Center of the reference point in the original plane wave
ref_center = [150,52];
% Center of the reference point in the recreated plane wave
repr_center = [150, 80];
pressure_at_recreation_center = p5(repr_center(1),repr_center(2))
% Area of error measurement
errorx = 30;
errory = 30;
% Shift the reproduced wave matrix to make comparison easier
shift_amnt = ref_center(2) - repr_center(2);
shifted_repr_pressure = circshift(p5, shift_amnt, 2);
shifted_repr_u = circshift(u5, shift_amnt, 2);
shifted_repr_v = circshift(v5, shift_amnt, 2);
shifted_repr_I5x = circshift(I5x_plane, shift_amnt, 2);
shifted_repr_I5y = circshift(I5y_plane, shift_amnt, 2);
% Calculate the error (squared)
%error_matrix_pressure = (shifted_repr_pressure - p1_plane).^2;   
%error_matrix_velocity = abs(shifted_repr_u - u1_plane).^2 + abs(shifted_repr_v - v1_plane).^2;
%error_matrix_intensity = abs(shifted_repr_Ix - I1x_plane).^2 + abs(shifted_repr_Iy - I1y_plane).^2;
% Percentage
error_matrix_pressure = abs((shifted_repr_pressure - p1_plane)./(p1_plane));   
error_matrix_velocity = abs((shifted_repr_u - u1_plane)./u1_plane) + abs(shifted_repr_v);
error_matrix_intensity = abs((shifted_repr_I5x - I1x_plane)./I1x_plane) + abs(shifted_repr_I5y);
% Velocity error
% Select an area from the error matrix
selected_error_pressure = error_matrix_pressure([(ref_center(1)-errorx):(ref_center(1)+errorx)],...
                                                [(ref_center(2)-errory):(ref_center(2)+errory)]);
selected_error_velocity = error_matrix_velocity([(ref_center(1)-errorx):(ref_center(1)+errorx)],...
                                                [(ref_center(2)-errory):(ref_center(2)+errory)]);
selected_error_intensity = error_matrix_intensity([(ref_center(1)-errorx):(ref_center(1)+errorx)],...
                                                  [(ref_center(2)-errory):(ref_center(2)+errory)]);

selected_size = max(size(selected_error_pressure));
center = ceil(selected_size/2);

% Plot
% Percentage limits in contour plots
error_limits = [0 0.05 0.1 0.2 0.5 0.7 1];

figure(10);
surf(Xplane,Yplane,p3);
view(2)
colorbar
caxis([-1,1])
title(sprintf("Plane wave, f=%d Hz, propagation angle=%d",fplane,dir1))

figure(11);
surf(Xplane,Yplane,p4);
view(2)
colorbar
caxis([-1,1])
title(sprintf("Plane wave, f=%d Hz, propagation angle=%d",fplane,dir2))

figure(12);
surf(Xplane,Yplane,p5);
view(2)
colorbar
caxis([-1,1])
title(sprintf("Combination of plane waves"))

figure(13);
%contourf(error_matrix_pressure, error_limits, 'ShowText', 'on');
surf(Xplane,Yplane,error_matrix_pressure);
view(2)
colorbar
%caxis([-1,1])
caxis([0,1])
title(sprintf("Error pressure"))

figure(14);
surf(Xplane,Yplane,error_matrix_velocity);
view(2)
colorbar
%caxis([-1,1])
caxis([0,1])
title(sprintf("Error velocity"))

figure(15);
surf(Xplane,Yplane,error_matrix_intensity);
view(2)
colorbar
%caxis([0,0.000001])
caxis([0,1])
title(sprintf("Error intensity"))

%close(figure(16));
figure(16);
[Ccont, hcont] = contourf(selected_error_pressure, error_limits, 'ShowText', 'on');
clabel(Ccont,hcont,'Color','w')
grid on
hold on
%surf(selected_error_pressure);
view(2)
plot3(center,center,100,'rx', 'LineWidth', 3, 'MarkerSize',10)
plot3([50 60],[4 4],[100 100], 'LineWidth',2,'Color','black')
plot3([30 60],[2 2],[100 100], 'LineWidth',2,'Color','black')
text(50,5,'0.33m')
text(30,3,'1m')
viscircles([center center],30)
hold off
colorbar
%caxis([0,0.1])
caxis([0,1])
title(sprintf("Error pressure, zoomed in"))
%axis([0 0 selected_size selected_size])

%close(figure(17));
figure(17);
[Ccont, hcont] = contourf(selected_error_velocity, error_limits, 'ShowText', 'on');
clabel(Ccont,hcont,'Color','w')
grid on
hold on
%surf(selected_error_velocity);
view(2)
plot3(center,center,100,'rx', 'LineWidth', 3, 'MarkerSize',10)
plot3([50 60],[4 4],[100 100], 'LineWidth',2,'Color','black')
plot3([30 60],[2 2],[100 100], 'LineWidth',2,'Color','black')
text(50,5,'0.33m')
text(30,3,'1m')
viscircles([center center],30)
hold off
colorbar
%caxis([0,0.01])
caxis([0,1])
title(sprintf("Error velocity, zoomed in"))
%axis([0 0 selected_size selected_size])

%close(figure(18));
figure(18);
[Ccont, hcont] = contourf(selected_error_intensity, error_limits, 'ShowText', 'on');
clabel(Ccont,hcont,'Color','w')
grid on
hold on
%surf(selected_error_intensity);
view(2)
plot3(center,center,100,'rx', 'LineWidth', 3, 'MarkerSize',10)
plot3([50 60],[4 4],[100 100], 'LineWidth',2,'Color','black')
plot3([30 60],[2 2],[100 100], 'LineWidth',2,'Color','black')
text(50,5,'0.33m')
text(30,3,'1m')
viscircles([center center],30)
hold off
colorbar
%caxis([0,0.000001])
caxis([0,1])
title(sprintf("Error intensity, zoomed in"))
%axis([0 0 selected_size selected_size])

figure(19);
surf(Xplane,Yplane,p1_plane);
view(2)
colorbar
title("1 plane wave propagating +x direction, 'pressure'");

figure(20);
quiver(u3,v3)
title("1 plane wave, angled direction 1, particle velocity");
axis([0 max(size(Xplane)) 0 max(size(Yplane))])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(21);
quiver(u4,v4)
title("1 plane wave, angled direction 2, particle velocity");
axis([0 max(size(Xplane)) 0 max(size(Yplane))])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(22);
quiver(u5,v5)
title("2 plane waves, angled directions, particle velocity");
axis([0 max(size(Xplane)) 0 max(size(Yplane))])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(23);
quiver(u5*u5*z0,v5*v5*z0)
title("2 plane waves, angled directions, 'intensity'");
axis([0 max(size(Xplane)) 0 max(size(Yplane))])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(24);
%quiver(Ix_plane,Iy_plane)
quiver(I5x_plane,I5y_plane)
title("2 plane waves, angled directions, 'active intensity'");
axis([0 max(size(Xplane)) 0 max(size(Yplane))])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(25);
quiver(real(u1_plane),zeros(max(size(Xplane))))
title("1 plane wave propagating +x direction, particle velocity");
axis([0 max(size(Xplane)) 0 max(size(Yplane))])
%set(gcf, 'Position', [2000, -200, 1000, 800])

figure(26);
quiver(I1x_plane,I1y_plane)
title("1 plane wave propagating +x direction, active intensity");
axis([0 max(size(Xplane)) 0 max(size(Yplane))])
%set(gcf, 'Position', [2000, -200, 1000, 800])

%% Spherical waves

f2 = 200;
k2 = 2*pi*f2/c;
E = 1;
F = 1;
% Center of recreation (meters away from the sources)
center = 2;
% Angle between the two sources
angle = 30;
% Position of the source in the meshgrid (10x10)
Xoff1 = 0;
Yoff1 = 4;
Xoff2 = 0;
Yoff2 = Yoff1 + 2*center*tan((angle/2)*pi/180);
%Xoff3 = 0;
%Yoff3 = Yoff2 + 2*center*tan((angle/2)*pi/180);
%Xoff4 = 0;
%Yoff4 = Yoff3 + 2*center*tan((angle/2)*pi/180);
r = sqrt((X-Xoff1).^2 + (Y-Yoff1).^2);
r = rescale(r,0,10);
r2 = sqrt((X-Xoff2).^2 + (Y-Yoff2).^2);
r2 = rescale(r2,0,10);
%r3 = sqrt((X-Xoff3).^2 + (Y-Yoff3).^2);
%r3 = rescale(r3,0,10);
%r4 = sqrt((X-Xoff4).^2 + (Y-Yoff4).^2);
%r4 = rescale(r4,0,10);
rsmall = r([51:100],[51:100]);
sp1 = E*exp(-1i*k2*r)./(r+1);
sp2 = -F*exp(-1i*k2*r2)./(r2+1);
%sp3 = F*exp(-1i*k2*r3)./(r3+1);
%sp4 = F*exp(-1i*k2*r4)./(r4+1);
%sp5 = (sp1 + sp2 + sp3 + sp4)/4;
sp3 = sp1 + sp2;
sp3real = real(sp3);
spsmall = E*exp(-1i*k2*rsmall)./(rsmall+1);

% Plot
figure(30);
surf(X,Y,real(sp1));
view(2)
colorbar

figure(31);
surf(X,Y,real(sp2));
view(2)
colorbar

figure(32);
surf(X,Y,real(sp3));
view(2)
colorbar


%% OLD stuff, velocities incorrect
% (tried to calculate velocities by substracting pressures between two 
%  time steps)
%{
%% Pressures
% One plane wave
p1_t1 = cos(-k*X + 2*pi*f*t1);
p1_t2 = cos(-k*X + 2*pi*f*t2);
% Two plane waves
p21_t1 = cos(-k*X + 2*pi*f*t1);
p22_t1 = cos(-k*Y + 2*pi*f*t1);
p21_t2 = cos(-k*X + 2*pi*f*t2);
p22_t2 = cos(-k*Y + 2*pi*f*t2);
p2_t1 = p21_t1 + p22_t1;
p2_t2 = p21_t2 + p22_t2;

%% Pressure velocities
point = [20,20];
x_basis = [1,0];
y_basis = [0,1];
u1_x = (p1_t2(point(1),point(2)) - p1_t1(point(1),point(2)))*x_basis
u1_y = [0,0]
u1 = u1_x + u1_y

% One plane wave
u11 = p1_t2 - p1_t1;

% Two plane waves
u2x = p21_t2 - p21_t1;
u2y = p22_t2 - p22_t1;
%}









