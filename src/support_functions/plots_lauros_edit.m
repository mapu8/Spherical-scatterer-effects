clear;

% acoustics
c = 343;
rho_0 = 1.2;
Z_0 = c*rho_0;

% playback diameter for computations, in meters
D = 4;
% grid resolution, in meters
d = 0.1;
x = -D:d:D;
y = -D:d:D;
[X,Y] = meshgrid(x,y);
% spherical coords
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X)));

%% point source expansion - pressure

t_vector = linspace(0,1,9); % time vector
frequencies = ones(1,1)*200;
source_positions = [[pi/2 0 1.2*D/2]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
[p_cmplx_original, ns_vector] = point_sources_pressure(t_vector(1),frequencies,source_positions,magnitudes,initial_phases,X,Y,D);
p_real_original = real(p_cmplx_original);
p_phase_original = angle(p_cmplx_original); % phase of pressure

% draw
p_real_original(R>D/2) = Inf;
p_phase_original(R>D/2) = Inf;
figure
subplot(1,2,1)
surf(X,Y,p_real_original); colorbar
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Pressure, f=%d Hz",frequencies(1)));

subplot(1,2,2)
surf(X,Y,p_phase_original); colorbar
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Phase, f=%d Hz",frequencies(1)));


%% spherical scatterer - points source - pressure

t = 0;
r_scat = 0.2; % radius of scatterer
frequencies = ones(1,3)*1000;
source_positions = [[pi/2 0 1.2*D/2];[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1 1 1];
initial_phases = [0 0 0];
[p_cmplx, ns_vector] = point_sources_scatterer_pressure(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,D,r_scat);
p_real = real(p_cmplx);
p_phase = angle(p_cmplx); % phase of pressure

% draw
p_real(R>D/2) = Inf;
p_phase(R>D/2) = Inf;
figure
subplot(1,2,1)
surf(X,Y,p_real); colorbar
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Pressure, f=%d Hz,\n scatterer radius=%.2f m",frequencies(1),r_scat));
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])

subplot(1,2,2)
surf(X,Y,p_phase); colorbar
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Phase, f=%d Hz,\n scatterer radius=%.2f m",frequencies(1),r_scat));
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])


%% point source expansion - velocity

t_vector = linspace(0,1,9); % time vector
frequencies = ones(1,1)*200;
source_positions = [[pi/2 0 1.2*D/2]];%[pi/4 0 1.2*D/2]];%[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
[v_cmplx_original, ns_vector] = point_sources_velocity(t_vector(1),frequencies,source_positions,magnitudes,initial_phases,X,Y,D);
v_real_original = real(v_cmplx_original);
v_mag_original = sqrt(sum(v_real_original.^2,3));

% draw
v_mag_original(R>D/2) = Inf;
figure
subplot(1,2,1)
surf(X,Y,v_mag_original); colorbar
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);

subplot(1,2,2)
vx_real = v_real_original(:,:,1); vx_real(R>D/2) = 0;
vy_real = v_real_original(:,:,2); vy_real(R>D/2) = 0;
vz_real = v_real_original(:,:,3); vz_real(R>D/2) = 0;
quiver(X,Y, vx_real, vy_real);
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);


%% spherical scatterer - points source - velocity

t_vector = linspace(0,1,9); % time vector
r_scat = 0.5; % radius of scatterer
frequencies = ones(1,2)*200;
source_positions = [[pi/2 0 1.2*D/2];[pi/4 0 1.2*D/2]]%;[3*pi/4 0 1.2*D/2]];
magnitudes = [1 1];
initial_phases = [0 0];
[v_cmplx, ns_vector] = point_sources_scatterer_velocity(t_vector(1),frequencies,source_positions,magnitudes,initial_phases,X,Y,D,r_scat);
v_real = real(v_cmplx);
v_mag = sqrt(sum(v_real.^2,3));

% draw
v_mag(R>D/2) = Inf;
v_mag(R<r_scat) = Inf;
figure
subplot(1,2,1)
surf(X,Y,v_mag); colorbar
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
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
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])


%% point source expansion - intensity

t_vector = linspace(0,1,9); % time vector
f = 200;
frequencies = ones(1,1)*f;
source_positions = [[pi/2 0 1.2*D/2]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
[p_cmplx, ns_vector] = point_sources_pressure(t_vector(1),frequencies,source_positions,magnitudes,initial_phases,X,Y,D);
[v_cmplx, ns_vector] = point_sources_velocity(t_vector(1),frequencies,source_positions,magnitudes,initial_phases,X,Y,D);
[i_cmplx] = point_sources_intensity(p_cmplx,v_cmplx);

i_real = real(i_cmplx);
i_mag = sqrt(sum(i_real.^2,3));

% draw
i_mag(R>D/2) = Inf;
figure
set(gcf,'Position',[100 100 2000 2000])
subplot(1,2,1)
surf(X,Y,i_mag); colorbar
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Magnitude of intensity, f=%d Hz",frequencies(1)));
xlabel("m");ylabel("m");


subplot(1,2,2)
ix_real = i_real(:,:,1); ix_real(R>D/2) = 0;
iy_real = i_real(:,:,2); iy_real(R>D/2) = 0;
iz_real = i_real(:,:,3); iz_real(R>D/2) = 0;
quiver(X,Y, ix_real, iy_real);
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Active intensity, f=%d Hz",frequencies(1)));
xlabel("m");ylabel("m");


%% spherical scatterer - points source - intensity

r_scat = 0.5; % radius of scatterer
f = 200;
t_vector = linspace(0,1/f,9); % time vector
frequencies = ones(1,1)*f;
source_positions = [[pi/2 0 1.2*D/2]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
[p_cmplx, ns_vector] = point_sources_scatterer_pressure(t_vector(1),frequencies,source_positions,magnitudes,initial_phases,X,Y,D,r_scat);
[v_cmplx, ns_vector] = point_sources_scatterer_velocity(t_vector(1),frequencies,source_positions,magnitudes,initial_phases,X,Y,D,r_scat);
[i_cmplx] = point_sources_intensity(p_cmplx,v_cmplx);

i_real = real(i_cmplx);
i_mag = sqrt(sum(i_real.^2,3));

% draw
i_mag(R>D/2) = Inf;
i_mag(R<r_scat) = Inf;
figure
set(gcf,'Position',[100 100 2000 2000])
subplot(1,2,1)
surf(X,Y,i_mag); colorbar
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Magnitude of intensity, f=%d Hz,\n scatterer radius=%.2f m",frequencies(1),r_scat));
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])
xlabel("m");ylabel("m");
%caxis([0 6e-6]);

subplot(1,2,2)
ix_real = i_real(:,:,1); ix_real(R>D/2) = 0; ix_real(R<r_scat) = 0;
iy_real = i_real(:,:,2); iy_real(R>D/2) = 0; iy_real(R<r_scat) = 0;
iz_real = i_real(:,:,3); iz_real(R>D/2) = 0; iz_real(R<r_scat) = 0;
quiver(X,Y, ix_real, iy_real);
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
axis([-1.2*D/2 1.2*D/2 -1.2*D/2 1.2*D/2 -100 100]);
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Active intensity, f=%d Hz,\n scatterer radius=%.2f m",frequencies(1),r_scat));
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])
xlabel("m");ylabel("m");



