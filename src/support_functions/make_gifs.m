%% Plot and save gifs and images

clear; clc;

% playback radius for computations, in meters
D = 1;
% grid resolution, in meters
d = 0.01;
x = -D:d:D;
y = -D:d:D;
[X,Y] = meshgrid(x,y);
% spherical coords
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X)));

c = 343;

% Plot properties
polar_font_size = 22;
plot_font_size = 22;
suptitle_font_size = 30;
title_font_size = 22;
legend_font_size = 22;

%% Point source gif - pressure

f = 200;
w = 2*pi*f;
k = w/c;
temp = linspace(0,1/f,16); t_vector = temp(1:length(temp)-1); % time vector
frequencies = ones(1,1)*f;
source_positions = [[pi/2 0 1.2*D]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
N = ceil(exp(1)*k*(D)/2); % maximum order of expansion

h = figure;
set(gcf,'Position',[100 100 2000 2000])
filename = sprintf("gifs/pressure_f%d.gif",frequencies(1));
for t = t_vector
    [p_cmplx, ns_vector] = point_sources_pressure(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D);
    p_real = real(p_cmplx);
    p_phase = angle(p_cmplx); % phase of pressure
    
    % draw
    p_real(R>D) = Inf;
    p_phase(R>D) = Inf;
    subplot(1,2,1)
    surf(X,Y,p_real); colorbar
    shading interp;
    % Draw sources
    for i=1:size(ns_vector,1)
        line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
    end
    axis equal
    axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
    view(2)
    caxis([-0.1,0.1])
    line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
    line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
    title(sprintf("Pressure, f=%d Hz",frequencies(1)));

    subplot(1,2,2)
    surf(X,Y,p_phase); colorbar
    shading interp;
    % Draw sources
    for i=1:size(ns_vector,1)
        line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
    end
    axis equal
    axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
    caxis([-3,3])
    view(2)
    line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
    line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
    title(sprintf("Phase, f=%d Hz",frequencies(1)));
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if t == 0 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
    end 
end

%% Point source gif with scatterer - pressure

f = 1000;
w = 2*pi*f;
k = w/c;
temp = linspace(0,1/f,16); t_vector = temp(1:length(temp)-1); % time vector
frequencies = ones(1,1)*f;
source_positions = [[pi/2 0 1.2*D]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
r_scat = 0.2; % radius of scatterer
N = ceil(exp(1)*k*(D)/2); % maximum order of expansion

h = figure;
set(gcf,'Position',[100 100 2000 2000])
filename = sprintf("gifs/pressure_f%d_r%.2f.gif",frequencies(1),r_scat);
for t = t_vector
    [p_cmplx, ns_vector] = point_sources_scatterer_pressure(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D,r_scat);
    p_real = real(p_cmplx);
    p_phase = angle(p_cmplx); % phase of pressure
    
    % draw
    p_real(R>D) = Inf;
    p_phase(R>D) = Inf;
    subplot(1,2,1)
    surf(X,Y,p_real); colorbar
    shading interp;
    % Draw sources
    for i=1:size(ns_vector,1)
        line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
    end
    axis equal
    axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
    view(2)
    caxis([-0.1,0.1])
    line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
    line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
    % draw scatterer
    rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])
    title(sprintf("Pressure, f=%d Hz,\n scatterer radius=%.2f m",frequencies(1),r_scat));

    subplot(1,2,2)
    surf(X,Y,p_phase); colorbar
    shading interp;
    % Draw sources
    for i=1:size(ns_vector,1)
        line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
    end
    axis equal
    axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
    caxis([-3,3])
    view(2)
    line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
    line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
    % draw scatterer
    rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])
    title(sprintf("Phase, f=%d Hz,\n scatterer radius=%.2f m",frequencies(1),r_scat));
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if t == 0 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
    end 
end


%% Point source gif - velocity

f = 200;
w = 2*pi*f;
k = w/c;
temp = linspace(0,1/f,16); t_vector = temp(1:length(temp)-1); % time vector
frequencies = ones(1,4)*f;
source_positions = [[pi/4 0.1959*pi 10*D];[-pi/4 -0.1959*pi 10*D];[3*pi/4 -0.1959*pi 10*D/2];[-3*pi/4 0.1959*pi 10*D]];
magnitudes = [0.0013 0.4987 0.0013 0.4987];
initial_phases = [0 0 0 0];
N = ceil(exp(1)*k*(D)/2); % maximum order of expansion

h = figure;
set(gcf,'Position',[100 100 2000 2000])
filename = sprintf("Thesis/gifs/test4_velocity_f%d.gif",frequencies(1));
for t = t_vector
    [v_cmplx_original, ns_vector] = point_sources_velocity(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D);
    v_real_original = real(v_cmplx_original);
    v_mag_original = sqrt(sum(v_real_original.^2,3));

    % draw
    v_mag_original(R>D) = Inf;
    subplot(1,2,1)
    surf(X,Y,v_mag_original); colorbar
    shading interp;
    % Draw sources
    for i=1:size(ns_vector,1)
        line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
    end
    axis equal
    axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
    caxis([0,3e-4])
    view(2)
    line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
    line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
    title(sprintf("Magnitude of velocity, f=%d Hz",frequencies(1)));

    subplot(1,2,2)
    vx_real = v_real_original(:,:,1); vx_real(R>D) = 0;
    vy_real = v_real_original(:,:,2); vy_real(R>D) = 0;
    vz_real = v_real_original(:,:,3); vz_real(R>D) = 0;
    quiver(X,Y, vx_real, vy_real,2);
    %quiver(X(1:2:end,1:2:end),Y(1:2:end,1:2:end),vx_real(1:2:end,1:2:end),vy_real(1:2:end,1:2:end))
    % Draw sources
    for i=1:size(ns_vector,1)
        line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
    end
    axis equal
    axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
    view(2)
    line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
    line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
    title(sprintf("Velocity, f=%d Hz",frequencies(1)));
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if t == 0 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
    end 
end


%% Point source gif with scatterer - velocity

f = 1000;
w = 2*pi*f;
k = w/c;
temp = linspace(0,1/f,16); t_vector = temp(1:length(temp)-1); % time vector
frequencies = ones(1,1)*f;
source_positions = [[pi/2 0 1.2*D]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
r_scat = 0.2; % radius of scatterer
N = ceil(exp(1)*k*(D)/2); % maximum order of expansion

h = figure;
set(gcf,'Position',[100 100 2000 2000])
filename = sprintf("gifs/velocity_f%d_r%.2f.gif",frequencies(1),r_scat);
for t = t_vector
    [v_cmplx, ns_vector] = point_sources_scatterer_velocity(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D,r_scat);
    v_real = real(v_cmplx);
    v_mag = sqrt(sum(v_real.^2,3));

    % draw
    v_mag(R>D) = Inf;
    v_mag(R<r_scat) = Inf;
    subplot(1,2,1)
    surf(X,Y,v_mag); colorbar
    shading interp;
    % Draw sources
    for i=1:size(ns_vector,1)
        line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
    end
    axis equal
    axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
    caxis([0,3e-4])
    view(2)
    line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
    line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
    title(sprintf("Magnitude of velocity, f=%d Hz,\n scatterer radius=%.2f m",frequencies(1),r_scat));
    % draw scatterer
    rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])

    subplot(1,2,2)
    vx_real = v_real(:,:,1); vx_real(R>D) = 0; vx_real(R<r_scat) = 0;
    vy_real = v_real(:,:,2); vy_real(R>D) = 0; vy_real(R<r_scat) = 0;
    vz_real = v_real(:,:,3); vz_real(R>D) = 0; vz_real(R<r_scat) = 0;
    quiver(X,Y, vx_real, vy_real,2);
    %quiver(X(1:2:end,1:2:end),Y(1:2:end,1:2:end),vx_real(1:2:end,1:2:end),vy_real(1:2:end,1:2:end))
    % Draw sources
    for i=1:size(ns_vector,1)
        line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
    end
    axis equal
    axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
    view(2)
    line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
    line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
    title(sprintf("Velocity, f=%d Hz,\n scatterer radius=%.2f m",frequencies(1),r_scat));
    % draw scatterer
    rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if t == 0 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
    end 
end


%% Point source image - pressure

f = 500;
t = 0;
w = 2*pi*f;
k = w/c;
image_path = "Thesis/images/";
frequencies = ones(1,1)*f;
source_positions = [[pi/2 0 2]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];

set(groot,'defaultFigurePosition',[2000 1000 1500 700]);
set(groot,'defaultFigurePaperType','a4');
set(groot,'defaultFigurePaperOrientation','portrait');

h = figure;
filename = sprintf("pressure_f%d",frequencies(1));

% Calculations
N = ceil(exp(1)*k*(D)/2); % maximum order of expansion
[p_cmplx, ns_vector] = point_sources_pressure(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D);
%[p_cmplx] = plane_wave_expansion_pressure(t,frequencies,magnitudes,initial_phases,[-pi/2 pi/2],X,Y,N,D);
p_real = real(p_cmplx);
p_phase = angle(p_cmplx); % phase of pressure

% Plot
ax = axes;
hsup = suptitle(sprintf("Pressure field from a point source, f = %d Hz",frequencies(1)));
set(hsup,'FontSize',suptitle_font_size);
subplot(1,2,1,ax)
surf(X,Y,p_real);
ax.FontSize = 20;
bar = colorbar;
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
axis([-1 1 -1 1 -100 100]);
view(2)
caxis([-0.1,0.1])
%line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
%line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Amplitude"),'FontSize',title_font_size);
xlabel("m"); ylabel("m");
xticks(-1:0.5:1); yticks(-1:0.5:1);

ax = axes;
subplot(1,2,2,ax)
surf(X,Y,p_phase*180/pi);
ax.FontSize = 20;
bar = colorbar;
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -180 180]);
axis([-1 1 -1 1 -180 180]);
caxis([-180,180])
view(2)
%line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
%line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Phase"),'FontSize',title_font_size);
xlabel("m"); ylabel("m"); ylabel(bar, "Degrees");
xticks(-1:0.5:1); yticks(-1:0.5:1);

% Save image and reset default figure properties
print(h, image_path + filename, '-dsvg')
set(groot,'defaultFigurePosition','remove');
set(groot,'defaultFigurePaperType','remove');
set(groot,'defaultFigurePaperOrientation','remove');


%% Point source image with scatterer - pressure

f = 1000;
r_scat = 0.2;
t = 0;
w = 2*pi*f;
k = w/c;
image_path = "Thesis/images/";
frequencies = ones(1,1)*f;
source_positions = [[pi/2 0 2]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
filename = sprintf("pressure_f%d_r%d",frequencies(1),r_scat*100);

set(groot,'defaultFigurePosition',[2000 1000 1500 700]);
set(groot,'defaultFigurePaperType','a4');
set(groot,'defaultFigurePaperOrientation','portrait');

% Calculations
N = ceil(exp(1)*k*(D)/2); % maximum order of expansion
[p_cmplx, ns_vector] = point_sources_scatterer_pressure(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D,r_scat);
p_cmplx(R<r_scat) = NaN;
p_real = real(p_cmplx);
p_phase = angle(p_cmplx); % phase of pressure

% Plot
h = figure;
ax = axes;
hsup = suptitle(sprintf("Pressure field from a point source, f = %d Hz, scatterer radius = %.2f m",frequencies(1),r_scat));
set(hsup,'FontSize',suptitle_font_size);
subplot(1,2,1,ax)
surf(X,Y,p_real);
ax.FontSize = plot_font_size;
bar = colorbar;
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
axis([-1 1 -1 1 -100 100]);
view(2)
%caxis([-0.1,0.1])
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])
title(sprintf("Amplitude"),'FontSize',title_font_size);
xlabel("m"); ylabel("m");
xticks(-1:0.5:1); yticks(-1:0.5:1);

ax = axes;
subplot(1,2,2,ax)
surf(X,Y,p_phase*180/pi);
ax.FontSize = plot_font_size;
bar = colorbar;
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -180 180]);
axis([-1 1 -1 1 -180 180]);
caxis([-180,180])
view(2)
line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])
title(sprintf("Phase"),'FontSize',title_font_size);
xlabel("m"); ylabel("m"); ylabel(bar, "Degrees");
xticks(-1:0.5:1); yticks(-1:0.5:1);

% Save image and reset default figure properties
print(h, image_path + filename, '-dsvg')
set(groot,'defaultFigurePosition','remove');
set(groot,'defaultFigurePaperType','remove');
set(groot,'defaultFigurePaperOrientation','remove');


%% Point source image - velocity

f = 1000;
c = 343;
w = 2*pi*f;
k = w/c;
t = 0;

% playback radius for computations, in meters
D = 1;
% grid resolution, in meters
d = 0.01;
x = -D:d:D;
y = -D:d:D;
[X,Y] = meshgrid(x,y);
% spherical coords
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X)));

image_path = "Thesis/images/";
frequencies = ones(1,1)*f;
source_positions = [[pi/2 0 2]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
filename = sprintf("velocity_f%d",frequencies(1));
N = ceil(exp(1)*k*(D)/2); % maximum order of expansion

[v_cmplx_original, ns_vector] = point_sources_velocity(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D);
v_real_original = real(v_cmplx_original);
v_mag_original = sqrt(sum(v_real_original.^2,3));
v_mag_original(R>D) = Inf;

% Plot
h = figure;
ax = axes;
hsup = suptitle(sprintf("Velocity field from a point source, f = %d Hz",frequencies(1)));
set(hsup,'FontSize',suptitle_font_size);
set(gcf,'Position',[2000 1000 1500 700])

subplot(1,2,1,ax)
surf(X,Y,v_mag_original);
ax.FontSize = plot_font_size;
colorbar;
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
axis([-1 1 -1 1 -100 100]);
%caxis([0,3e-4])
view(2)
%line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
%line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Magnitude of velocity"),'FontSize',title_font_size);
xlabel("m"); ylabel("m");
xticks(-1:0.5:1); yticks(-1:0.5:1);

ax = axes;
subplot(1,2,2,ax)
vx_real = v_real_original(:,:,1); vx_real(R>D) = 0;
vy_real = v_real_original(:,:,2); vy_real(R>D) = 0;
vz_real = v_real_original(:,:,3); vz_real(R>D) = 0;
hop_size = 16;
scaling = 1;
quiver(X(1:hop_size:end), Y(1:hop_size:end), vx_real(1:hop_size:end), vy_real(1:hop_size:end), scaling)
%quiver(ax,X,Y, vx_real, vy_real,2);
ax.FontSize = plot_font_size;
grid on
%quiver(X(1:2:end,1:2:end),Y(1:2:end,1:2:end),vx_real(1:2:end,1:2:end),vy_real(1:2:end,1:2:end))
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',2); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
axis([-1 1 -1 1 -100 100]);
view(2)
%line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
%line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Velocity"),'FontSize',title_font_size);
xlabel("m"); ylabel("m");
xticks(-1:0.5:1); yticks(-1:0.5:1);

% Save image and reset default figure properties
set(h,'Position',[2000 1000 1500 700]);
set(h,'PaperType','a4');
set(h,'PaperOrientation','portrait');
print(h, image_path + filename, '-dsvg')


%% Point source image with scatterer - velocity

f = 1000;
r_scat = 0.2;
c = 343;
w = 2*pi*f;
k = w/c;
t = 0;

% playback radius for computations, in meters
D = 1;
% grid resolution, in meters
d = 0.01;
x = -D:d:D;
y = -D:d:D;
[X,Y] = meshgrid(x,y);
% spherical coords
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X)));

image_path = "Thesis/images/";
frequencies = ones(1,1)*f;
source_positions = [[pi/2 0 2]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
filename = sprintf("velocity2_f%d_r%d",frequencies(1),r_scat*100);
N = ceil(exp(1)*k*(D)/2); % maximum order of expansion

[v_cmplx, ns_vector] = point_sources_scatterer_velocity(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D,r_scat);
v_real = real(v_cmplx);
v_mag = sqrt(sum(v_real.^2,3));
v_mag(R>D) = Inf;
v_mag(R<r_scat) = Inf;

% Plot
h = figure;
ax = axes;
hsup = suptitle(sprintf("Velocity field from a point source, f = %d Hz, scatterer radius = %.2f m",frequencies(1),r_scat));
set(hsup,'FontSize',suptitle_font_size);
set(gcf,'Position',[2000 1000 1500 700])

subplot(1,2,1,ax)
surf(X,Y,v_mag);
ax.FontSize = plot_font_size;
colorbar;
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
axis([-1 1 -1 1 -100 100]);
%caxis([0,3e-4])
view(2)
%line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
%line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Magnitude of velocity"),'FontSize',title_font_size);
xlabel("m"); ylabel("m");
xticks(-1:0.5:1); yticks(-1:0.5:1);
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])

ax = axes;
subplot(1,2,2,ax)
vx_real = v_real(:,:,1); vx_real(R>D) = 0; vx_real(R<r_scat) = 0;
vy_real = v_real(:,:,2); vy_real(R>D) = 0; vy_real(R<r_scat) = 0;
vz_real = v_real(:,:,3); vz_real(R>D) = 0; vz_real(R<r_scat) = 0;
hop_size = 16;
scaling = 1;
quiver(X(1:hop_size:end), Y(1:hop_size:end), vx_real(1:hop_size:end), vy_real(1:hop_size:end), scaling)
%quiver(ax,X,Y, vx_real, vy_real,2);
ax.FontSize = plot_font_size;
grid on
%quiver(X(1:2:end,1:2:end),Y(1:2:end,1:2:end),vx_real(1:2:end,1:2:end),vy_real(1:2:end,1:2:end))
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
axis([-1 1 -1 1 -100 100]);
view(2)
%line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
%line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Velocity"),'FontSize',title_font_size);
xlabel("m"); ylabel("m");
xticks(-1:0.5:1); yticks(-1:0.5:1);
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])


% Save image and reset default figure properties
set(h,'Position',[2000 1000 1500 700]);
set(h,'PaperType','a4');
set(h,'PaperOrientation','portrait');
print(h, image_path + filename, '-dsvg')


%% Point source image - intensity

f = 1000;
c = 343;
w = 2*pi*f;
k = w/c;
t = 0;

% playback radius for computations, in meters
D = 1;
% grid resolution, in meters
d = 0.01;
x = -D:d:D;
y = -D:d:D;
[X,Y] = meshgrid(x,y);
% spherical coords
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X)));

image_path = "Thesis/images/";
frequencies = ones(1,1)*f;
source_positions = [[pi/2 0 2]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
filename = sprintf("intensity_f%d",frequencies(1));
N = ceil(exp(1)*k*(D)/2); % maximum order of expansion

[v_cmplx_original, ns_vector] = point_sources_velocity(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D);
v_real_original = real(v_cmplx_original);
v_mag_original = sqrt(sum(v_real_original.^2,3));
v_mag_original(R>D) = Inf;
[p_cmplx, ns_vector] = point_sources_pressure(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D);

intensity = 0.5*real(p_cmplx.*conj(v_cmplx_original));
intensity_mag = sqrt(sum(intensity.^2,3));

% Plot
h = figure;
ax = axes;
hsup = suptitle(sprintf("Intensity field from a point source, f = %d Hz",frequencies(1)));
set(hsup,'FontSize',suptitle_font_size);
set(gcf,'Position',[100 100 1500 700])

subplot(1,2,1,ax)
surf(X,Y,intensity_mag);
ax.FontSize = plot_font_size;
colorbar;
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
axis([-1 1 -1 1 -100 100]);
%caxis([0,3e-4])
view(2)
%line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
%line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Magnitude of intensity"),'FontSize',title_font_size);
xlabel("m"); ylabel("m");
xticks(-1:0.5:1); yticks(-1:0.5:1);

ax = axes;
subplot(1,2,2,ax)
vx_real = intensity(:,:,1); vx_real(R>D) = 0;
vy_real = intensity(:,:,2); vy_real(R>D) = 0;
vz_real = intensity(:,:,3); vz_real(R>D) = 0;
hop_size = 16;
scaling = 2;
quiver(X(1:hop_size:end), Y(1:hop_size:end), vx_real(1:hop_size:end), vy_real(1:hop_size:end), scaling)
%quiver(ax,X,Y, vx_real, vy_real,2);
ax.FontSize = plot_font_size;
grid on
%quiver(X(1:2:end,1:2:end),Y(1:2:end,1:2:end),vx_real(1:2:end,1:2:end),vy_real(1:2:end,1:2:end))
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
axis([-1 1 -1 1 -100 100]);
view(2)
%line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
%line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Intensity"),'FontSize',title_font_size);
xlabel("m"); ylabel("m");
xticks(-1:0.5:1); yticks(-1:0.5:1);

% Save image and reset default figure properties
set(h,'Position',[2000 1000 1500 700]);
set(h,'PaperType','a4');
set(h,'PaperOrientation','portrait');
print(h, image_path + filename, '-dsvg')

%% Point source image with scatterer - intensity

f = 1000;
r_scat = 0.2;
c = 343;
w = 2*pi*f;
k = w/c;
t = 0;

% playback radius for computations, in meters
D = 1;
% grid resolution, in meters
d = 0.01;
x = -D:d:D;
y = -D:d:D;
[X,Y] = meshgrid(x,y);
% spherical coords
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X)));

image_path = "Thesis/images/";
frequencies = ones(1,1)*f;
source_positions = [[pi/2 0 2]];%[pi/4 0 1.2*D/2];[3*pi/4 0 1.2*D/2]];
magnitudes = [1];
initial_phases = [0];
filename = sprintf("intensity_f%d_r%d",frequencies(1),r_scat*100);
N = ceil(exp(1)*k*(D)/2); % maximum order of expansion

[v_cmplx_original, ns_vector] = point_sources_scatterer_velocity(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D,r_scat);
v_real_original = real(v_cmplx_original);
v_mag_original = sqrt(sum(v_real_original.^2,3));
v_mag_original(R>D) = Inf;
[p_cmplx, ns_vector] = point_sources_scatterer_pressure(t,frequencies,source_positions,magnitudes,initial_phases,X,Y,N,D,r_scat);
intensity = 0.5*real(p_cmplx.*conj(v_cmplx_original));
intensity_mag = sqrt(sum(intensity.^2,3));

% Plot
h = figure;
ax = axes;
hsup = suptitle(sprintf("Intensity field from a point source, f = %d Hz, scatterer radius = %.2f m",frequencies(1),r_scat));
set(hsup,'FontSize',suptitle_font_size);
set(gcf,'Position',[2000 1000 1500 700])

subplot(1,2,1,ax)
surf(X,Y,intensity_mag);
ax.FontSize = plot_font_size;
colorbar;
shading interp;
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
axis([-1 1 -1 1 -100 100]);
%caxis([0,3e-4])
view(2)
%line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
%line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Magnitude of intensity"),'FontSize',title_font_size);
xlabel("m"); ylabel("m");
xticks(-1:0.5:1); yticks(-1:0.5:1);
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])

ax = axes;
subplot(1,2,2,ax)
vx_real = intensity(:,:,1); vx_real(R>D) = 0;
vy_real = intensity(:,:,2); vy_real(R>D) = 0;
vz_real = intensity(:,:,3); vz_real(R>D) = 0;
hop_size = 16;
scaling = 2;
quiver(X(1:hop_size:end), Y(1:hop_size:end), vx_real(1:hop_size:end), vy_real(1:hop_size:end), scaling)
%quiver(ax,X,Y, vx_real, vy_real,2);
ax.FontSize = plot_font_size;
grid on
%quiver(X(1:2:end,1:2:end),Y(1:2:end,1:2:end),vx_real(1:2:end,1:2:end),vy_real(1:2:end,1:2:end))
% Draw sources
for i=1:size(ns_vector,1)
    line(ns_vector(i,1),ns_vector(i,2),2,'linestyle','none','marker','o','markerfacecolor','k','markersize',10); % draw source
end
axis equal
%axis([-1.2*D 1.2*D -1.2*D 1.2*D -100 100]);
axis([-1 1 -1 1 -100 100]);
view(2)
%line([0 D/5], [0 0], [2 2], 'Color','k','LineWidth', 1); % draw axes
%line([0 0], [0 D/5], [2 2], 'Color','k','LineWidth', 1);
title(sprintf("Intensity"),'FontSize',title_font_size);
xlabel("m"); ylabel("m");
xticks(-1:0.5:1); yticks(-1:0.5:1);
% draw scatterer
rectangle('Position',[-r_scat -r_scat 2*r_scat 2*r_scat],'Curvature',[1 1],'FaceColor',[0 0 0])


% Save image and reset default figure properties
set(h,'Position',[2000 1000 1500 700]);
set(h,'PaperType','a4');
set(h,'PaperOrientation','portrait');
print(h, image_path + filename, '-dsvg')


%% Plot loudspeaker setups

image_path = "Thesis/images/";
filename = sprintf("loudspeaker_setups");

% Helper function to convert directions from Matlab's azimuth-elevation to
% azimuth-inclination
aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

% T-design, 1st order, 4 loudspeakers
[tdesign_vecs_first_order, dirs_tdesign_ele_first_order] = getTdesign(2*1); % get the loudspeaker positions
dirs_tdesign_incl_first_order = aziElev2aziIncl(dirs_tdesign_ele_first_order); % convert to azi-incl
Ktdesign_first_order = size(dirs_tdesign_incl_first_order,1); % design size

% Dodeca, 3rd order, 20 loudspeakers
[u_dodeca, dirs_dodeca_ele, mesh_dodeca] = platonicSolid('dodeca');
dirs_dodeca_incl = aziElev2aziIncl(dirs_dodeca_ele); % convert to azi-incl
Kdodeca = size(dirs_dodeca_incl,1); % design size

% T-design, 4th order, 36 loudspeakers
[tdesign_vecs, dirs_tdesign_ele] = getTdesign(2*4); % get the loudspeaker positions
dirs_tdesign_incl = aziElev2aziIncl(dirs_tdesign_ele); % convert to azi-incl
Ktdesign = size(dirs_tdesign_incl,1); % design size

set(groot,'defaultFigurePosition',[2000 1000 1600 600]);
set(groot,'defaultFigurePaperType','a4');
set(groot,'defaultFigurePaperOrientation','portrait');

% Plot the loudspeaker setups
h = figure;
ax = axes;
hsup = suptitle(sprintf("Virtual loudspeaker setups"));
set(hsup,'FontSize',24);
           subplot(1,3,1,ax); plotSphFunctionTriangle(ones(length(dirs_tdesign_incl_first_order),1), dirs_tdesign_incl_first_order, 'real', ax); ax.FontSize = 20; view(-10,10); axis([-1 1 -1 1 -1 1]); title(sprintf('T-design, 1st order, %d loudspeakers',Ktdesign_first_order),'FontSize',20)
ax = axes; subplot(1,3,2,ax); plotSphFunctionTriangle(ones(length(dirs_dodeca_incl),1), dirs_dodeca_incl, 'real', ax); ax.FontSize = 20; view(-10,10); axis([-1 1 -1 1 -1 1]); title(sprintf('Dodeca, 3rd order, %d loudspeakers',Kdodeca),'FontSize',20)
ax = axes; subplot(1,3,3,ax); plotSphFunctionTriangle(ones(length(dirs_tdesign_incl),1), dirs_tdesign_incl, 'real', ax); ax.FontSize = 20; view(-10,10); axis([-1 1 -1 1 -1 1]); title(sprintf('T-design, 4th order, %d loudspeakers',Ktdesign),'FontSize',20)

% Save image and reset default figure properties
%print(h, image_path + filename, '-dsvg')
set(groot,'defaultFigurePosition','remove');
set(groot,'defaultFigurePaperType','remove');
set(groot,'defaultFigurePaperOrientation','remove');


%% Plane wave

f = 200;
c = 343;
rho0 = 1.225;
z0 = rho0 * c;
k = 2*pi*f/c;
x = linspace(0,10,100);
[X,Y]=meshgrid(x,x);
A = 1;
B = 1;
temp = linspace(0,1,16); t_vector = temp(1:length(temp)-1); % time vector

h = figure;
set(gcf,'Position',[100 100 2000 2000])
filename = 'test.gif';
for t = t_vector
    
    % temporal term
    temporal_term = exp(1i*2*pi*t);
    
    frequencies = ones(1,2)*200;
    amplitudes = ones(1,2)*1;
    directions = [30 -30];
    plane_waves = tilted_plane_waves(X,Y,frequencies,amplitudes,directions,true);
    p3c = plane_waves(:,:,1)*temporal_term;
    p4c = plane_waves(:,:,2)*temporal_term;
    
    % Velocities
    [u3c,v3c] = pol2cart(deg2rad(60),p3c/z0);
    [u4c,v4c] = pol2cart(deg2rad(directions(2)),p4c/z0);
    
    % draw
    quiver(real(u3c),real(v3c))
    title("1 plane wave propagating +x direction, particle velocity");
    axis([0 100 0 100])
    %set(gcf, 'Position', [2000, -200, 1000, 800])
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if t == 0 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append'); 
    end 
end


%% GIF encoder

filename = sprintf("Thesis/gifs/encoder_test.gif");

% Helper function to convert directions from Matlab's azimuth-elevation to
% azimuth-inclination
aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

% Eigenmike angles, in [azimuth, elevation] form
mic_dirs = ...
    [0    21;
    32     0;
     0   -21;
   328     0;
     0    58;
    45    35;
    69     0;
    45   -35;
     0   -58;
   315   -35;
   291     0;
   315    35;
    91    69;
    90    32;
    90   -31;
    89   -69;
   180    21;
   212     0;
   180   -21;
   148     0;
   180    58;
   225    35;
   249     0;
   225   -35;
   180   -58;
   135   -35;
   111     0;
   135    35;
   269    69;
   270    32;
   270   -32;
   271   -69];
mic_dirs_rad = mic_dirs*pi/180;
mic_dirs_incl = aziElev2aziIncl(mic_dirs_rad);
% flip the mic directions to correspond to the DOAs of each capsule
mic_dirs_incl_inv = mic_dirs_incl - pi;
nMics = size(mic_dirs,1);

% Variables
c = 343;
fs = 40000;
Lfilt = 800;
R = 0.042; % radius of EigenMike
nBins = Lfilt/2 + 1;
f_max = 20000;
kR_max = 2*pi*f_max*R/c;
array_order = ceil(2*kR_max);
f = (0:Lfilt/2)'*fs/Lfilt;
kR = 2*pi*f*R/c;

% Type and order of approximation
arrayType = 'rigid';
N_order = 30;

% Spherical harmonics for the mic array
sht_order = floor(sqrt(nMics)-1); % approximate for uniformly arranged mics
%Y_mics = sqrt(4*pi) * getSH(sht_order, mic_dirs_incl, 'real'); % real SH matrix for microphones
Y_mics = getSH(sht_order, mic_dirs_incl, 'real'); % real SH matrix for microphones

% Responses from DOAs
doas_deg = [0 0];
doas_rad = doas_deg*pi/180;
doas_incl = aziElev2aziIncl(doas_rad);
dops_rad = doas_rad - pi;
%dops_incl = aziElev2aziIncl(dops_rad);
dops_incl = doas_incl - pi;
Ndoas = size(doas_deg,1);
[~, H_array_sim] = simulateSphArray(Lfilt, mic_dirs_rad, doas_rad, arrayType, R, array_order, fs);

% Apply a plain SHT on the microphone responses without equalization.
M_mic2sh_sht = (1/nMics)*Y_mics;

% Apply single channel regularized inversion
maxG_dB = 15; % maximum allowed amplification
[H_filt_radinv, ~] = arraySHTfiltersTheory_radInverse(R, nMics, sht_order, Lfilt, fs, maxG_dB);
[H_filt_softLim, ~] = arraySHTfiltersTheory_softLim(R, nMics, sht_order, Lfilt, fs, maxG_dB);

% combine the per-order filters with the SHT matrix for evaluation of full filter matrix
pseudo_Y_mics = pinv(Y_mics);
a_nm_radinv = zeros(size(Y_mics,2),size(H_array_sim,1),size(H_array_sim,3));
a_nm_softLim = zeros(size(Y_mics,2),size(H_array_sim,1),size(H_array_sim,3));
for f_indx = 1:nBins
    %pseudo_Y_mics = pinv(squeeze(B(:,:,f_indx)));
    %a_nm_radinv(:,f_indx,:) = diag(replicatePerOrder(H_filt(f_indx,:),2))*M_mic2sh_sht'*H_array_sim(f_indx,:,:)';
    a_nm_radinv(:,f_indx,:) = diag(replicatePerOrder(H_filt_radinv(f_indx,:),2))*pseudo_Y_mics*H_array_sim(f_indx,:,:).';
    a_nm_softLim(:,f_indx,:) = diag(replicatePerOrder(H_filt_softLim(f_indx,:),2))*pseudo_Y_mics*H_array_sim(f_indx,:,:).';
    %a_nm_radinv(:,f_indx,:) = diag(replicatePerOrder(1./b_N(f_indx,:),2))*pseudo_Y_mics*H_array_sim(f_indx,:,:)';
    %a_nm_radinv(:,f_indx,:) = diag(replicatePerOrder(W_inv(f_indx,:),2))*pseudo_Y_mics*H_array_sim(f_indx,:,:)';
    %a_nm_radinv(:,f_indx,:) = squeeze(M_mic2sh_radinv(:,:,f_indx))*squeeze(H_array_sim(f_indx,:,:))';
end

% Flip the directional SH functions (The microphone response gets flipped for some
% unknown reason, i.e. a plane wave coming from +x is perceived to come
% from -x, when calculating the SH weights in the encoding part. As a
% quick-fix, we just flip the SH functions here to get the responses from 
% the right directions.)
%a_nm_radinv(2:end,:,:) = -real(a_nm_radinv(2:end,:,:));
%a_nm_softLim(2:end,:,:) = -real(a_nm_softLim(2:end,:,:));

% Reproduce the captured sound field in some virtual loudspeaker array
% Loudspeaker setup
N_tdesign = 4;
N_speakers = N_tdesign;
[~, speaker_dirs] = getTdesign(2*N_tdesign);
source_distance = 2;
speaker_dirs_dist = cat(2,speaker_dirs,ones(size(speaker_dirs,1),1)*source_distance);

% Decoder
decoding_method = 'sad';
[D_speakers, N_tdesign] = ambiDecoder(speaker_dirs*180/pi, decoding_method, 1, N_tdesign);

% Frequency vector
f_increment = fs/Lfilt;
f_vec = (0:Lfilt/2)*fs/Lfilt;
%f_selected = [f_increment f_increment f_increment 5*f_increment 10*f_increment 50*f_increment];
f_selected = [10*f_increment];
[~,f_indexes_selected] = ismembertol(f_selected,f_vec);

% Grid
D = 0.5;
d = 0.01;
x = -D:d:D;
y = -D:d:D;
[X,Y] = meshgrid(x,y);

% Pressure matrix
M = zeros(size(X,1),size(X,2),length(f_selected),Ndoas);
M_ideal = zeros(size(X,1),size(X,2),length(f_selected),Ndoas);
M_reference = zeros(size(X,1),size(X,2),length(f_selected),Ndoas);

% Gains for decoding only
%doa_ele = [-pi 0];
%doa_incl = aziElev2aziIncl(doa_ele); % azimuth/inclination
Y_ideal = getSH(N_tdesign, doas_incl, 'real');
%Y_ideal = getSH(N_tdesign, dops_incl, 'real');
gains_ideal = D_speakers*Y_ideal';

% Calculate the pressure fields
% Loudspeaker gains
%gains = D_speakers*squeeze(a_nm(:,:,doa_indx));
%gains = D_speakers*squeeze(a_nm_radinv(:,:,doa_indx));
gains = D_speakers*squeeze(a_nm_softLim(:,:,doa_indx));
%gains = D_speakers*squeeze(a_nm_test(:,:,doa_indx));
gains_abs = abs(gains);
gains_real = real(gains);
gains_phase = angle(gains);

tic
h = figure;
%set(gcf,'Position',[1500 100 2500 800])
set(gcf,'Position',[2000 100 300 800])
for t_div = 0:1/8:1-1/8
    f_indx = 1;
    for f = f_selected
        fprintf('Frequency %d/%d. ',f,f_selected(end))
        toc
        t = 0 + (1/f)*t_div;
        w = 2*pi*f; % angular frequency
        k = w/c; % wavenumber
        m0 = 1; % magnitude
        phi0 = [mod((source_distance*f/c),1)*2*pi]; % initial phase
        N_max = ceil(exp(1)*k*(D)/2); % maximum order of expansion

        pressure = zeros(size(X,1));
        pressure_ideal = zeros(size(X,1));

        % Calculate the pressure field
        for i=1:max(size(speaker_dirs))
            source_position = speaker_dirs_dist(i,:);
            
            % Calculations for all the points in the grid
            %[temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains_real(i,f_indexes_selected(f_indx)),gains_phase(i,f_indexes_selected(f_indx)),X,Y,N_max,D);
            %[temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains_real(i,f_indexes_selected(f_indx)),0,X,Y,N_max,D);
            [temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains(i,f_indexes_selected(f_indx)),phi0,X,Y,N_speakers,D);
            %[temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains_abs(i,f_indexes_selected(f_indx)),gains_phase(i,f_indexes_selected(f_indx)),X,Y,N_max,D);
            %[temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains_abs(i,f_indexes_selected(f_indx)),0,X,Y,N_max,D);
            pressure = pressure + temp_pressure;
            
            % Decoding only
            [temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains_ideal(i),phi0,X,Y,N_speakers,D);
            pressure_ideal = pressure_ideal + temp_pressure;
        end
        
        % Reference field
        pressure_reference = plane_wave_expansion_pressure(t,f,1,0,dops_incl,X,Y,N_max,D);
        
        % Store the pressure field
        M(:,:,f_indx,doa_indx) = pressure;
        M_ideal(:,:,f_indx,doa_indx) = pressure_ideal;
        M_reference(:,:,f_indx,doa_indx) = pressure_reference;
        f_indx = f_indx + 1;
    end
    
    % Plot
    subplot(3,1,1)
    surf(X,Y,((real(M(:,:,1,1)))))
    view(2)
    shading interp
    caxis([-0.7 0.7])
    subplot(3,1,2)
    surf(X,Y,((real(M_ideal(:,:,1,1)))))
    view(2)
    shading interp
    caxis([-0.05 0.05])
    subplot(3,1,3)
    surf(X,Y,((real(M_reference(:,:,1,1)))))
    view(2)
    shading interp
    caxis([-1 1])
    
    %{
    subplot(3,6,1)
    surf(X,Y,((real(M(:,:,1,1)))))
    view(2)
    shading interp
    caxis([-4 4])
    subplot(3,6,2)
    surf(X,Y,((real(M(:,:,2,1)))))
    view(2)
    shading interp
    caxis([-4 4])
    subplot(3,6,3)
    surf(X,Y,((real(M(:,:,3,1)))))
    view(2)
    shading interp
    caxis([-4 4])
    subplot(3,6,4)
    surf(X,Y,((real(M(:,:,4,1)))))
    view(2)
    shading interp
    caxis([-4 4])
    subplot(3,6,5)
    surf(X,Y,((real(M(:,:,5,1)))))
    view(2)
    shading interp
    caxis([-4 4])
    subplot(3,6,6)
    surf(X,Y,((real(M(:,:,6,1)))))
    view(2)
    shading interp
    caxis([-4 4])
    subplot(3,6,7)
    surf(X,Y,((real(M_ideal(:,:,1,1)))))
    view(2)
    shading interp
    caxis([-0.05 0.05])
    subplot(3,6,8)
    surf(X,Y,((real(M_ideal(:,:,2,1)))))
    view(2)
    shading interp
    caxis([-0.05 0.05])
    subplot(3,6,9)
    surf(X,Y,((real(M_ideal(:,:,3,1)))))
    view(2)
    shading interp
    caxis([-0.05 0.05])
    subplot(3,6,10)
    surf(X,Y,((real(M_ideal(:,:,4,1)))))
    view(2)
    shading interp
    caxis([-0.05 0.05])
    subplot(3,6,11)
    surf(X,Y,((real(M_ideal(:,:,5,1)))))
    view(2)
    shading interp
    caxis([-0.05 0.05])
    subplot(3,6,12)
    surf(X,Y,((real(M_ideal(:,:,6,1)))))
    view(2)
    shading interp
    caxis([-0.05 0.05])
    subplot(3,6,13)
    surf(X,Y,((real(M_reference(:,:,1,1)))))
    view(2)
    shading interp
    caxis([-1 1])
    subplot(3,6,14)
    surf(X,Y,((real(M_reference(:,:,2,1)))))
    view(2)
    shading interp
    caxis([-1 1])
    subplot(3,6,15)
    surf(X,Y,((real(M_reference(:,:,3,1)))))
    view(2)
    shading interp
    caxis([-1 1])
    subplot(3,6,16)
    surf(X,Y,((real(M_reference(:,:,4,1)))))
    view(2)
    shading interp
    caxis([-1 1])
    subplot(3,6,17)
    surf(X,Y,((real(M_reference(:,:,5,1)))))
    view(2)
    shading interp
    caxis([-1 1])
    subplot(3,6,18)
    surf(X,Y,((real(M_reference(:,:,6,1)))))
    view(2)
    shading interp
    caxis([-1 1])
    %}
    
    % GIF
    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if t == 0 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
    end 
end



%% Encoding gif 2

workfolder = 'Thesis';
file = 'Thesis/saved_files/test.mat';
% Grid
D = 0.5;
d = 0.01;
x = -D:d:D;
y = -D:d:D;
[X,Y] = meshgrid(x,y);
f_vec_test = [500];
doa_vec = [pi pi/2; pi/2 pi/2; pi/2 pi/2];
doa_vec_elev = [pi 0; pi/4 0];
doa_vec_elev = doa_vec_elev';
%dop_vec = [0 0; -pi/2 0];
decoding_method = 'sad'; %  'SAD','MMD','EPAD','ALLRAD','CSAD'
%[~, speaker_dirs] = getTdesign(2*4);
%[~, speaker_dirs, ~] = platonicSolid('dodeca');
dtu_setup = load('Thesis/saved_files/DTU_ls_dirs_deg.mat');
dtu_setup = dtu_setup.ls_dirs_deg;
speaker_dirs = dtu_setup*pi/180;
source_distance = 2;
indexes = [];
max_distance = 0;
N_diffuse_reps = 0;
delay = 0;% + (1/500)*1/4;
N_speakers = 0;
head_diameter = 0.2;
head_pos = [0 0 0];
encoding_regularization = 'radinv_dtu_moa_2d_sh';
setting = 'dtu_eigenmike_encoding';

temp = linspace(0,1,16); t_vector = temp(1:length(temp)-1); % time vector

h = figure;
filename = sprintf("Thesis/gifs/dtumic2D_sh2D_encoding.gif");
for t = t_vector
    
    delay = 0 + (1/500)*t;
    savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec_test,doa_vec_elev,decoding_method,speaker_dirs,source_distance,indexes,max_distance,N_diffuse_reps,delay,N_speakers,encoding_regularization)

    % Load the file
    M = load(file);
    M = M.M;
    
    % Plot
    subplot(1,2,1)
    surf(X,Y,((real(M(:,:,1)))))
    caxis([-0.4 0.4])
    %quiver(X,Y, real(M(:,:,1,1)), real(M(:,:,1,2)))
    view(2)
    shading interp
    %caxis([-4 4])
    subplot(1,2,2)
    surf(X,Y,((real(M(:,:,2)))))
    caxis([-0.4 0.4])
    %quiver(X,Y, real(M(:,:,2,1)), real(M(:,:,2,2)))
    view(2)
    shading interp
    %caxis([-0.05 0.05])
    
    drawnow 
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if t == 0 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1, 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','DelayTime',0.1,'WriteMode','append'); 
    end 
end




