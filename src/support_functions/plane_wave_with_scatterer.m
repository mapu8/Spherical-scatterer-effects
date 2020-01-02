%% Analyze the effect of a scatterer to the plane wave pressure, velocity and active intensity

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

%% Calculate the impulse responses from the sphere

% Microphone settings
R = 0.1;        % Radius of the sphere in the middle.
Rmax = 0.5;     % Radius of the measurement area.
Nmics = 40;     % Number of mics in an array extending from the center 
                % to the edge of the measurement area.
mic_azi = 0;        % Azimuth angle of the mics (rad)
mic_elev = 0;       % Elevation angle of the mics (rad)
% Directions of arrival of the plane waves (rad)
doa_azi = (0:1:359)'*pi/180;
doa_dirs = [doa_azi zeros(size(doa_azi))];

% Meshgrid
x = linspace(0,1,100);
[X,Y]=meshgrid(x,x); 

% Plane waves
%{
Nwaves = 2;
c = 343;
rho0 = 1.225;
z0 = rho0 * c;
k = 2*pi*f/c;
A = 1;
frequencies = ones(1,Nwaves)*f;
amplitudes = ones(1,Nwaves)*A;
propagation_directions = [180 90];
DOAs = propagation_directions - 180;
plane_waves = tilted_plane_waves(X,Y,frequencies,amplitudes,propagation_directions,false);
%}

% Microphone settings
N_order = 40; % Order of expansion
% Multiple mics with increasing distance from the sphere
micDistances = linspace(R,Rmax,Nmics)';
mic_dirs_sph = [ones(Nmics,1)*mic_azi ones(Nmics,1)*mic_elev micDistances];

% Sample rate for the sphericalScatterer()
fs = 40000;
Lfilt = 2000;
Nfft = Lfilt;
frange = (0:Nfft/2)*fs/Nfft;
fgap = fs/(2*size(frange,2));  % Frequency gap between two matrix elements in 
                               % the H_mic_sph impulse responses returned
                               % by sphericalScatterer().

% Compute simulated impulse responses
[~, H_mic_sph] = sphericalScatterer(mic_dirs_sph, doa_dirs, R, N_order, Lfilt, fs);

%% Pick single frequencies from the impulse response

% Selected frequencies
frequencies = [200, 500, 2000, 3000];
individual_plots = [200, 2000];

% Transform the impulse responses from spherical to cartesian coordinates
Nfrequencies = max(size(frequencies));
freq_indexes = round(frequencies./fgap);
H_single_frequencies = H_mic_sph(freq_indexes,:,:);
for n=1:size(H_single_frequencies,1)
    M(n,:,:) = fromSPHtoCART(X, micDistances, squeeze(H_single_frequencies(n,:,:)), mic_azi, doa_azi);
end

% Determine the form of the subplots
if (Nfrequencies <= 3)
    subplotform = [1,3];
else
    subplotform = [2, round(Nfrequencies/2)]; 
end

% Plot
% Waveforms
figure(35)
for i=1:Nfrequencies
    subplot(subplotform(1),subplotform(2),i)
    surf(X,Y,real(squeeze(M(i,:,:))));
    view(2)
    shading flat
    caxis([-1 1]), colorbar, colormap(jet)
    title(sprintf("Plane wave field with a spherical scatterer, f=%d Hz", frequencies(i)))
end

% SPL
figure(36)
for i=1:Nfrequencies
    subplot(subplotform(1),subplotform(2),i)
    surf(X,Y,20*log10(abs(squeeze(M(i,:,:)))));
    view(2)
    shading flat
    caxis([-20 10]), colorbar, colormap(jet)
    title(sprintf("Plane wave field with a spherical scatterer, f=%d Hz, SPL", frequencies(i)))
end

% Individual plots
Nplots = 36;
for i=1:max(size(individual_plots))
    index = find(frequencies==individual_plots(i));
    if(isempty(index))
       sprintf("Couldn't find frequency %d from the calculated frequencies.", individual_plots(i))
       continue
    end
    Nplots = Nplots + 1;
    
    figure(Nplots)
    subplot(1,2,1)
    surf(X,Y,real(squeeze(M(index,:,:))));
    view(2)
    shading flat
    caxis([-1 1]), colorbar, colormap(jet)
    title(sprintf("Plane wave field with a spherical scatterer, f=%d Hz", frequencies(index)))
    
    subplot(1,2,2)
    surf(X,Y,20*log10(abs(squeeze(M(index,:,:)))));
    view(2)
    shading flat
    caxis([-20 10]), colorbar, colormap(jet)
    title(sprintf("Plane wave field with a spherical scatterer, f=%d Hz, SPL", frequencies(index)))
end


