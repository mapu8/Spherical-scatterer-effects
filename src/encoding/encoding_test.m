%% Test encoding
%
% Record plane waves with virtual EigenMike and use the recorded signals as
% input signals for the virtual loudspeaker arrays.

clear all; clc;

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

% Tetrahedron microphone array
%[~,mic_dirs_rad] = getTdesign(2);
%mic_dirs = mic_dirs_rad*180/pi;
%mic_dirs_incl = aziElev2aziIncl(mic_dirs_rad);

nMics = size(mic_dirs,1);

% Variables
c = 343;
fs = 48000;
Lfilt = 960;
R = 0.042; % radius of EigenMike
nBins = Lfilt/2 + 1;
f_max = 10000;
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
doas_deg = [0 0; 90 0; -180 0]; % [azi ele]
%doas_deg = mic_dirs(1,:);
doas_rad = doas_deg*pi/180;
doas_incl = aziElev2aziIncl(doas_rad);
%dops_rad = doas_rad - pi;
%dops_incl = aziElev2aziIncl(dops_rad);
%dops_incl = doas_incl - pi;

[n_doa_x,n_doa_y,n_doa_z] = sph2cart(doas_rad(:,1),doas_rad(:,2),1);
[dops_rad(:,1), dops_rad(:,2), ~] = cart2sph(-n_doa_x,-n_doa_y,-n_doa_z);
dops_incl = aziElev2aziIncl(dops_rad);

Ndoas = size(doas_deg,1);
[h_array_sim, H_array_sim_full] = simulateSphArray(Lfilt, mic_dirs_rad, doas_rad, arrayType, R, array_order, fs);
%H_array_sim = real(H_array_sim);

% Spherical harmonics for DOAs
Y_doas = sqrt(4*pi) * getSH(sht_order, doas_incl, 'real');

% Sound field coefficients on the surface of a rigid sphere
b_N = sphModalCoeffs(sht_order, kR, arrayType, []);

% Apply a plain SHT on the microphone responses without equalization.
M_mic2sh_sht = (1/nMics)*Y_mics;

% Pseudo-inversion of the microphone SH functions
pseudo_Y_mics = pinv(Y_mics);

% Apply single channel regularized inversion
maxG_dB = 15; % maximum allowed amplification
[H_filt_radinv_full, ~] = arraySHTfiltersTheory_radInverse(R, nMics, sht_order, Lfilt, fs, maxG_dB);
[H_filt_softLim_full, ~] = arraySHTfiltersTheory_softLim(R, nMics, sht_order, Lfilt, fs, maxG_dB);
[H_filt_regLS_full, ~] = arraySHTfiltersTheory_regLS(R, mic_dirs_rad, sht_order, Lfilt, fs, maxG_dB);

% Select frequencies up to f_max
H_array_sim = H_array_sim_full(1:floor(f_max*Lfilt/fs)+1,:,:);
H_filt_radinv = H_filt_radinv_full(1:floor(f_max*Lfilt/fs)+1,:);
H_filt_softLim = H_filt_softLim_full(1:floor(f_max*Lfilt/fs)+1,:);
H_filt_regLS = H_filt_regLS_full(:,:,1:floor(f_max*Lfilt/fs)+1);

% combine the per-order filters with the SHT matrix for evaluation of full filter matrix
a_nm_radinv = zeros(size(Y_mics,2),size(H_array_sim,1),size(H_array_sim,3));
a_nm_softLim = zeros(size(Y_mics,2),size(H_array_sim,1),size(H_array_sim,3));
a_nm_regLS = zeros(size(Y_mics,2),size(H_array_sim,1),size(H_array_sim,3));
for doa_indx = 1:Ndoas
    for f_indx = 1:floor(f_max*Lfilt/fs)+1
        %a_nm_radinv(:,f_indx,:) = diag(replicatePerOrder(H_filt(f_indx,:),2))*M_mic2sh_sht'*H_array_sim(f_indx,:,:)';
        a_nm_radinv(:,f_indx,doa_indx) = diag(replicatePerOrder(H_filt_radinv(f_indx,:),2))*pseudo_Y_mics*H_array_sim(f_indx,:,doa_indx).';
        a_nm_softLim(:,f_indx,doa_indx) = diag(replicatePerOrder(H_filt_softLim(f_indx,:),2))*pseudo_Y_mics*H_array_sim(f_indx,:,doa_indx).';
        a_nm_regLS(:,f_indx,doa_indx) = H_filt_regLS(:,:,f_indx)*H_array_sim(f_indx,:,doa_indx).';
        %a_nm_radinv(:,f_indx,:) = diag(replicatePerOrder(1./b_N(f_indx,:),2))*pseudo_Y_mics*H_array_sim(f_indx,:,:)';
        %a_nm_radinv(:,f_indx,:) = diag(replicatePerOrder(W_inv(f_indx,:),2))*pseudo_Y_mics*H_array_sim(f_indx,:,:)';
        %a_nm_radinv(:,f_indx,:) = squeeze(M_mic2sh_radinv(:,:,f_indx))*squeeze(H_array_sim(f_indx,:,:))';
    end
end

%% Reproduce the captured sound field in some virtual loudspeaker array

% Loudspeaker setup
N_tdesign = 4;
%a_nm_radinv = a_nm_radinv(1:16,:);
%a_nm_softLim = a_nm_softLim(1:16,:);

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
%f_selected = [f_increment f_increment 5*f_increment 10*f_increment 50*f_increment 100*f_increment];
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
Y_ideal = getSH(N_tdesign, doas_incl, 'real');
gains_ideal_all = D_speakers*Y_ideal';
% Compensate distance attenuation
gains_ideal_all = 4*pi*gains_ideal_all*source_distance;

% Calculate the pressure fields
for doa_indx = 1:Ndoas
    % Loudspeaker gains
    %gains = D_speakers*squeeze(a_nm(:,:,doa_indx));
    %gains = D_speakers*squeeze(a_nm_radinv(:,:,doa_indx));
    %gains = D_speakers*squeeze(a_nm_softLim(:,:,doa_indx));
    gains = D_speakers*squeeze(a_nm_regLS(:,:,doa_indx))*sqrt(4*pi);
    %gains = D_speakers*squeeze(a_nm_test(:,:,doa_indx));
    gains_abs = abs(gains);
    gains_real = real(gains);
    gains_phase = angle(gains);

    % Compensate the distance attenuation
    gains = gains*source_distance;
    gains_ideal = gains_ideal_all(:,doa_indx);
    
    tic
    f_indx = 1;
    for f = f_selected
        fprintf('Frequency %d/%d. ',f,f_selected(end))
        toc
        t = 0 + (1/f)/8;
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
            %[temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains_real(i,f_indexes_selected(f_indx)),gains_phase(i,f_indexes_selected(f_indx)),X,Y,N_speakers,D);
            %[temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains_real(i,f_indexes_selected(f_indx)),0,X,Y,N_speakers,D);
            [temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains(i,f_indexes_selected(f_indx)),phi0,X,Y,N_speakers,D);
            %[temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains_abs(i,f_indexes_selected(f_indx)),gains_phase(i,f_indexes_selected(f_indx)),X,Y,N_speakers,D);
            %[temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains_abs(i,f_indexes_selected(f_indx)),0,X,Y,N_speakers,D);
            pressure = pressure + temp_pressure;
            
            % Decoding only
            [temp_pressure, ~] = point_sources_pressure(t,f,source_position,gains_ideal(i),phi0,X,Y,N_speakers,D);
            %temp_pressure = plane_wave_expansion_pressure(t,f,gains_ideal(i),ph0,aziElev2aziIncl(speaker_dirs(i,:))-pi,X,Y,N_speakers,D);
            pressure_ideal = pressure_ideal + temp_pressure;
            
            % dop [azimuth inclination] in radians
            %plane_wave_expansion_pressure(t,f,m0,phi0,dop,X,Y,N,D,indexes)
        end
        
        % Reference field
        pressure_reference = plane_wave_expansion_pressure(t,f,1,0,doas_rad(doa_indx,:),X,Y,N_max,D);
        
        % Store the pressure field
        M(:,:,f_indx,doa_indx) = pressure;
        M_ideal(:,:,f_indx,doa_indx) = pressure_ideal;
        M_reference(:,:,f_indx,doa_indx) = pressure_reference;
        f_indx = f_indx + 1;
    end
end

% Spit out values for debugging
max_encoding = max(max(real(M)));
max_ideal = max(max(real(M_ideal)));
max_ref = max(max(real(M_reference)));
text = sprintf("Encoding: %.2f\nDecoding: %.2f\nReference: %.2f\n",round(max_encoding,3),round(max_ideal,3),round(max_ref,3));
disp(text)


%% Plot
f_indx = 2;
doa_indx = 3;
figure

% Plot
subplot(3,1,1)
surf(X,Y,((real(M(:,:,1,doa_indx)))))
view(2)
shading interp
%caxis([-4 4])
subplot(3,1,2)
surf(X,Y,((real(M_ideal(:,:,1,doa_indx)))))
view(2)
shading interp
%caxis([-0.05 0.05])
subplot(3,1,3)
surf(X,Y,((real(M_reference(:,:,1,doa_indx)))))
view(2)
shading interp
%caxis([-1 1])

%{
subplot(3,6,1)
surf(X,Y,((real(M(:,:,1,1)))))
view(2)
shading interp
subplot(3,6,2)
surf(X,Y,((real(M(:,:,2,1)))))
view(2)
shading interp
subplot(3,6,3)
surf(X,Y,((real(M(:,:,3,1)))))
view(2)
shading interp
subplot(3,6,4)
surf(X,Y,((real(M(:,:,4,1)))))
view(2)
shading interp
subplot(3,6,5)
surf(X,Y,((real(M(:,:,5,1)))))
view(2)
shading interp
subplot(3,6,6)
surf(X,Y,((real(M(:,:,6,1)))))
view(2)
shading interp
subplot(3,6,7)
surf(X,Y,((real(M_ideal(:,:,1,1)))))
view(2)
shading interp
subplot(3,6,8)
surf(X,Y,((real(M_ideal(:,:,2,1)))))
view(2)
shading interp
subplot(3,6,9)
surf(X,Y,((real(M_ideal(:,:,3,1)))))
view(2)
shading interp
subplot(3,6,10)
surf(X,Y,((real(M_ideal(:,:,4,1)))))
view(2)
shading interp
subplot(3,6,11)
surf(X,Y,((real(M_ideal(:,:,5,1)))))
view(2)
shading interp
subplot(3,6,12)
surf(X,Y,((real(M_ideal(:,:,6,1)))))
view(2)
shading interp
subplot(3,6,13)
surf(X,Y,((real(M_reference(:,:,1,1)))))
view(2)
shading interp
subplot(3,6,14)
surf(X,Y,((real(M_reference(:,:,2,1)))))
view(2)
shading interp
subplot(3,6,15)
surf(X,Y,((real(M_reference(:,:,3,1)))))
view(2)
shading interp
subplot(3,6,16)
surf(X,Y,((real(M_reference(:,:,4,1)))))
view(2)
shading interp
subplot(3,6,17)
surf(X,Y,((real(M_reference(:,:,5,1)))))
view(2)
shading interp
subplot(3,6,18)
surf(X,Y,((real(M_reference(:,:,6,1)))))
view(2)
shading interp
%}


