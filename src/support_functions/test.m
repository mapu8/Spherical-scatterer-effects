
% Helper function to convert directions from Matlab's azimuth-elevation to
% azimuth-inclination
aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

%% Test initial phase and gain boost

% Seem to work fine. Phase differences between plane wave and point source
% implementation with shorter distances from the loudspeaker (large
% distance --> point source acts as a plane wave).

f = 4000;
c = 343;
m0_test = 1;
r = 2;
t_test = 0;
w = 2*pi*f; % angular frequency
k = w/c; % wavenumber
source_position_test = [0 0 r];
% Grid
D_test = 0.2;
d_test = 0.01;
x_test = -D_test:d_test:D_test;
y_test = -D_test:d_test:D_test;
z_test = -D_test:d_test:D_test;
[X_test,Y_test] = meshgrid(x_test,y_test);

% Indexes
% Center point
center_x_test = ceil(size(X_test,1)/2);
center_y_test = ceil(size(Y_test,1)/2);

% Indexes
% Shift from the head position
head_diameter_test = 0.2;
head_pos_test = [pi 0 0.1]; % head position from the center [azi ele r]
head_direction_test = pi/2;  % azimuth pointing direction
r_scat_test = head_diameter_test/2;
[head_x_test,head_y_test,head_z_test] = sph2cart(head_pos_test(1),head_pos_test(2),head_pos_test(3));
head_x_amnt_test = round(-head_x_test/d);
head_y_amnt_test = round(-head_y_test/d);
offset_test = [0 0];
offset_amnt_test = round(offset_test/d);
points_test = [(offset_amnt_test(:,1) + center_x_test + head_x_amnt_test) (offset_amnt_test(:,2) + center_y_test + head_y_amnt_test)];
% Distance from the center (when the head is located in the center)
dist_p_test = sqrt((offset_test(:,1) + head_x_test).^2+(offset_test(:,2) + head_y_test).^2);
max_dist_p_test = max(dist_p_test);
% Distance between points
dist_p2_p1_test = sqrt((offset_test(1,1)-offset_test(:,1)).^2+(offset_test(1,2)-offset_test(:,2)).^2);
% Indexes for single points in the field
indexes_test = sub2ind(size(X_test), points_test(:,2), points_test(:,1));
indexes_center = sub2ind(size(X_test), center_y_test, center_x_test);
N_max_test = ceil(exp(1)*k*(max_dist_p_test)/2); % maximum order of expansion

r_new = r + head_pos_test(3);
f_vec_test = 50:50:10000;
vec_test_repr_noscat = [];
i = 0;
for f = f_vec_test
    i = i + 1;
    w = 2*pi*f; % angular frequency
    k = w/c; % wavenumber
    source_position_test_new = [0 0 r_new];

    % Without scatterer
    phi0_test_repr_noscat = mod((r*f/c),1)*2*pi;
    gain_boost_test_noscat = m0_test / real(m0_test*exp(1i*mod((r*f/c),1)*2*pi)*exp(1i*(w*t_test*1 - k*sqrt(r^2)))/(4*pi*sqrt(r^2)));
    phi0_test_ref_noscat = 0; % initial phase
    
    % With scatterer
    phi0_test_repr_scat = mod((r_new*f/c),1)*2*pi;
    gain_boost_test_scat = m0_test / real(m0_test*exp(1i*mod((r_new*f/c),1)*2*pi)*exp(1i*(w*t_test*1 - k*sqrt(r_new^2)))/(4*pi*sqrt(r_new^2)));
    dop_ele_test = [source_position_test(1)-pi source_position_test(2)];
    dop_incl_test = [dop_ele_test(1) pi/2-dop_ele_test(2)];
    [dop_x,dop_y,dop_z] = sph2cart(dop_ele_test(1),dop_ele_test(2),1);
    dop_cart_test = [dop_x dop_y dop_z];
    head_cart_test = [head_x_test head_y_test head_z_test];
    wave_distance_test = dot(head_cart_test,dop_cart_test)/norm(dop_cart_test);
    phi0_test_ref_scat = -mod((wave_distance_test*f/c),1)*2*pi; % initial phase
    
    
    [temp_pressure_repr_noscat, ~] = point_sources_pressure(t_test,f,source_position_test,gain_boost_test_noscat,phi0_test_repr_noscat,X_test,Y_test,N_max_test,D_test,indexes_center);
    [temp_pressure_ref_noscat] = plane_wave_expansion_pressure(t_test,f,m0_test,phi0_test_ref_noscat,dop_incl_test,X_test,Y_test,N_max_test,D_test,indexes_center);

    [temp_pressure_repr, ~] = point_sources_scatterer_pressure(t_test,f,source_position_test_new,gain_boost_test_scat,phi0_test_repr_scat,X_test,Y_test,N_max_test,D_test,r_scat_test,indexes_test);
    [temp_pressure_ref] = plane_wave_scatterer_pressure(t_test,f,m0_test,phi0_test_ref_scat,dop_incl_test,X_test,Y_test,N_max_test,D_test,r_scat_test,indexes_test);

    % No scatterer point source
    vec_test_repr_noscat(i,1) = temp_pressure_repr_noscat;
    vec_test_repr_noscat(i,2) = abs(temp_pressure_repr_noscat);
    
    % No scatterer plane wave
    vec_test_repr_noscat(i,3) = temp_pressure_ref_noscat;
    vec_test_repr_noscat(i,4) = abs(temp_pressure_ref_noscat);
    
    % Scatterer point source
    vec_test_repr_noscat(i,5) = temp_pressure_repr;
    vec_test_repr_noscat(i,6) = abs(temp_pressure_repr);

    % Scatterer plane wave
    vec_test_repr_noscat(i,7) = temp_pressure_ref;
    vec_test_repr_noscat(i,8) = abs(temp_pressure_ref);
end


%% Test Triton script

workfolder = "Thesis/";
outputfile = "Thesis/saved_files/test.mat";
%delay_distance=0.02;
N_diffuse_repetitions = 20;
setting = "dodeca";
rng_seed = 1;
savePressureFieldsTriton(workfolder,outputfile,N_diffuse_repetitions,setting,rng_seed);%,delay_distance);

% Load the output matrix
M_test = load(outputfile);
M_test = M_test.M;
test1 = squeeze(real(M_test(:,:,4,1)));

% Grid
D = 1.1;
d = 0.05;
x = -D:d:D;
y = -D:d:D;
z = -D:d:D;
[X,Y] = meshgrid(x,y);
% Spherical coordinates
[Azi,Elev,R] = cart2sph(X,Y,zeros(size(X)));

% Test plot
figure
surf(X,Y,test1);
view(2)
shading interp


%% Test VBAP

f_vec_test = 50:50:200;
dop_azi_vec = (0:1:360)*pi/180;
file_vbap = 'Thesis/saved_files/vbap_dodeca.mat';
M_vbap = load(file_vbap);
M_vbap = M_vbap.M;
pressure = squeeze((abs(M_vbap(1,:,:,:,:))));
pressure = circshift(pressure,ceil(size(pressure,1)/2),1);


figure
surf(f_vec_test,dop_azi_vec*180/pi,pressure);
view(2);
shading interp;
xlim([0 10000]);
ylim([0 360]);
xticks([0 100 200 500 1000 2000 3000 5000 10000]);
yticks([0 90 180 270 360]);
ylabel("DOA (degrees)");
xlabel("f [Hz]");
set(gca,'Xscale','log')


%% Test velocity function
t_test = 0;
frequencies_test = [1000];
c = 343;
w_test = 2*pi*frequencies_test(1); % angular frequency
k_test = w_test/c; % wavenumber
magnitudes_test = [1];
initial_phases_test = [0];
source_positions_test = [0 0 0.1];
% Grid
D_test = 1.1;
d_test = 0.05;
x_test = -D_test:d_test:D_test;
y_test = -D_test:d_test:D_test;
z_test = -D_test:d_test:D_test;
[X_test,Y_test] = meshgrid(x_test,y_test);
N_test = ceil(exp(1)*k_test*(D_test)/2);
%[v_cmplx, position_vectors] = point_sources_scatterer_velocity(t_test,frequencies_test,source_positions_test,magnitudes_test,initial_phases_test,X_test,Y_test,N_test,D_test,r_scat_test);
head_diameter_test = 0.2;
head_pos_test = [pi 0 0.1]; % head position from the center [azi ele r]
head_direction_test = pi/2;  % azimuth pointing direction
r_scat_test = head_diameter_test/2;
f_vec_test = 50:50:10000;
dop_azi_vec_test = (0:5:360)*pi/180; % azimuth angles for directions of propagation
doa_azi_vec_test = dop_azi_vec_test - pi;
dop_vec_test = dop_azi_vec_test;
dop_vec_test(2,:) = 0;
[~, speaker_dirs_test, ~] = platonicSolid('dodeca');
source_distance_test = 2;
file_test = 'Thesis/saved_files/test.mat';
saveVelocityFieldsReproduction(file_test,head_diameter_test,head_pos_test,X_test,Y_test,D_test,d_test,f_vec_test,dop_vec_test,'sad',speaker_dirs_test,source_distance_test)

%% Test ambidecoder order
N_tdesign_test = 4;
% azi/ele radians
[~, speaker_dirs_test] = getTdesign(2*N_tdesign_test);
%[~, speaker_dirs_test, ~] = platonicSolid('dodeca');
source_distance_test = 2;
decoding_method_test = 'sad';
% azi/ele degrees
[D_test1, N_test] = ambiDecoder(speaker_dirs_test(:,[1 2])*180/pi, decoding_method_test,1)

% azi/ele degrees
getLayoutAmbisonicOrder(speaker_dirs_test*180/pi)


%% Random calculations

% Mixed-ordered Ambisonics system, from: 
% An Introduction to Higher Order Ambisonic by Florian Hollerweger, corrected version October 2008
Mh = 5; % horizontal order
Mv = 1:3; % vertical order
Nmixed = (2*Mh+1) + ((Mv+1).^2 - (2*Mv+1))


%% Test angles

angl_doa = [45 0];
angl_doa_rad = angl_doa*pi/180;
angl_doa_incl = aziElev2aziIncl(angl_doa_rad);
angl_dop_rad = angl_doa_rad - pi;
angl_dop_incl = angl_doa_incl - pi;
angl_dop_rad_proper = [angl_doa(1)-pi 0];
angl_dop_incl_proper = [angl_doa_incl(1)-pi angl_doa_incl(2)];
angl_dop_incl_zenith_flip = [angl_doa_incl(1) angl_doa_incl(2)-pi];
angl_dop_incl_from_elev = aziElev2aziIncl(angl_dop_rad);

angl_dop2 = [0 0];
angl_dop2_rad = angl_dop2*pi/180;
angl_dop2_incl = aziElev2aziIncl(angl_dop2_rad);
angl_doa2 = [180 0];
angl_doa2_rad = angl_doa2*pi/180;
angl_doa2_incl = aziElev2aziIncl(angl_doa2_rad);
angl_doa2_from_dop = angl_doa2_incl - pi;

[n_dop_x,n_dop_y,n_dop_z] = sph2cart(angl_dop2_rad(1),angl_dop2_rad(2),1);
[doas_rad(1), doas_rad(2), ~] = cart2sph(-n_dop_x,-n_dop_y,-n_dop_z);
doas_rad_incl = aziElev2aziIncl([doas_rad(1) doas_rad(2)]);

Y_real = getSH(1, angl_doa2_incl, 'real')
Y_from_dop = getSH(1, angl_doa2_from_dop, 'real')
Y_from_sph2cart = getSH(1, doas_rad_incl, 'real')




%%
n_dop_rad = round([cos(angl_dop_rad(1))*cos((angl_dop_rad(2))) sin(angl_dop_rad(1))*cos((angl_dop_rad(2))) sin((angl_dop_rad(2)))],4)

n_dop_incl_proper = round([cos(angl_dop_incl_proper(1))*sin(abs(angl_dop_incl_proper(2))) sin(angl_dop_incl_proper(1))*sin(abs(angl_dop_incl_proper(2))) cos(abs(angl_dop_incl_proper(2)))],4)
n_dop_incl_zenith_flip = round([cos(angl_dop_incl_zenith_flip(1))*sin(abs(angl_dop_incl_zenith_flip(2))) sin(angl_dop_incl_zenith_flip(1))*sin(abs(angl_dop_incl_zenith_flip(2))) cos(abs(angl_dop_incl_zenith_flip(2)))],4)
n_dop_incl = round([cos(angl_dop_incl(1))*sin(abs(angl_dop_incl(2))) sin(angl_dop_incl(1))*sin(abs(angl_dop_incl(2))) cos(abs(angl_dop_incl(2)))],4)
n_dop_incl_noabs = round([cos(angl_dop_incl(1))*sin((angl_dop_incl(2))) sin(angl_dop_incl(1))*sin((angl_dop_incl(2))) cos((angl_dop_incl(2)))],4)
n_dop_incl_from_elev = round([cos(angl_dop_incl_from_elev(1))*sin(abs(angl_dop_incl_from_elev(2))) sin(angl_dop_incl_from_elev(1))*sin(abs(angl_dop_incl_from_elev(2))) cos(abs(angl_dop_incl_from_elev(2)))],4)
n_doa_incl = round([cos(angl_doa_incl(1))*sin(abs(angl_doa_incl(2))) sin(angl_doa_incl(1))*sin(abs(angl_doa_incl(2))) cos(abs(angl_doa_incl(2)))],4)


try
    close test_fig
catch e
end
test_fig = figure;

plot3([0 n_doa_incl(1)],[0 n_doa_incl(2)],[0 n_doa_incl(3)],'DisplayName','DOA')
hold on
%plot3([0 n_dop_incl_proper(1)],[0 n_dop_incl_proper(2)],[0 n_dop_incl_proper(3)],'DisplayName','proper')
plot3([0 n_dop_incl_zenith_flip(1)],[0 n_dop_incl_zenith_flip(2)],[0 n_dop_incl_zenith_flip(3)],'DisplayName','zenith flip')
plot3([0 n_dop_incl(1)],[0 n_dop_incl(2)],[0 n_dop_incl(3)],'DisplayName','incl')
plot3([0 n_dop_incl_noabs(1)],[0 n_dop_incl_noabs(2)],[0 n_dop_incl_noabs(3)],'DisplayName','incl no abs')
%plot3([0 n_dop_incl_from_elev(1)],[0 n_dop_incl_from_elev(2)],[0 n_dop_incl_from_elev(3)],'DisplayName','from elev')
grid on
legend


%% Test 4pi*source_distance normalization + fixed doa calculation

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
decoding_method = 'allrad_normal'; %  'SAD','MMD','EPAD','ALLRAD','CSAD'
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
N_speakers = 5;
head_diameter = 0.2;
head_pos = [0 0 0];
encoding_regularization = 'radinv_dtu_moa';
setting = 'dtu_eigenmike_encoding';

%savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec_test,doa_vec_elev,decoding_method,speaker_dirs,source_distance)
savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec_test,doa_vec_elev,decoding_method,speaker_dirs,source_distance,indexes,max_distance,N_diffuse_reps,delay,N_speakers)
%savePressureFieldsReferenceNoScat(file,X,Y,D,d,f_vec_test,doa_vec_elev,indexes,max_distance,N_diffuse_reps,delay)
%savePressureFieldsReference(file,head_diameter,head_pos,X,Y,D,d,f_vec_test,doa_vec_elev,indexes,max_distance,N_diffuse_reps,delay)
%savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec_test,doa_vec_elev,decoding_method,speaker_dirs,source_distance,indexes,max_distance,N_diffuse_reps,delay,N_speakers)
%savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec_test,doa_vec_elev,speaker_dirs,source_distance,N_speakers,indexes,max_distance,N_diffuse_reps,delay)
%savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec_test,doa_vec_elev,speaker_dirs,source_distance,N_speakers,indexes,max_distance,delay)

%saveVelocityFieldsReference(file,head_diameter,head_pos,X,Y,D,d,f_vec_test,doa_vec_elev,indexes,max_distance,delay)
%saveVelocityFieldsReferenceNoScat(file,X,Y,D,d,f_vec_test,doa_vec_elev,indexes,max_distance,delay)
%saveVelocityFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec_test,doa_vec_elev,decoding_method,speaker_dirs,source_distance,indexes,max_distance,delay,N_speakers)
%saveVelocityFieldsReproductionNoScat(file,X,Y,D,d,f_vec_test,doa_vec_elev,decoding_method,speaker_dirs,source_distance,indexes,max_distance,delay,N_speakers)
%saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec_test,doa_vec_elev,speaker_dirs,source_distance,N_speakers,indexes,max_distance,delay)
%saveVelocityFieldsVBAPNoScat(file,X,Y,D,d,f_vec_test,doa_vec_elev,speaker_dirs,source_distance,N_speakers,indexes,max_distance,delay)

%savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec_test,doa_vec_elev,decoding_method,speaker_dirs,source_distance,indexes,max_distance,N_diffuse_reps,delay,N_speakers,encoding_regularization)
%savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec_test,doa_vec_elev,decoding_method,speaker_dirs,source_distance,indexes,max_distance,N_diffuse_reps,delay,N_speakers,encoding_regularization)

%savePressureFieldsTriton(workfolder,file,N_diffuse_reps,setting,1,0,encoding_regularization)

% Load the file
M = load(file);
M = M.M;

max1 = real(M(ceil(size(X,1)/2),ceil(size(X,2)/2),1));
max2 = real(M(ceil(size(X,1)/2),ceil(size(X,2)/2),2));
text = sprintf("Center1: %.3f\nCenter2: %.3f\n",round(max1,4),round(max2,4));
disp(text)

% Plot 
figure

% Plot
subplot(1,2,1)
surf(X,Y,((real(M(:,:,1)))))
%quiver(X,Y, real(M(:,:,1,1)), real(M(:,:,1,2)))
view(2)
shading interp
%caxis([-4 4])
subplot(1,2,2)
surf(X,Y,((real(M(:,:,2)))))
%quiver(X,Y, real(M(:,:,2,1)), real(M(:,:,2,2)))
view(2)
shading interp
%caxis([-0.05 0.05])


%% DTU setup

dtu_setup = load('Thesis/saved_files/DTU_ls_dirs_deg.mat');
dtu_setup = dtu_setup.ls_dirs_deg;
dtu_rad = dtu_setup*pi/180;

getLayoutAmbisonicOrder(dtu_setup)

temp = round(dtu_setup(:,1));
temp(temp>180) = temp(temp>180)-360;
ls_dirs_deg = cat(2,temp,round(dtu_setup(:,2)));

%save('Thesis/saved_files/DTU_ls_dirs_deg_-180_180','ls_dirs_deg');

%%
test_anm = load('Thesis/saved_files/temp_anm.mat');
test_anm = test_anm.temp_anm;

test_anm_ones = load('Thesis/saved_files/temp_anm_ones.mat');
test_anm_ones = test_anm_ones.temp_anm;

%% Test velocity

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
doa_vec_elev = [pi 0; pi/4 0; 0 0];
doa_vec_elev = doa_vec_elev';
%dop_vec = [0 0; -pi/2 0];
decoding_method = 'allrad'; %  'SAD','MMD','EPAD','ALLRAD','CSAD'
%[~, speaker_dirs] = getTdesign(2*4);
%[~, speaker_dirs, ~] = platonicSolid('dodeca');
dtu_setup = load('Thesis/saved_files/DTU_ls_dirs_deg.mat');
dtu_setup = dtu_setup.ls_dirs_deg;
speaker_dirs = dtu_setup*pi/180;
source_distance = 2;
indexes = [20 20 50 100]';
max_distance = 0;
N_diffuse_reps = 0;
delay = 0;% + (1/500)*1/4;
N_speakers = 0;
head_diameter = 0.2;
head_pos = [0 0 0];
encoding_regularization = 'radinv';
setting = 'ref_velocity';

saveVelocityFieldsReference(file,head_diameter,head_pos,X,Y,D,d,f_vec_test,doa_vec_elev,indexes,max_distance,delay)
%savePressureFieldsTriton(workfolder,file,N_diffuse_reps,setting,1,0,encoding_regularization)


%% Caclulate effective frequency range in off-axis location (sound source in the front, off-axis location to the side)

source_distance = 2;
off_axis = 0.1;
speaker_angl_diff = 22.5; % Angle between two speakers
speaker_angl_diff_rad = speaker_angl_diff*pi/180;
c = 343;

% Distances from the sources to the new point
r1 = sqrt((source_distance*cos(speaker_angl_diff_rad/2))^2 + (source_distance*sin(speaker_angl_diff_rad/2) - off_axis)^2);
r2 = sqrt((source_distance*cos(speaker_angl_diff_rad/2))^2 + (source_distance*sin(speaker_angl_diff_rad/2) + off_axis)^2);
dr = abs(r1-r2)

% Time difference
dt = dr/c;

% Freuquency when the time delay is half-wavelenght
f_dip = 1/(2*dt)


%% Save DTU mic directions

mic_dirs_azi_elev_rad = ...
    [0.0                   1.40194351305251;
    -3.14159265358979      1.40194351305251;
    0.523598775598299      0.964876234040517;
    1.57079632679490       0.964876234040517;
    2.61799387799149       0.964876234040517;
    -2.61799387799150      0.964876234040517;
    -1.57079632679490      0.964876234040517;
    0.523598775598299      0.964876234040517;
    0.628318530717959      0.502232586957315;
    1.25663706143592       0.502232586957315;
    1.88495559215388       0.502232586957315;
    2.51327412287183       0.502232586957315;
    -3.14159265358979      0.502232586957315;
    -2.51327412287183      0.502232586957315;
    -1.88495559215388      0.502232586957315;
    -1.25663706143592      0.502232586957315;
    -0.628318530717959     0.502232586957315;
    0.0                    0.502232586957315;
    0.589048622548086      0.0;
    0.981747704246811      0.0;
    1.37444678594553       0.0;
    1.76714586764426       0.0;
    2.15984494934298       0.0;
    2.55254403104171       0.0;
    2.94524311274043       0.0;
    -2.94524311274043      0.0;
    -2.55254403104171      0.0;
    -2.15984494934298      0.0;
    -1.76714586764426      0.0;
    -1.37444678594553      0.0;
    -0.981747704246811     0.0;
    -0.589048622548086     0.0;
    -0.196349540849362     0.0;
    0.196349540849362      0.0;
    0.314159265358979      -0.502232586957315;
    0.942477796076938      -0.502232586957315;
    1.57079632679490       -0.502232586957315;
    2.19911485751286       -0.502232586957315;
    2.82743338823081       -0.502232586957315;
    -2.82743338823081      -0.502232586957315;
    -2.19911485751286      -0.502232586957315;
    -1.57079632679490      -0.502232586957315;
    -0.942477796076938     -0.502232586957315;
    -0.314159265358979     -0.502232586957315;
    0.0                    -0.964876234040517;
    1.04719755119660       -0.964876234040517;
    2.09439510239320       -0.964876234040517;
    3.14159265358979       -0.964876234040517;
    -2.09439510239320      -0.964876234040517;
    -1.04719755119660      -0.964876234040517;
    -1.57079632679490      -1.40194351305251;
    1.57079632679490       -1.40194351305251];
save('Thesis/saved_files/DTU_mic_dirs.mat','mic_dirs_azi_elev_rad');


%% Test condition number k(B*B) = ||B*B||_2 ||(B*B)^-1||_2

fs = 48000;
Lfilt = 960;
maxG_dB = 15;
doa_vec_elev = [pi 0; pi/4 0; 0 0];
doa_vec_elev = doa_vec_elev';
encoding_regularization = 'radinv';

encodeWithEigenMike(doa_vec_elev,encoding_regularization,fs,Lfilt,maxG_dB);

























