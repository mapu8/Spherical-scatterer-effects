%% OBSELETE: Use savePressureFieldsTriton.m instead.
% Some useful functions at the end of this file, after the "Calculate ..."
% sections.

%%
% Calculate the directivity patterns for two points near the ear. Includes
% also calculations for diffuse field analysis. Save the calculated
% matrixes into files for later analysis in directivity_patterns.m.

% NOTE: Here inclination angle is the angle from the +z-axis,
%       elevation angle is the angle from xy-plane towards z-axis

clear;clc;

% Variables
head_diameter = 0.2;
head_pos = [pi 0 0.1]; % head position from the center [azi ele r]
head_direction = pi/2;  % azimuth pointing direction
r_scat = head_diameter/2;
D = head_diameter/2 + head_pos(3); % playback radius for calculations, in meters
d = 0.01; % grid resolution, in meters
file_p_ref = 'Thesis/saved_files/reference_pressure50_10k_1degree.mat';
file_p_dodeca = 'Thesis/saved_files/reproduced_pressure50_10k_1degree_dodeca.mat';
file_p_tdesign = 'Thesis/saved_files/reproduced_pressure50_10k_1degree_tdesign.mat';
file_p_ref_small = 'Thesis/saved_files/reference_pressure50_10k_small.mat';
file_p_dodeca_small = 'Thesis/saved_files/reproduced_pressure50_10k_dodeca_small.mat';
file_p_tdesign_small = 'Thesis/saved_files/reproduced_pressure50_10k_tdesign_small.mat';
file_v_ref = 'Thesis/saved_files/reference_velocity50_10k.mat';
file_v_dodeca = 'Thesis/saved_files/reproduced_velocity50_10k_dodeca.mat';
file_v_tdesign = 'Thesis/saved_files/reproduced_velocity50_10k_tdesign.mat';
file_v_ref_small = 'Thesis/saved_files/reference_velocity50_10k_small.mat';
file_v_dodeca_small = 'Thesis/saved_files/reproduced_velocity50_10k_dodeca_small.mat';
file_v_tdesign_small = 'Thesis/saved_files/reproduced_velocity50_10k_tdesign_small.mat';
file_diffuse_ref = 'Thesis/saved_files/diffuse_ref.mat';
file_diffuse_dodeca = 'Thesis/saved_files/diffuse2_dodeca.mat';
file_diffuse_tdesign = 'Thesis/saved_files/diffuse2_tdesign.mat';
file_diffuse_tdesign_first_order = 'Thesis/saved_files/diffuse2_tdesign_first_order.mat';
file_vbap_dodeca = 'Thesis/saved_files/vbap_dodeca.mat';
file_vbap_dodeca_noscat = 'Thesis/saved_files/vbap_dodeca_noscat.mat';
file_vbap_tdesign = 'Thesis/saved_files/vbap_tdesign.mat';
file_vbap_tdesign_noscat = 'Thesis/saved_files/vbap_tdesign_noscat.mat';
file_vbap_2D = 'Thesis/saved_files/vbap_2D.mat';
file_vbap_2D_noscat = 'Thesis/saved_files/vbap_2D_noscat.mat';
file_test = 'Thesis/saved_files/first_order_tdesign.mat';
f_vec = 50:50:10000; % frequencies
dop_azi_vec = (0:1:360)*pi/180; % azimuth angles for directions of propagation
doa_azi_vec = dop_azi_vec - pi;
dop_vec = dop_azi_vec;
dop_vec(2,:) = 0;

% Acoustics
c = 343;
rho_0 = 1.2;
Z_0 = c*rho_0;

if (head_pos(3) > D)
   error('Head/spherical scatterer is outside the reproduction field.'); 
end

% Grid
x = -D:d:D;
y = -D:d:D;
z = -D:d:D;
[X,Y] = meshgrid(x,y);

% Center point
center_x = ceil(size(X,1)/2);
center_y = ceil(size(Y,1)/2);

% Shift from the head position
[head_x,head_y,head_z] = sph2cart(head_pos(1),head_pos(2),head_pos(3));
head_x_amnt = round(-head_x/d);
head_y_amnt = round(-head_y/d);

% Single points picked from the field, offset from the center in [x y]
% in meters, first point should be origin (0 0)
offset_no_scat = [0 0; 0 0.02; 0 0.05; 0 0.1; 0.02 0; 0.05 0; 0.1 0; -0.2 0];
offset(:,1) = offset_no_scat(:,1) + r_scat + head_pos(3)*cos(head_pos(1));
offset(:,2) = offset_no_scat(:,2) + head_pos(3)*sin(head_pos(1));
offset_amnt_no_scat = round(offset_no_scat/d);
offset_amnt = round(offset/d);
% Points
points_no_scat = [(offset_amnt_no_scat(:,1) + center_x) (offset_amnt_no_scat(:,2) + center_y)];
points = [(offset_amnt(:,1) + center_x + head_x_amnt) (offset_amnt(:,2) + center_y + head_y_amnt)];
% Distance from the center
dist_p_no_scat = sqrt((offset_no_scat(:,1)).^2+(offset_no_scat(:,2)).^2);
max_dist_p_no_scat = max(dist_p_no_scat);
dist_p = sqrt((offset(:,1) + head_x).^2+(offset(:,2) + head_y).^2);
max_dist_p = max(dist_p);

% Indexes for single poinst in the field
indexes_no_scat = sub2ind(size(X), points_no_scat(:,2), points_no_scat(:,1));
indexes = sub2ind(size(X), points(:,2), points(:,1));

% Table of lambda values
lambda_title = {'f','wave_length','lambda/2','modulo(wave_length/0.02,1)','modulo(wave_length/0.05,1)','modulo(wave_length/0.1,1)','(lambda/2)/0.02','(lambda/2)/0.05','(lambda/2)/0.1'};
lambda_matrix = [f_vec' c./f_vec' ((c./f_vec)/2)' mod((c./f_vec')/0.02,1) mod((c./f_vec')/0.05,1) mod((c./f_vec')/0.1,1) (((c./f_vec)/2)/0.02)' (((c./f_vec)/2)/0.05)' (((c./f_vec)/2)/0.1)'];
lambda_table = [lambda_title;num2cell(lambda_matrix)];


%% Calculate the reference pressure matrix and store it to a file
% NOTE: no need to run, if the matrix is already stored

%savePressureFieldsReference(file_p_ref,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec);
savePressureFieldsReference(file_p_ref,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,indexes,max_dist_p);


%% Calculate the reproduced pressure matrix and store it to a file
% NOTE: no need to run, if the matrix is already stored

%[~, speaker_dirs, ~] = platonicSolid('dodeca');
N_tdesign = 2;
[~, speaker_dirs] = getTdesign(2*N_tdesign);

source_distance = 2;
decoding_method = 'sad';
%savePressureFieldsReproduction(file_p_dodeca,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,decoding_method,speaker_dirs,source_distance)
%savePressureFieldsReproduction(file_p_tdesign,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,decoding_method,speaker_dirs,source_distance)
savePressureFieldsReproduction(file_test,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p)

%N_tdesign = 4;
%[~, speaker_dirs] = getTdesign(2*N_tdesign);
%savePressureFieldsReproduction(file_p_tdesign,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p)

%% Calculate the reference velocity matrix and store it to a file
% NOTE: no need to run, if the matrix is already stored

%saveVelocityFieldsReference(file_v_ref,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec);
saveVelocityFieldsReference(file_test,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,test_indx);


%% Calculate the reproduced velocity matrix and store it to a file
% NOTE: no need to run, if the matrix is already stored

[~, speaker_dirs, ~] = platonicSolid('dodeca');
%N_tdesign = 4;
%[~, speaker_dirs] = getTdesign(2*N_tdesign);
source_distance = 2;
decoding_method = 'sad';
%saveVelocityFieldsReproduction(file_v_dodeca,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,decoding_method,speaker_dirs,source_distance)
%saveVelocityFieldsReproduction(file_v_tdesign,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,decoding_method,speaker_dirs,source_distance)
saveVelocityFieldsReproduction(file_test,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,decoding_method,speaker_dirs,source_distance,test_indx)


%% Calculate the reference diffuse field response and store it to a file
% NOTE: no need to run, if the matrix is already stored

[~, pw_dirs] = getTdesign(21);
pw_dirs = pw_dirs';
N_diffuse_reps = 1000;

savePressureFieldsReference(file_diffuse_ref,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,indexes,max_dist_p,N_diffuse_reps);


%% Calculate the reproduced diffuse field response and store it to a file
% NOTE: no need to run, if the matrix is already stored

[~, speaker_dirs, ~] = platonicSolid('dodeca');
%N_tdesign = 4;
%[~, speaker_dirs] = getTdesign(2*N_tdesign);

[~, pw_dirs] = getTdesign(21);
pw_dirs = pw_dirs';
N_diffuse_reps = 1000;
source_distance = 2;
decoding_method = 'sad';

savePressureFieldsReproduction(file_diffuse_dodeca,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,N_diffuse_reps);
%savePressureFieldsReproduction(file_diffuse_tdesign,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,N_diffuse_reps);


%% Calculate the reproduced pressure fields using only VBAP panning
% NOTE: no need to run, if the matrix is already stored

source_distance = 2;
[~, speaker_dirs_dodeca, ~] = platonicSolid('dodeca');
%savePressureFieldsVBAP(file_vbap_dodeca,head_diameter,head_pos,X,Y,D,d,50:50:200,dop_vec,speaker_dirs_dodeca,source_distance,indexes(1),max_dist_p)
savePressureFieldsVBAPNoScat(file_vbap_dodeca_noscat,X,Y,D,d,50:50:200,dop_vec,speaker_dirs_dodeca,source_distance,indexes_no_scat(1),max_dist_p_no_scat)

%N_tdesign = 4;
%[~, speaker_dirs_tdesign] = getTdesign(2*N_tdesign);
%savePressureFieldsVBAP(file_vbap_tdesign,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,speaker_dirs_tdesign,source_distance,indexes,max_dist_p)
%savePressureFieldsVBAP(file_vbap_tdesign_noscat,X,Y,D,d,f_vec,dop_vec,speaker_dirs_tdesign,source_distance,indexes_no_scat(1),max_dist_p_no_scat)

%speaker_dirs_2D(:,1) = 0:(2*pi/16):2*pi-2*pi/16;
%speaker_dirs_2D(:,2) = 0;
%savePressureFieldsVBAP(file_vbap_2D,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,speaker_dirs_2D,source_distance,indexes,max_dist_p)
%savePressureFieldsVBAPNoScat(file_vbap_2D_noscat,X,Y,D,d,f_vec,dop_vec,speaker_dirs_2D,source_distance,indexes_no_scat,max_dist_p_no_scat)


%% Combine Triton results
% Combine multiple files into one

% Output files
file_diffuse_ref = 'Thesis/saved_files_all/diffuse_reference_center.mat';
file_diffuse_dodeca = 'Thesis/saved_files_all/diffuse_dodeca_center.mat';
file_diffuse_tdesign = 'Thesis/saved_files_all/diffuse_tdesign_center.mat';
file_diffuse_tdesign_first_order = 'Thesis/saved_files_all/diffuse_tdesign_first_order_center.mat';
file_diffuse_listening_room_allrad = 'Thesis/saved_files_all/diffuse_listening_room_allrad2.mat';
file_diffuse_listening_room_sad = 'Thesis/saved_files_all/diffuse_listening_room_sad3.mat';

N_start_ref = 0;
N_end_ref = 49;
N_start_dodeca = 200;
N_end_dodeca = 299;
N_start_tdesign = 300;
N_end_tdesign = 399;
N_start_tdesign_first_order = 200;
N_end_tdesign_first_order = 299;
N_start_listening_room_allrad = 400;
N_end_listening_room_allrad = 499;
N_start_listening_room_sad = 600;
N_end_listening_room_sad = 699;
ref_prefix = 'Thesis/saved_files_all/diffuse_center/diffuse_reference_center';
dodeca_prefix = 'Thesis/saved_files_all/diffuse_center/diffuse_dodeca_center';
tdesign_prefix = 'Thesis/saved_files_all/diffuse_center/diffuse_tdesign_center';
tdesign_first_order_prefix = 'Thesis/saved_files_all/diffuse_center/diffuse_tdesign_first_order_center';
listening_room_allrad_prefix = 'Thesis/saved_files/diffuse_listening_room_allrad';
listening_room_sad_prefix = 'Thesis/saved_files/diffuse_listening_room_sad2';
combineResults(ref_prefix,N_start_ref,N_end_ref,file_diffuse_ref);
combineResults(dodeca_prefix,N_start_dodeca,N_end_dodeca,file_diffuse_dodeca);
combineResults(tdesign_prefix,N_start_tdesign,N_end_tdesign,file_diffuse_tdesign);
combineResults(tdesign_first_order_prefix,N_start_tdesign_first_order,N_end_tdesign_first_order,file_diffuse_tdesign_first_order);
%combineResults(listening_room_allrad_prefix,N_start_listening_room_allrad,N_end_listening_room_allrad,file_diffuse_listening_room_allrad);
%combineResults(listening_room_sad_prefix,N_start_listening_room_sad,N_end_listening_room_sad,file_diffuse_listening_room_sad);


%% Trim the saved files, so they can fit to github
% NOTE: no need to run, if the matrix is already stored

% Load the saved fields
M_p_ref = load(file_p_ref);
M_p_ref = M_p_ref.M;
M_p_dodeca = load(file_p_dodeca);
M_p_dodeca = M_p_dodeca.M;
M_p_tdesign = load(file_p_tdesign);
M_p_tdesign = M_p_tdesign.M;
%M_v_ref = load(file_v_ref);
%M_v_ref = M_v_ref.M;
%M_v_dodeca = load(file_v_dodeca);
%M_v_dodeca = M_v_dodeca.M;
%M_v_tdesign = load(file_v_tdesign);
%M_v_tdesign = M_v_tdesign.M;

% Center point and new matrix size
center_x = ceil(size(M_p_ref,2)/2);
center_y = ceil(size(M_p_ref,1)/2);
cut_size = 10; % offset from center point

% Cut the matrices
M_p_ref = M_p_ref([center_y-cut_size:center_y+cut_size],[center_x-cut_size:center_x+cut_size],:,:);
M_p_dodeca = M_p_dodeca([center_y-cut_size:center_y+cut_size],[center_x-cut_size:center_x+cut_size],:,:);
M_p_tdesign = M_p_tdesign([center_y-cut_size:center_y+cut_size],[center_x-cut_size:center_x+cut_size],:,:);
%M_v_ref = M_v_ref([center_y-cut_size:center_y+cut_size],[center_x-cut_size:center_x+cut_size],:,:,:);
%M_v_dodeca = M_v_dodeca([center_y-cut_size:center_y+cut_size],[center_x-cut_size:center_x+cut_size],:,:,:);
%M_v_tdesign = M_v_tdesign([center_y-cut_size:center_y+cut_size],[center_x-cut_size:center_x+cut_size],:,:,:);

% Save to a new file
save(file_p_ref_small, 'M_p_ref');
save(file_p_dodeca_small, 'M_p_dodeca');
save(file_p_tdesign_small, 'M_p_tdesign');
%save(file_v_ref_small, 'M_v_ref');
%save(file_v_dodeca_small, 'M_v_dodeca');
%save(file_v_tdesign_small, 'M_v_tdesign');


%% Select single points from the pressure field and store them into a file

% File paths
file_p_ref_full = 'Thesis/saved_files/non_diffuse_reference_0.mat';
file_p_dodeca_full = 'Thesis/saved_files/non_diffuse_dodeca_1.mat';
file_p_tdesign_full = 'Thesis/saved_files/non_diffuse_tdesign_2.mat';
file_p_ref_selected = 'Thesis/saved_files/reference_pressure50_10k_selected_points.mat';
file_p_dodeca_selected = 'Thesis/saved_files/dodeca_pressure50_10k_selected_points.mat';
file_p_tdesign_selected = 'Thesis/saved_files/tdesign_pressure50_10k_selected_points.mat';

% Single points picked from the field, offset from the center in [x y]
% in meters, first point should be origin (0 0)
offset = [0 0; 0 0.02; 0 0.05; 0 0.1; 0.02 0; 0.05 0; 0.1 0; -0.2 0];
offset_amnt = round(offset/d);
N_point_distances = size(offset,1); % number of different different distances between two points
% Points
points = [(offset_amnt(:,1) + center_x) (offset_amnt(:,2) + center_y)];
% Distance from the center
dist_p = sqrt((offset(:,1)).^2+(offset(:,2)).^2);
max_dist_p = max(dist_p);
% Distance between points
dist_p2_p1 = sqrt((offset(1,1)-offset(:,1)).^2+(offset(1,2)-offset(:,2)).^2);

% Indexes for single points in the field
indexes = sub2ind(size(X), points(:,2), points(:,1));

% Load the full matrices
M_p_ref = load(file_p_ref_full);
M_p_ref = M_p_ref.M;
M_p_dodeca = load(file_p_dodeca_full);
M_p_dodeca = M_p_dodeca.M;
M_p_tdesign = load(file_p_tdesign_full);
M_p_tdesign = M_p_tdesign.M;

% New matrices
M_p_ref_new = zeros(size(indexes,1),1,size(M_p_ref,3),size(M_p_ref,4));
M_p_dodeca_new = zeros(size(indexes,1),1,size(M_p_ref,3),size(M_p_ref,4));
M_p_tdesign_new = zeros(size(indexes,1),1,size(M_p_ref,3),size(M_p_ref,4));

% Select the points from the matrices
for f_index = 1:size(M_p_ref,4)
    for doa_index = 1:size(M_p_ref,3)
        temp_ref = M_p_ref(:,:,doa_index,f_index);
        temp_dodeca = M_p_dodeca(:,:,doa_index,f_index);
        temp_tdesign = M_p_tdesign(:,:,doa_index,f_index);
        M_p_ref_new(:,:,doa_index,f_index) = temp_ref(indexes);
        M_p_dodeca_new(:,:,doa_index,f_index) = temp_dodeca(indexes);
        M_p_tdesign_new(:,:,doa_index,f_index) = temp_tdesign(indexes);
    end
end

% Save the matrices
M = M_p_ref_new;
save(file_p_ref_selected, 'M');
M = M_p_dodeca_new;
save(file_p_dodeca_selected, 'M');
M = M_p_tdesign_new;
save(file_p_tdesign_selected, 'M');


%% Sum diffuse
% Sum plane waves used in diffuse field calculations (unnecessary, if the
% summation was already done when calculating the diffuse field)

% Output prefix
output_prefix = 'Thesis/saved_files/diffuse_dodeca_vbap_d02_summed';

N_start = 100;
N_end = 199;
prefix = 'Thesis/saved_files/diffuse_dodeca_vbap_d02';


for i = N_start:N_end
    filename = sprintf("%s_%d.mat",prefix,i);
    outputfile = sprintf("%s_%d.mat",output_prefix,i);
    M_temp = load(filename);
    M_temp = M_temp.M;
    M = sum(M_temp,3);
    save(outputfile, 'M');
end













