% Analyze and plot the pressure errors of the calculated pressures at mic positions. The
% pressure matrices are calculated with savePressureFieldsXXX.m functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 18/12/2019
% Lauros Pajunen, lauros.pajunen@alumni.aalto.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% This is a supplementary code for [1]. The plots in [1] can be reproduced
% by running this script.
%
%
% The simulated 5D matrices are constructed with the following dimensions:
% 1D - point in the field, 8 points with indexes:
%      1: center (right ear position when head in the field)
%      2: 2cm front
%      3: 5cm front
%      4: 10cm front
%      5: 2cm right side
%      6: 5cm  right side
%      7: 10cm right side
%      8: 10cm left side (left ear position when head in the field)
% 2D - just 1
% 3D - used DOP (361 angles accross the horizontal plane, for example)
% 4D - frequencies (200 frequencies, if 50:50:10000)
% 5D - loudspeaker setup/reproduction (reference, dodeca, t-design,...,etc.)
%
% In total the matrices are sized something like M(8,1,361,200,6).
% Example: M(1,1,91,50,2) --> outputs the calculated pressure values for
% the center point (next to right ear) from 91-1=90 DOP angle (DOA=270) in
% unit circle notation, with frequency of 2500 Hz, using dodeca
% reproduction setup. Check savePressureFieldsXXX.m for more detail on how
% the matrices are constructed.

% NOTE: Here inclination angle is the angle from the +z-axis,
%       elevation angle is the angle from xy-plane towards z-axis

clear;clc;

% Path to 'data' folder (default is that the code is run
% in Spherical-scatterer-effects folder, i.e. pwd = '../Spherical-scatterer-effects')
data_path = "data/";

% Frequency and DOA vectors used in the simulations
f_vec = 50:50:10000; % frequencies
dop_azi_vec = (0:1:360)*pi/180; % azimuth angles for directions of propagation
dop_vec = dop_azi_vec;
dop_vec(2,:) = 0;

% Plot properties
global polar_font_size plot_font_size suptitle_font_size title_font_size legend_font_size
polar_font_size = 22;
plot_font_size = 22;
suptitle_font_size = 30;
title_font_size = 22;
legend_font_size = 22;

% Paths to files containing the simulated data
% Reference
file_p_ref_center = data_path + 'simulated_responses/non_diffuse_reference_center_4.mat';
file_p_ref_noscat = data_path + 'simulated_responses/non_diffuse_reference_noscat_204.mat';

% Dodecahedron, 20 speakers, N=3
file_p_dodeca_center = data_path + 'simulated_responses/non_diffuse_dodeca5_center_203.mat';
file_p_dodeca_noscat = data_path + 'simulated_responses/non_diffuse_dodeca3_noscat_205.mat';
file_vbap_dodeca_center = data_path + 'simulated_responses/vbap_dodeca_center_216.mat';
file_vbap_dodeca_noscat = data_path + 'simulated_responses/vbap_dodeca_noscat_211.mat';
file_dodeca_encoding_radinv_center = data_path + 'simulated_responses/non_diffuse_dodeca_eigenmike_encoding2_radinv_center_201.mat';
file_dodeca_encoding_radinv_noscat = data_path + 'simulated_responses/non_diffuse_dodeca_eigenmike_encoding2_radinv_noscat_205.mat';

% T-design, 36 speakers, N=4
file_p_tdesign_center = data_path + 'simulated_responses/non_diffuse_tdesign4_center_204.mat';
file_p_tdesign_noscat = data_path + 'simulated_responses/non_diffuse_tdesign4_noscat_206.mat';
file_vbap_tdesign_center = data_path + 'simulated_responses/vbap_tdesign_center_217.mat';
file_vbap_tdesign_noscat = data_path + 'simulated_responses/vbap_tdesign_noscat_212.mat';
file_tdesign_encoding_radinv_center = data_path + 'simulated_responses/non_diffuse_tdesign_eigenmike_encoding2_radinv_center_301.mat';
file_tdesign_encoding_radinv_noscat = data_path + 'simulated_responses/non_diffuse_tdesign_eigenmike_encoding2_radinv_noscat_305.mat';

% T-design, 4 speakers (tetrahedron), N=1
file_p_tdesign_first_order_center = data_path + 'simulated_responses/non_diffuse_tdesign_first_order_center_7.mat';
file_p_tdesign_first_order_noscat = data_path + 'simulated_responses/non_diffuse_tdesign_first_order_noscat_207.mat';
file_vbap_tdesign_first_order_center = data_path + 'simulated_responses/vbap_tdesign_first_order_center_218.mat';
file_vbap_tdesign_first_order_noscat = data_path + 'simulated_responses/vbap_tdesign_first_order_noscat_215.mat';
file_tdesign_first_order_encoding_radinv_center = data_path + 'simulated_responses/non_diffuse_tdesign_first_order_eigenmike_encoding2_radinv_center_401.mat';
file_tdesign_first_order_encoding_radinv_noscat = data_path + 'simulated_responses/non_diffuse_tdesign_first_order_eigenmike_encoding2_radinv_noscat_405.mat';

% Audio Visual Immersion Lab setup in Technical University of Denmark (DTU), 64 speakers, N=5
file_p_dtu_center_sad = data_path + 'simulated_responses/non_diffuse_dtu_center_sad_508.mat';
file_p_dtu_noscat_sad = data_path + 'simulated_responses/non_diffuse_dtu_noscat_sad_509.mat';
file_p_dtu_center_allrad = data_path + 'simulated_responses/non_diffuse_dtu_center_allrad_5thorder_518.mat';
file_p_dtu_noscat_allrad = data_path + 'simulated_responses/non_diffuse_dtu_noscat_allrad_5thorder_519.mat';
file_p_dtu_center_allrad_normal = data_path + 'simulated_responses/non_diffuse_dtu_center_allrad_normal_5thorder_514.mat';
file_p_dtu_noscat_allrad_normal = data_path + 'simulated_responses/non_diffuse_dtu_noscat_allrad_normal_5thorder_515.mat';
file_p_dtu_center_mmd = data_path + 'simulated_responses/non_diffuse_dtu_center_mmd_5thorder_514.mat';
file_p_dtu_noscat_mmd = data_path + 'simulated_responses/non_diffuse_dtu_noscat_mmd_5thorder_515.mat';
file_vbap_dtu_center = data_path + 'simulated_responses/vbap_dtu_center_510.mat';
file_vbap_dtu_noscat = data_path + 'simulated_responses/vbap_dtu_noscat_511.mat';
file_dtu_encoding_radinv_center_allrad = data_path + 'simulated_responses/non_diffuse_dtu_eigenmike_encoding3_radinv_center_allrad_501.mat';
file_dtu_encoding_radinv_noscat_allrad = data_path + 'simulated_responses/non_diffuse_dtu_eigenmike_encoding3_radinv_noscat_allrad_505.mat';
file_dtu_encoding_radinv_center_sad = data_path + 'simulated_responses/non_diffuse_dtu_eigenmike_encoding2_radinv_center_sad_501.mat';
file_dtu_encoding_radinv_noscat_sad = data_path + 'simulated_responses/non_diffuse_dtu_eigenmike_encoding2_radinv_noscat_sad_505.mat';
file_dtu_encoding_dtumic2D_radinv_center_allrad = data_path + 'simulated_responses/non_diffuse_dtu_dtumicMOA_2D_encoding_radinv_center_allrad2_502.mat';
file_dtu_encoding_dtumic3D_radinv_center_allrad = data_path + 'simulated_responses/non_diffuse_dtu_dtumic3D_encoding_radinv_center_allrad_5thorder_503.mat';
file_dtu_encoding_dtumic3D_radinv_center_mmd = data_path + 'simulated_responses/non_diffuse_dtu_dtumic3D_encoding_radinv_center_mmd_503.mat';
file_dtu_encoding_dtumic3D_radinv_center_allrad_normal = data_path + 'simulated_responses/non_diffuse_dtu_dtumic3D_encoding_radinv_center_allrad_normal_501.mat';
file_dtu_encoding_dtumicMOA_radinv_center_allrad = data_path + 'simulated_responses/non_diffuse_dtu_dtumicMOA_encoding_radinv_center_allrad_no4pi_502.mat';
file_dtu_encoding_dtumic2D_radinv_noscat_allrad = data_path + 'simulated_responses/non_diffuse_dtu_dtumicMOA_2D_encoding_radinv_noscat_allrad2_506.mat';
file_dtu_encoding_dtumic3D_radinv_noscat_allrad = data_path + 'simulated_responses/non_diffuse_dtu_dtumic3D_encoding_radinv_noscat_allrad_5thorder_502.mat';
file_dtu_encoding_dtumic3D_radinv_noscat_mmd = data_path + 'simulated_responses/non_diffuse_dtu_dtumic3D_encoding_radinv_noscat_mmd_502.mat';
file_dtu_encoding_dtumic3D_radinv_noscat_allrad_normal = data_path + 'simulated_responses/non_diffuse_dtu_dtumic3D_encoding_radinv_noscat_allrad_normal_500.mat';
file_dtu_encoding_dtumicMOA_radinv_noscat_allrad = data_path + 'simulated_responses/non_diffuse_dtu_dtumicMOA_encoding_radinv_noscat_allrad_no4pi_506.mat';

% Elevated circle
file_vbap_dtu_center_vertical = data_path + 'simulated_responses/vbap_dtu_vertical180_510.mat';
file_vbap_dtu_noscat_vertical = data_path + 'simulated_responses/vbap_dtu_noscat_vertical180_511.mat';
file_p_ref_center_vertical = data_path + 'simulated_responses/non_diffuse_reference_center_vertical180_1.mat';
file_p_ref_noscat_vertical = data_path + 'simulated_responses/non_diffuse_reference_noscat_vertical180_2.mat';
file_p_dtu_center_allrad_vertical = data_path + 'simulated_responses/non_diffuse_dtu_center_allrad_5thorder_vertical180_518.mat';
file_p_dtu_noscat_allrad_vertical = data_path + 'simulated_responses/non_diffuse_dtu_noscat_allrad_5thorder_vertical180_519.mat';


%% Load the simulated responses

% Reference
M_p_ref_noscat = load(file_p_ref_noscat);
M_p_ref_noscat = M_p_ref_noscat.M;
M_p_ref_center = load(file_p_ref_center);
M_p_ref_center = M_p_ref_center.M;

% Dodeca
M_p_dodeca_noscat = load(file_p_dodeca_noscat);
M_p_dodeca_noscat = M_p_dodeca_noscat.M;
M_p_dodeca_center = load(file_p_dodeca_center);
M_p_dodeca_center = M_p_dodeca_center.M;
M_vbap_dodeca_center = load(file_vbap_dodeca_center);
M_vbap_dodeca_center = M_vbap_dodeca_center.M;
M_vbap_dodeca_noscat = load(file_vbap_dodeca_noscat);
M_vbap_dodeca_noscat = M_vbap_dodeca_noscat.M;
M_dodeca_encoding_radinv_center = load(file_dodeca_encoding_radinv_center);
M_dodeca_encoding_radinv_center = M_dodeca_encoding_radinv_center.M;
M_dodeca_encoding_radinv_noscat = load(file_dodeca_encoding_radinv_noscat);
M_dodeca_encoding_radinv_noscat = M_dodeca_encoding_radinv_noscat.M;

% T-design 4th order
M_p_tdesign_noscat = load(file_p_tdesign_noscat);
M_p_tdesign_noscat = M_p_tdesign_noscat.M;
M_p_tdesign_center = load(file_p_tdesign_center);
M_p_tdesign_center = M_p_tdesign_center.M;
M_vbap_tdesign_center = load(file_vbap_tdesign_center);
M_vbap_tdesign_center = M_vbap_tdesign_center.M;
M_vbap_tdesign_noscat = load(file_vbap_tdesign_noscat);
M_vbap_tdesign_noscat = M_vbap_tdesign_noscat.M;
M_tdesign_encoding_radinv_center = load(file_tdesign_encoding_radinv_center);
M_tdesign_encoding_radinv_center = M_tdesign_encoding_radinv_center.M;
M_tdesign_encoding_radinv_noscat = load(file_tdesign_encoding_radinv_noscat);
M_tdesign_encoding_radinv_noscat = M_tdesign_encoding_radinv_noscat.M;

% T-design 1st order
M_p_tdesign_first_order_noscat = load(file_p_tdesign_first_order_noscat);
M_p_tdesign_first_order_noscat = M_p_tdesign_first_order_noscat.M;
M_p_tdesign_first_order_center = load(file_p_tdesign_first_order_center);
M_p_tdesign_first_order_center = M_p_tdesign_first_order_center.M;
M_vbap_tdesign_first_order_center = load(file_vbap_tdesign_first_order_center);
M_vbap_tdesign_first_order_center = M_vbap_tdesign_first_order_center.M;
M_vbap_tdesign_first_order_noscat = load(file_vbap_tdesign_first_order_noscat);
M_vbap_tdesign_first_order_noscat = M_vbap_tdesign_first_order_noscat.M;
M_tdesign_first_order_encoding_radinv_center = load(file_tdesign_first_order_encoding_radinv_center);
M_tdesign_first_order_encoding_radinv_center = M_tdesign_first_order_encoding_radinv_center.M;
M_tdesign_first_order_encoding_radinv_noscat = load(file_tdesign_first_order_encoding_radinv_noscat);
M_tdesign_first_order_encoding_radinv_noscat = M_tdesign_first_order_encoding_radinv_noscat.M;

% DTU
M_p_dtu_center_sad = load(file_p_dtu_center_sad);
M_p_dtu_center_sad = M_p_dtu_center_sad.M;
M_p_dtu_noscat_sad = load(file_p_dtu_noscat_sad);
M_p_dtu_noscat_sad = M_p_dtu_noscat_sad.M;
M_p_dtu_center_allrad = load(file_p_dtu_center_allrad);
M_p_dtu_center_allrad = M_p_dtu_center_allrad.M;
M_p_dtu_noscat_allrad = load(file_p_dtu_noscat_allrad);
M_p_dtu_noscat_allrad = M_p_dtu_noscat_allrad.M;
M_p_dtu_center_allrad_normal = load(file_p_dtu_center_allrad_normal);
M_p_dtu_center_allrad_normal = M_p_dtu_center_allrad_normal.M;
M_p_dtu_noscat_allrad_normal = load(file_p_dtu_noscat_allrad_normal);
M_p_dtu_noscat_allrad_normal = M_p_dtu_noscat_allrad_normal.M;
M_p_dtu_center_mmd = load(file_p_dtu_center_mmd);
M_p_dtu_center_mmd = M_p_dtu_center_mmd.M;
M_p_dtu_noscat_mmd = load(file_p_dtu_noscat_mmd);
M_p_dtu_noscat_mmd = M_p_dtu_noscat_mmd.M;
M_vbap_dtu_center = load(file_vbap_dtu_center);
M_vbap_dtu_center = M_vbap_dtu_center.M;
M_vbap_dtu_noscat = load(file_vbap_dtu_noscat);
M_vbap_dtu_noscat = M_vbap_dtu_noscat.M;

% Elevated circle
M_p_reference_center_vertical = load(file_p_ref_center_vertical);
M_p_reference_center_vertical = M_p_reference_center_vertical.M;
M_p_reference_noscat_vertical = load(file_p_ref_noscat_vertical);
M_p_reference_noscat_vertical = M_p_reference_noscat_vertical.M;
M_vbap_dtu_center_vertical = load(file_vbap_dtu_center_vertical);
M_vbap_dtu_center_vertical = M_vbap_dtu_center_vertical.M;
M_vbap_dtu_noscat_vertical = load(file_vbap_dtu_noscat_vertical);
M_vbap_dtu_noscat_vertical = M_vbap_dtu_noscat_vertical.M;
M_p_dtu_center_allrad_vertical = load(file_p_dtu_center_allrad_vertical);
M_p_dtu_center_allrad_vertical = M_p_dtu_center_allrad_vertical.M;
M_p_dtu_noscat_allrad_vertical = load(file_p_dtu_noscat_allrad_vertical);
M_p_dtu_noscat_allrad_vertical = M_p_dtu_noscat_allrad_vertical.M;

% Simulations with encoding
M_dtu_encoding_radinv_center_allrad = load(file_dtu_encoding_radinv_center_allrad);
M_dtu_encoding_radinv_center_allrad = M_dtu_encoding_radinv_center_allrad.M;
M_dtu_encoding_radinv_noscat_allrad = load(file_dtu_encoding_radinv_noscat_allrad);
M_dtu_encoding_radinv_noscat_allrad = M_dtu_encoding_radinv_noscat_allrad.M;
M_dtu_encoding_radinv_center_sad = load(file_dtu_encoding_radinv_center_sad);
M_dtu_encoding_radinv_center_sad = M_dtu_encoding_radinv_center_sad.M;
M_dtu_encoding_radinv_noscat_sad = load(file_dtu_encoding_radinv_noscat_sad);
M_dtu_encoding_radinv_noscat_sad = M_dtu_encoding_radinv_noscat_sad.M;
M_dtu_encoding_dtumic2D_radinv_center_allrad = load(file_dtu_encoding_dtumic2D_radinv_center_allrad);
M_dtu_encoding_dtumic2D_radinv_center_allrad = M_dtu_encoding_dtumic2D_radinv_center_allrad.M;
M_dtu_encoding_dtumic3D_radinv_center_allrad = load(file_dtu_encoding_dtumic3D_radinv_center_allrad);
M_dtu_encoding_dtumic3D_radinv_center_allrad = M_dtu_encoding_dtumic3D_radinv_center_allrad.M;
M_dtu_encoding_dtumic3D_radinv_center_allrad_normal = load(file_dtu_encoding_dtumic3D_radinv_center_allrad_normal);
M_dtu_encoding_dtumic3D_radinv_center_allrad_normal = M_dtu_encoding_dtumic3D_radinv_center_allrad_normal.M;
M_dtu_encoding_dtumic3D_radinv_center_mmd = load(file_dtu_encoding_dtumic3D_radinv_center_mmd);
M_dtu_encoding_dtumic3D_radinv_center_mmd = M_dtu_encoding_dtumic3D_radinv_center_mmd.M;
M_dtu_encoding_dtumicMOA_radinv_center_allrad = load(file_dtu_encoding_dtumicMOA_radinv_center_allrad);
M_dtu_encoding_dtumicMOA_radinv_center_allrad = M_dtu_encoding_dtumicMOA_radinv_center_allrad.M;
M_dtu_encoding_dtumic2D_radinv_noscat_allrad = load(file_dtu_encoding_dtumic2D_radinv_noscat_allrad);
M_dtu_encoding_dtumic2D_radinv_noscat_allrad = M_dtu_encoding_dtumic2D_radinv_noscat_allrad.M;
M_dtu_encoding_dtumic3D_radinv_noscat_allrad = load(file_dtu_encoding_dtumic3D_radinv_noscat_allrad);
M_dtu_encoding_dtumic3D_radinv_noscat_allrad = M_dtu_encoding_dtumic3D_radinv_noscat_allrad.M;
M_dtu_encoding_dtumic3D_radinv_noscat_allrad_normal = load(file_dtu_encoding_dtumic3D_radinv_noscat_allrad_normal);
M_dtu_encoding_dtumic3D_radinv_noscat_allrad_normal = M_dtu_encoding_dtumic3D_radinv_noscat_allrad_normal.M;
M_dtu_encoding_dtumic3D_radinv_noscat_mmd = load(file_dtu_encoding_dtumic3D_radinv_noscat_mmd);
M_dtu_encoding_dtumic3D_radinv_noscat_mmd = M_dtu_encoding_dtumic3D_radinv_noscat_mmd.M;
M_dtu_encoding_dtumicMOA_radinv_noscat_allrad = load(file_dtu_encoding_dtumicMOA_radinv_noscat_allrad);
M_dtu_encoding_dtumicMOA_radinv_noscat_allrad = M_dtu_encoding_dtumicMOA_radinv_noscat_allrad.M;

% Concatenate matrices
M_p_noscat = cat(5,M_p_ref_noscat,M_p_dodeca_noscat,M_p_tdesign_noscat,M_p_tdesign_first_order_noscat,M_p_dtu_noscat_sad,M_p_dtu_noscat_allrad,M_p_dtu_noscat_allrad_normal,M_p_dtu_noscat_mmd);
M_p_center = cat(5,M_p_ref_center,M_p_dodeca_center,M_p_tdesign_center,M_p_tdesign_first_order_center,M_p_dtu_center_sad,M_p_dtu_center_allrad,M_p_dtu_center_allrad_normal,M_p_dtu_center_mmd);
M_vbap_center = cat(5,M_vbap_dodeca_center,M_vbap_tdesign_center,M_vbap_tdesign_first_order_center,M_vbap_dtu_center);
M_vbap_noscat = cat(5,M_vbap_dodeca_noscat,M_vbap_tdesign_noscat,M_vbap_tdesign_first_order_noscat,M_vbap_dtu_noscat);
M_encoding_radinv_center = cat(5,M_dodeca_encoding_radinv_center,M_tdesign_encoding_radinv_center,M_tdesign_first_order_encoding_radinv_center,M_dtu_encoding_radinv_center_allrad,M_dtu_encoding_radinv_center_sad,M_dtu_encoding_dtumic2D_radinv_center_allrad,M_dtu_encoding_dtumic3D_radinv_center_allrad,M_dtu_encoding_dtumicMOA_radinv_center_allrad,M_dtu_encoding_dtumic3D_radinv_center_allrad_normal,M_dtu_encoding_dtumic3D_radinv_center_mmd);
M_encoding_radinv_noscat = cat(5,M_dodeca_encoding_radinv_noscat,M_tdesign_encoding_radinv_noscat,M_tdesign_first_order_encoding_radinv_noscat,M_dtu_encoding_radinv_noscat_allrad,M_dtu_encoding_radinv_noscat_sad,M_dtu_encoding_dtumic2D_radinv_noscat_allrad,M_dtu_encoding_dtumic3D_radinv_noscat_allrad,M_dtu_encoding_dtumicMOA_radinv_noscat_allrad,M_dtu_encoding_dtumic3D_radinv_noscat_allrad_normal,M_dtu_encoding_dtumic3D_radinv_noscat_mmd);
M_p_dtu_vertical = cat(5,M_p_dtu_center_allrad_vertical,M_p_dtu_noscat_allrad_vertical,M_vbap_dtu_center_vertical,M_vbap_dtu_noscat_vertical,M_p_reference_center_vertical,M_p_reference_noscat_vertical);


%% Calculate sound pressure errors, e = 20log(|p_repr/p_ref|), and differences, delta_e = e_head - e_no_head.
%  The above equations are eq. (17) and (18) in [1], respectively.

%{
Indexes for the error matrices, e.g. M(:,:,n)

db_error_spl & db_error_center_noscat_diff
1: reference
2: dodeca
3: t-design 4th order
4: t-design 1st order
5: DTU SAD
6: DTU AllRAD amplitude norm
7: DTU AllRAD power norm
8: DTU MMD

db_error_spl_vbap
1: dodeca
2: t-design 4th order
3: t-design 1st order
4: DTU

db_error_spl_encoding_radinv
1: dodeca
2: t-design 4th order
3: t-design 1st order
4: DTU EigenMike AllRAD
5: DTU EigenMike SAD
6: DTU 2D dtumic AllRAD
7: DTU 3D dtumic AllRAD
8: DTU MOA dtumic AllRAD
9: DTU 3D dtumic AllRAD power norm
10: DTU 3D dtumic MMD

db_error_spl_vertical & db_error_center_noscat_diff_vertical
1: AllRAD center
2: AllRAD noscat
3: VBAP center
4: VBAP noscat
5: ref center
6: ref noscat
%}

% Errors
db_error_spl_noscat = squeeze(20*log10(abs(M_p_noscat(7,:,:,:,:)./M_p_noscat(7,:,:,:,1))));
db_error_spl_center = squeeze(20*log10(abs(M_p_center(1,:,:,:,:)./M_p_center(1,:,:,:,1))));
db_error_spl_encoding_radinv_center = squeeze(20*log10(abs(M_encoding_radinv_center(1,:,:,:,:)./M_p_center(1,:,:,:,1))));
db_error_spl_encoding_radinv_noscat = squeeze(20*log10(abs(M_encoding_radinv_noscat(7,:,:,:,:)./M_p_noscat(7,:,:,:,1))));
db_error_spl_vbap_center = squeeze(20*log10(abs(M_vbap_center(1,:,:,:,:)./M_p_center(1,:,:,:,1))));
db_error_spl_vbap_noscat = squeeze(20*log10(abs(M_vbap_noscat(7,:,:,:,:)./M_p_noscat(7,:,:,:,1))));
db_error_spl_vertical = squeeze(20*log10(abs(M_p_dtu_vertical(1,:,:,:,:)./M_p_dtu_vertical(1,:,:,:,5))));
db_error_spl_vertical(:,:,[2 4]) = squeeze(20*log10(abs(M_p_dtu_vertical(7,:,:,:,[2 4])./M_p_dtu_vertical(7,:,:,:,6))));

% Shift the fields, so that the plots run from 0-360 degrees
db_error_spl_center = circshift(db_error_spl_center,ceil(size(db_error_spl_center,1)/2),1);
db_error_spl_noscat = circshift(db_error_spl_noscat,ceil(size(db_error_spl_noscat,1)/2),1);
db_error_spl_encoding_radinv_center = circshift(db_error_spl_encoding_radinv_center,ceil(size(db_error_spl_encoding_radinv_center,1)/2),1);
db_error_spl_encoding_radinv_noscat = circshift(db_error_spl_encoding_radinv_noscat,ceil(size(db_error_spl_encoding_radinv_noscat,1)/2),1);
db_error_spl_vbap_center = circshift(db_error_spl_vbap_center,ceil(size(db_error_spl_vbap_center,1)/2),1);
db_error_spl_vbap_noscat = circshift(db_error_spl_vbap_noscat,ceil(size(db_error_spl_vbap_noscat,1)/2),1);
db_error_spl_vertical = circshift(db_error_spl_vertical,ceil(size(db_error_spl_vertical,1)/2),1);

% Calculate the difference between centered head and noscat
db_error_center_noscat_diff = db_error_spl_center - db_error_spl_noscat;
db_error_center_noscat_diff_vbap = db_error_spl_vbap_center - db_error_spl_vbap_noscat;
db_error_center_noscat_diff_vertical = db_error_spl_vertical(:,:,[1 3]) - db_error_spl_vertical(:,:,[2 4]);


%% Plot f/DOA (only 3rd, 4th, and DTU)

plot_browser = 0; % 0 or 1 to open plot browser for plots
print_version = 1; % 0 or 1, make labels and plots bigger for printed versions


% Sound pressure error, head in the center
figure_sound_pressure_error_against_freq_and_doa_center = figure;
if (plot_browser)
    plotbrowser('on');
end
ax = axes;
if (print_version)
    %h = suptitle(sprintf("Sound pressure error, head in the center\n\n"));
    %set(h,'FontSize',80);
else
    h = suptitle(sprintf("Sound pressure error, head in the center\n\n"));
    set(h,'FontSize',getGlobalSuptitleFontSize());
end
           subplot(8,1,1,ax); plot_spl(ax,"Dodecahedron, 20 speakers, HOA, N=3, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_center(:,:,2)),[],"0",[],print_version,[],3);
ax = axes; subplot(8,1,2,ax); plot_spl(ax,"Dodecahedron, 20 speakers, VBAP, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_vbap_center(:,:,1)),[],"0",[],print_version,[],3);
ax = axes; subplot(8,1,3,ax); plot_spl(ax,"T-design, 36 speakers, HOA, N=4, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_center(:,:,3)),[],"0",[],print_version,[],4); 
ax = axes; subplot(8,1,4,ax); plot_spl(ax,"T-design, 36 speakers, VBAP, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_vbap_center(:,:,2)),[],"0",[],print_version,[],4);
ax = axes; subplot(8,1,5,ax); plot_spl(ax,"DTU, 64 speakers, HOA, N=5, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_center(:,:,7)),[],"0",[],print_version,[],5); 
ax = axes; subplot(8,1,6,ax); plot_spl(ax,"DTU, 64 speakers, VBAP, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_vbap_center(:,:,4)),[],"0",[],print_version,[],5);
ax = axes; subplot(8,1,7,ax); plot_spl(ax,"DTU, 64 speakers, HOA, N=5, elevation",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_vertical(:,:,1)),[],"0",[],print_version,[],5);
ax = axes; subplot(8,1,8,ax); plot_spl(ax,"DTU, 64 speakers, VBAP, elevation",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_vertical(:,:,3)),[],"0",[],print_version,[],5);


% Sound pressure error, no scatterer
figure_sound_pressure_error_against_freq_and_doa_noscat = figure;
if (plot_browser)
    plotbrowser('on');
end
ax = axes;
if (print_version)
    %h = suptitle(sprintf("Sound pressure error, no scatterer\n\n"));
    %set(h,'FontSize',80);
else
    h = suptitle(sprintf("Sound pressure error, no scatterer\n\n"));
    set(h,'FontSize',getGlobalSuptitleFontSize());
end
           subplot(8,1,1,ax); plot_spl(ax,"Dodecahedron, 20 speakers, HOA, N=3, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_noscat(:,:,2)),[],"0",[],print_version,[],3);
ax = axes; subplot(8,1,2,ax); plot_spl(ax,"Dodecahedron, 20 speakers, VBAP, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_vbap_noscat(:,:,1)),[],"0",[],print_version,[],3);
ax = axes; subplot(8,1,3,ax); plot_spl(ax,"T-design, 36 speakers, HOA, N=4, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_noscat(:,:,3)),[],"0",[],print_version,[],4);
ax = axes; subplot(8,1,4,ax); plot_spl(ax,"T-design, 36 speakers, VBAP, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_vbap_noscat(:,:,2)),[],"0",[],print_version,[],4);
ax = axes; subplot(8,1,5,ax); plot_spl(ax,"DTU, 64 speakers, HOA, N=5, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_noscat(:,:,7)),[],"0",[],print_version,[],5);
ax = axes; subplot(8,1,6,ax); plot_spl(ax,"DTU, 64 speakers, VBAP, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_vbap_noscat(:,:,4)),[],"0",[],print_version,[],5);
ax = axes; subplot(8,1,7,ax); plot_spl(ax,"DTU, 64 speakers, HOA, N=5, elevation",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_vertical(:,:,2)),[],"0",[],print_version,[],5);
ax = axes; subplot(8,1,8,ax); plot_spl(ax,"DTU, 64 speakers, VBAP, elevation",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_vertical(:,:,4)),[],"0",[],print_version,[],5);


% Difference between centered head and noscat
figure_sound_pressure_error_against_freq_and_doa_difference = figure;
if (plot_browser)
    plotbrowser('on');
end
ax = axes;
if (print_version)
    %h = suptitle(sprintf("Difference between head in the field and no head\n\n"));
    %set(h,'FontSize',80);
else
    h = suptitle(sprintf("Difference between head in the field and no head\n\n"));
    set(h,'FontSize',getGlobalSuptitleFontSize());
end
           subplot(8,1,1,ax); plot_spl(ax,"Dodecahedron, 20 speakers, HOA, N=3, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_center_noscat_diff(:,:,2)),[],"0",[],print_version);
ax = axes; subplot(8,1,2,ax); plot_spl(ax,"Dodecahedron, 20 speakers, VBAP, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_center_noscat_diff_vbap(:,:,1)),[],"0",[],print_version);
ax = axes; subplot(8,1,3,ax); plot_spl(ax,"T-design, 36 speakers, HOA, N=4, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_center_noscat_diff(:,:,3)),[],"0",[],print_version);
ax = axes; subplot(8,1,4,ax); plot_spl(ax,"T-design, 36 speakers, VBAP, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_center_noscat_diff_vbap(:,:,3)),[],"0",[],print_version);
ax = axes; subplot(8,1,5,ax); plot_spl(ax,"DTU, 64 speakers, HOA, N=5, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_center_noscat_diff(:,:,7)),[],"0",[],print_version);
ax = axes; subplot(8,1,6,ax); plot_spl(ax,"DTU, 64 speakers, VBAP, azimuth",f_vec,dop_azi_vec*180/pi,squeeze(db_error_center_noscat_diff_vbap(:,:,4)),[],"0",[],print_version);
ax = axes; subplot(8,1,7,ax); plot_spl(ax,"DTU, 64 speakers, HOA, N=5, elevation",f_vec,dop_azi_vec*180/pi,squeeze(db_error_center_noscat_diff_vertical(:,:,1)),[],"0",[],print_version);
ax = axes; subplot(8,1,8,ax); plot_spl(ax,"DTU, 64 speakers, VBAP, elevation",f_vec,dop_azi_vec*180/pi,squeeze(db_error_center_noscat_diff_vertical(:,:,2)),[],"0",[],print_version);


% Sound pressure error, head in the center with encoding
figure_sound_pressure_error_center_encoding = figure;
if (plot_browser)
    plotbrowser('on');
end
ax = axes;
if (print_version)
    %h = suptitle(sprintf("Sound pressure error, head in the center\n\n"));
    %set(h,'FontSize',80);
else
    h = suptitle(sprintf("Sound pressure error, head in the center, encoding\n\n"));
    set(h,'FontSize',getGlobalSuptitleFontSize());
end
           subplot(6,1,1,ax); plot_spl(ax,"Dodecahedron, 20 speakers, N=3",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_center(:,:,2)),[],"0",[],print_version,[],3);
ax = axes; subplot(6,1,2,ax); plot_spl(ax,"Dodecahedron, with Eigenmike encoding",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_encoding_radinv_center(:,:,1)),[],"0",[],print_version,[],3);
ax = axes; subplot(6,1,3,ax); plot_spl(ax,"T-design, 36 speakers, N=4",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_center(:,:,3)),[],"0",[],print_version,[],4); 
ax = axes; subplot(6,1,4,ax); plot_spl(ax,"T-design, with Eigenmike encoding",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_encoding_radinv_center(:,:,2)),[],"0",[],print_version,[],4); 
ax = axes; subplot(6,1,5,ax); plot_spl(ax,"DTU, 64 speakers, N=5",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_center(:,:,7)),[],"0",[],print_version,[],5);
ax = axes; subplot(6,1,6,ax); plot_spl(ax,"DTU, with DTU microphone encoding",f_vec,dop_azi_vec*180/pi,squeeze(db_error_spl_encoding_radinv_center(:,:,9)),[],"0",[],print_version,[],5);


%% Save pressure error figures
f_doa_position = [100 100 500 500];
paper_position = [0.25 0.25 5 9];
paper_position_centimeters = [0.25 0.25 21/2 25];
image_path = "images/";

set(figure_sound_pressure_error_against_freq_and_doa_center, 'PaperOrientation', 'portrait'); set(figure_sound_pressure_error_against_freq_and_doa_center, 'PaperType', 'a4'); set(figure_sound_pressure_error_against_freq_and_doa_center, 'PaperPositionMode', 'auto'); set(figure_sound_pressure_error_against_freq_and_doa_center, 'PaperUnits', 'centimeters'); set(figure_sound_pressure_error_against_freq_and_doa_center,'PaperPosition',paper_position_centimeters); set(figure_sound_pressure_error_against_freq_and_doa_center,'Position',f_doa_position);      print(figure_sound_pressure_error_against_freq_and_doa_center, image_path + "pressure_error_center", '-dpng', '-r150')
set(figure_sound_pressure_error_against_freq_and_doa_noscat, 'PaperOrientation', 'portrait'); set(figure_sound_pressure_error_against_freq_and_doa_noscat, 'PaperType', 'a4'); set(figure_sound_pressure_error_against_freq_and_doa_noscat, 'PaperPositionMode', 'auto'); set(figure_sound_pressure_error_against_freq_and_doa_noscat, 'PaperUnits', 'centimeters'); set(figure_sound_pressure_error_against_freq_and_doa_noscat,'PaperPosition',paper_position_centimeters); set(figure_sound_pressure_error_against_freq_and_doa_noscat,'Position',f_doa_position);      print(figure_sound_pressure_error_against_freq_and_doa_noscat, image_path + "pressure_error_noscat", '-dpng', '-r150')
set(figure_sound_pressure_error_against_freq_and_doa_difference, 'PaperOrientation', 'portrait'); set(figure_sound_pressure_error_against_freq_and_doa_difference, 'PaperType', 'a4'); set(figure_sound_pressure_error_against_freq_and_doa_difference, 'PaperPositionMode', 'auto'); set(figure_sound_pressure_error_against_freq_and_doa_difference, 'PaperUnits', 'centimeters'); set(figure_sound_pressure_error_against_freq_and_doa_difference,'PaperPosition',paper_position_centimeters); set(figure_sound_pressure_error_against_freq_and_doa_difference,'Position',f_doa_position);      print(figure_sound_pressure_error_against_freq_and_doa_difference, image_path + "pressure_error_difference", '-dpng', '-r150')
set(figure_sound_pressure_error_center_encoding, 'PaperOrientation', 'portrait'); set(figure_sound_pressure_error_center_encoding, 'PaperType', 'a4'); set(figure_sound_pressure_error_center_encoding, 'PaperPositionMode', 'auto'); set(figure_sound_pressure_error_center_encoding, 'PaperUnits', 'centimeters'); set(figure_sound_pressure_error_center_encoding,'PaperPosition',paper_position_centimeters); set(figure_sound_pressure_error_center_encoding,'Position',f_doa_position);      print(figure_sound_pressure_error_center_encoding, image_path + "pressure_error_encoding", '-dpng', '-r150')


%% Functions

function [] = plot_spl(ax,stitle,X,Y,Z,color_axis,bar_label,error_lims_pressure_db,print_version,error_limits,N)
    surf(ax,X,Y,Z);
    % Draw the x dB error line
    if (exist('error_lims_pressure_db') && ~isempty(error_lims_pressure_db))
        hold on
        for limit_index = 1:size(error_lims_pressure_db,2)
            limit = error_limits(limit_index);
            line(error_lims_pressure_db(:,limit_index),Y,100*ones(size(Y)),'color','r','linewidth',1);
        end
    end
    % Draw the Ambisonic frequency limit: f_lim <= (c*N)/(2pi*r)
    if (exist('N') && N ~= 0)
        f_line = 343*N/(2*pi*0.1);
        line(f_line*ones(size(Y)),Y,100*ones(size(Y)),'color','w','linewidth',1)
    end
    if (exist('print_version') && print_version ~= 0)
        ax.FontSize = 8;
        title(stitle,'FontSize',8,'Interpreter','latex');
        set(ax,'TickLabelInterpreter','latex')
    else
        ax.FontSize = getGlobalPlotFontSize();
        title(stitle,'FontSize',getGlobalTitleFontSize());
    end
    view(2);
    shading interp;
    if (exist('color_axis') && ~isempty(color_axis))
        caxis(color_axis);
    else
        caxis([-10 10]);
    end
    xlim([50 10000]);
    ylim([0 360]);
    xticks([0 100 1000 10000]);
    xticklabels({'0','$10^2$','$10^3$','$10^4$'})
    yticks([0 90 180 270 360]);
    yticklabels({'-90, r','0, f','90, l','180, b','-90, r'})
    if (contains(stitle,'elevation'))
        yticklabels({'0, r','-90, d','0, l','90, u','0, r'})
    end
    set(ax,'Xscale','log')
    set(ax,'TickDir','out')
    ax.TickLength = [0.02, 0.02]; % Make tick marks longer.
    if (exist('print_version') && print_version ~= 0)
        if (stitle == "T-design, 4th order" || stitle == "2D, 32 loudspeakers" || stitle == "2D VBAP, 32 speakers" || contains(stitle,'DTU, 64 speakers, N=5, with') || contains(stitle,'Vertical, no sphere') || contains(stitle,'VBAP, no sphere') || contains(stitle,'DTU, 64 speakers, VBAP, elevation') || contains(stitle,'DTU, with DTU'))
            bar = colorbar;
            if(exist('bar_label') && bar_label ~= "0")
                ylabel(bar, bar_label,'Interpreter','latex');
            else
                ylabel(bar, "Error (dB)",'Interpreter','latex');
            end
            set(bar,'TickLabelInterpreter','latex')
            y = ylabel("Direction of Arrival (degrees)",'Interpreter','latex');
            xlabel("Frequency (Hz)",'Interpreter','latex');
            set(get(bar,'Label'),'FontSize',9)
            set(y,'FontSize',9);
            set(bar,'Position',[0.93, 0.115, 0.015, 0.8])
            set(get(bar,'Label'),'Position',[2.5 0])
            set(y,'Units','normalized','Position',[-0.09, 8.5, 0])% [0.25 0.25 21/2 25];
            if (contains(stitle,'DTU, 64 speakers, VBAP, elevation'))
                set(y,'Units','normalized','Position',[-0.12, 6, 0])% [0.25 0.25 21/2 25];
            elseif (contains(stitle,'Vertical, no sphere') || contains(stitle,'VBAP, no sphere'))
                set(y,'Units','normalized','Position',[-0.1, 2.5, 0])% [0.25 0.25 21/2 25];
            elseif (contains(stitle,'DTU, 64 speakers, N=5, with') || contains(stitle,'DTU, with DTU'))
                set(y,'Units','normalized','Position',[-0.12, 4  , 0])% [0.25 0.25 21/2 25];
            end
        end
    else
        if (stitle == "T-design, 4th order" || stitle == "2D, 16 loudspeakers" || stitle == "2D VBAP, 16 speakers" || stitle == "DTU, AllRAD, N=7")% || stitle == "DTU, SAD, N=7")% || stitle == "DTU, N=7")
            bar = colorbar;
            if(exist('bar_label') && bar_label ~= "0")
                ylabel(bar, bar_label);
            else
                ylabel(bar, "Error (dB)");
            end
            y = ylabel("Direction of Arrival (degrees)");
            xlabel("Frequency (Hz)");
            set(bar,'Position',[0.94, 0.05, 0.01, 0.9])
            set(get(bar,'Label'),'Position',[3 0])
            set(y,'Units','Normalized','Position',[-0.07, 6.4, 0])
        end
    end
end

function value = getGlobalPlotFontSize()
    global plot_font_size
    value = plot_font_size;
end
function value = getGlobalSuptitleFontSize()
    global suptitle_font_size
    value = suptitle_font_size;
end
function value = getGlobalTitleFontSize()
    global title_font_size
    value = title_font_size;
end
function value = getGlobalLegendFontSize()
    global legend_font_size
    value = legend_font_size;
end


%% References
%{
    [1] Pajunen, L., Politis, A., Pulkki, V., Vaalgamaa, M., Strömmer, S., 2020,
        Effects of rigid spherical scatterer on spatial audio reproduction quality.
        Submitted for Audio Engineering Society Convention 148.

%}































