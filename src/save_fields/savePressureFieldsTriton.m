function [] = savePressureFieldsTriton(workfolder,outputfile,N_diffuse_reps,setting,rng_seed,delay_distance,encoding_regularization)
%SAVEPRESSUREFIELDSTRITON Calculate pressure fields and save them to a
%                         matrix (call this function from a shell script).
%
%   Call the correct savePressureFieldsSETTING.m function defined in the
%   'setting' argument. This function is called from a shell script that is
%   executed in Triton high-computing cluster.
%
% ARGUMENTS:
% workfolder - the root folder with all the necessary matlab libraries
%              (mainly from https://github.com/polarch and this library)
% outputfile - the calculated pressure matrix is stored into this file
% N_diffuse_reps - number of repetitions for the diffuse field calculations
%                  (use 0 for non-diffuse field calculations)
% setting - reference or reproduction loudspeaker setting name
% rng_seed - seed for the random number generator, this should be unique
%            for each diffuse field script launched by the shell script
% delay_distance - the distance from the center the delay time should be
%                  adjusted to, in meters, (0.02, 0.05 or 0.10 meters)
%                  OPTIONAL (used in XXX_delay setups)
% encoding_regularization - {'RADINV','REGLS','SOFTLIM','NONE'} for a
%                           regularization filter when virtual EigenMike is
%                           used, OPTIONAL

try
    % Add folders to path
    addpath(genpath(workfolder));

    % Seed the random number generator
    rng(rng_seed);
    
    % Path to 'data' folder (default is that the code is run
    % in Spherical-scatterer-effects folder, i.e. pwd = '../Spherical-scatterer-effects')
    data_path = "data/";

    outputfile
    N_diffuse_reps

    % Variables
    %head_diameter = 0.2;
    head_diameter = 0.2;
    %head_pos = [pi 0 head_diameter/2]; % head position from the center [azi ele r]
    head_pos = [pi 0 0]; % head position from the center [azi ele r]
    head_direction = pi/2;  % azimuth pointing direction
    r_scat = head_diameter/2;
    %D = head_diameter/2 + head_pos(3) + 0.1; % playback radius for calculations, in meters
    D = 1 + head_pos(3); % playback radius for calculations, in meters
    d = 0.01; % grid resolution, in meters
    f_vec = 50:50:10000; % frequencies
    %f_vec = [100 500 1000 4000];
    %f_vec = 500;
    %dop_azi_vec = (0:90:360)*pi/180;
    doa_azi_vec = (-180:1:180)*pi/180; % azimuth angles for directions of arrivals
    doa_vec = doa_azi_vec;
    doa_vec(2,:) = 0;
    
    % Elevated circle
    doa_azi_vec = zeros(1,361);
    doa_azi_vec(1:91) = (180-90)*pi/180;
    doa_azi_vec(92:271) = (0-90)*pi/180;
    doa_azi_vec(272:361) = (180-90)*pi/180;
    doa_ele_vec = [0:1:90 89:-1:-90 -89:1:0]*pi/180;
    doa_vec = doa_azi_vec;
    doa_vec(2,:) = doa_ele_vec;
    
    file = sprintf("%s/%s",pwd,outputfile);
    c = 343;
    
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
    offset_no_scat = [0 0; 0 0.02; 0 0.05; 0 0.1; 0.02 0; 0.05 0; 0.1 0; -0.1 0];
    offset_amnt_no_scat = round(offset_no_scat/d);
    offset_amnt = round(offset/d);
    % Points
    points_no_scat = [(offset_amnt_no_scat(:,1) + center_x) (offset_amnt_no_scat(:,2) + center_y)]
    points = [(offset_amnt(:,1) + center_x + head_x_amnt) (offset_amnt(:,2) + center_y + head_y_amnt)]
    % Distance from the center
    dist_p_no_scat = sqrt((offset_no_scat(:,1)).^2+(offset_no_scat(:,2)).^2);
    max_dist_p_no_scat = max(dist_p_no_scat);
    dist_p = sqrt((offset(:,1) + head_x).^2+(offset(:,2) + head_y).^2);
    max_dist_p = max(dist_p);

    % Indexes for single poinst in the field
    indexes_no_scat = sub2ind(size(X), points_no_scat(:,2), points_no_scat(:,1));
    indexes = sub2ind(size(X), points(:,2), points(:,1));

    % Plane wave directions for diffuse field calculations
    [~, pw_dirs] = getTdesign(21);
    pw_dirs = pw_dirs';

    switch (lower(setting))
        case 'ref'
            % Non-diffuse
            %savePressureFieldsReference(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec);
            savePressureFieldsReference(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,indexes,max_dist_p);

        case 'dodeca'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            %savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p)

        case 'tdesign'
            N_tdesign = 4;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            %savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,0,0,N_tdesign)

        case 'tdesign_first_order'
            N_tdesign = 1;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            %savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p)

        case 'ref_diffuse'
            % Diffuse
            savePressureFieldsReference(file,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,indexes,max_dist_p,N_diffuse_reps)

        case 'dodeca_diffuse'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            decoding_method = 'sad';
            % Diffuse
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,N_diffuse_reps);

        case 'tdesign_diffuse'
            N_tdesign = 4;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            % Diffuse
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,N_diffuse_reps,0,N_tdesign);

        case 'tdesign_first_order_diffuse'
            N_tdesign = 1;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            % Diffuse
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,N_diffuse_reps);

        case 'dodeca_vbap_diffuse'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            N_speakers = getLayoutAmbisonicOrder(speaker_dirs*180/pi);
            % Diffuse
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p,N_diffuse_reps);

        case 'tdesign_vbap_diffuse'
            N_tdesign = 4;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            % Diffuse
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p,N_diffuse_reps);

        case 'tdesign_first_order_vbap_diffuse'
            N_tdesign = 1;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            % Diffuse
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p,N_diffuse_reps);

        case '2d_vbap_16_diffuse'
            hop_size = 2*pi/16;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((16-1)/2);
            % Diffuse
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p,N_diffuse_reps)
        
        case '2d_vbap_32_diffuse'
            hop_size = 2*pi/32;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((32-1)/2);
            % Diffuse
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p,N_diffuse_reps)

        case 'ref_noscat'
            % Non-diffuse
            %savePressureFieldsReferenceNoScat(file,X,Y,D,d,f_vec,doa_vec);
            savePressureFieldsReferenceNoScat(file,X,Y,D,d,f_vec,doa_vec,indexes_no_scat,max_dist_p_no_scat);

        case 'dodeca_noscat'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            %savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat)

        case 'tdesign_noscat'
            N_tdesign = 4;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            %savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat,0,0,N_tdesign)

        case 'tdesign_first_order_noscat'
            N_tdesign = 1;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            %savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat)

        case 'dodeca_vbap'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            N_speakers = getLayoutAmbisonicOrder(speaker_dirs*180/pi);
            %savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)

        case 'tdesign_vbap'
            N_tdesign = 4;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            %savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)

        case 'tdesign_first_order_vbap'
            N_tdesign = 1;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            %savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)
    
        case '2d_vbap_16'
            hop_size = 2*pi/16;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((16-1)/2);
            %savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)
        
        case '2d_vbap_32'
            hop_size = 2*pi/32;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((32-1)/2);
            %savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)

        case 'dodeca_vbap_noscat'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            N_speakers = getLayoutAmbisonicOrder(speaker_dirs*180/pi);
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            %savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat)

        case 'tdesign_vbap_noscat'
            N_tdesign = 4;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            %savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat)

        case 'tdesign_first_order_vbap_noscat'
            N_tdesign = 1;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            %savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat)

        case '2d_vbap_16_noscat'
            hop_size = 2*pi/16;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((16-1)/2);
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            %savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat)
        
        case '2d_vbap_32_noscat'
            hop_size = 2*pi/32;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((32-1)/2);
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            %savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat)

        case 'ref_velocity'
            % Non-diffuse
            %saveVelocityFieldsReference(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec);
            saveVelocityFieldsReference(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,indexes,max_dist_p);

        case 'dodeca_velocity'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            %saveVelocityFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            saveVelocityFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p)

        case 'tdesign_velocity'
            N_tdesign = 4;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            %saveVelocityFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            saveVelocityFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,0,N_tdesign)

        case 'tdesign_first_order_velocity'
            N_tdesign = 1;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            %saveVelocityFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            saveVelocityFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p)

        case 'ref_velocity_noscat'
            % Non-diffuse
            %saveVelocityFieldsReferenceNoScat(file,X,Y,D,d,f_vec,doa_vec);
            saveVelocityFieldsReferenceNoScat(file,X,Y,D,d,f_vec,doa_vec,indexes_no_scat,max_dist_p_no_scat);

        case 'dodeca_velocity_noscat'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            %saveVelocityFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            saveVelocityFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat)

        case 'tdesign_velocity_noscat'
            N_tdesign = 4;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            %saveVelocityFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            saveVelocityFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat,0,N_tdesign)

        case 'tdesign_first_order_velocity_noscat'
            N_tdesign = 1;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            %saveVelocityFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            saveVelocityFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat)

        case 'dodeca_vbap_velocity'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            N_speakers = getLayoutAmbisonicOrder(speaker_dirs*180/pi);
            %saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)

        case 'tdesign_vbap_velocity'
            N_tdesign = 4;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            %saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)

        case 'tdesign_first_order_vbap_velocity'
            N_tdesign = 1;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            %saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)
    
        case '2d_vbap_16_velocity'
            hop_size = 2*pi/16;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((16-1)/2);
            %saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)
        
        case '2d_vbap_32_velocity'
            hop_size = 2*pi/32;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((32-1)/2);
            %saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)

        case 'dodeca_vbap_velocity_noscat'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            N_speakers = getLayoutAmbisonicOrder(speaker_dirs*180/pi);
            %saveVelocityFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            saveVelocityFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat)

        case 'tdesign_vbap_velocity_noscat'
            N_tdesign = 4;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            %saveVelocityFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            saveVelocityFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat)

        case 'tdesign_first_order_vbap_velocity_noscat'
            N_tdesign = 1;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            %saveVelocityFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            saveVelocityFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat)

        case '2d_vbap_16_velocity_noscat'
            hop_size = 2*pi/16;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((16-1)/2);
            %saveVelocityFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,N_speakers,source_distance)
            saveVelocityFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat)
        
        case '2d_vbap_32_velocity_noscat'
            hop_size = 2*pi/32;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((32-1)/2);
            %saveVelocityFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,N_speakers,source_distance)
            saveVelocityFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat)
            
        case 'listening_room_allrad'
            speaker_dirs = load(data_path + 'simulated_responses/listening_room_speaker_directions.mat');
            speaker_dirs = speaker_dirs.dirs;
            speaker_dirs(:,1:2) = speaker_dirs(:,1:2)*pi/180;
            source_distance = 2;
            decoding_method = 'allrad';
            % Non-diffuse
            %savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p)
   
        case 'listening_room_allrad_noscat'
            speaker_dirs = load(data_path + 'simulated_responses/listening_room_speaker_directions.mat');
            speaker_dirs = speaker_dirs.dirs;
            speaker_dirs(:,1:2) = speaker_dirs(:,1:2)*pi/180;
            source_distance = 2;
            decoding_method = 'allrad';
            % Non-diffuse
            %savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat)
        
        case 'listening_room_allrad_diffuse'
            speaker_dirs = load(data_path + 'simulated_responses/listening_room_speaker_directions.mat');
            speaker_dirs = speaker_dirs.dirs;
            speaker_dirs(:,1:2) = speaker_dirs(:,1:2)*pi/180;
            source_distance = 2;
            decoding_method = 'allrad';
            % Diffuse
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,N_diffuse_reps);
        
        case 'listening_room_sad'
            speaker_dirs = load(data_path + 'simulated_responses/listening_room_speaker_directions.mat');
            speaker_dirs = speaker_dirs.dirs;
            speaker_dirs(:,1:2) = speaker_dirs(:,1:2)*pi/180;
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            %savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p)
   
        case 'listening_room_sad_noscat'
            speaker_dirs = load(data_path + 'simulated_responses/listening_room_speaker_directions.mat');
            speaker_dirs = speaker_dirs.dirs;
            speaker_dirs(:,1:2) = speaker_dirs(:,1:2)*pi/180;
            source_distance = 2;
            decoding_method = 'sad';
            % Non-diffuse
            %savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat)
        
        case 'listening_room_sad_diffuse'
            speaker_dirs = load(data_path + 'simulated_responses/listening_room_speaker_directions.mat');
            speaker_dirs = speaker_dirs.dirs;
            speaker_dirs(:,1:2) = speaker_dirs(:,1:2)*pi/180;
            source_distance = 2;
            decoding_method = 'sad';
            % Diffuse
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,pw_dirs,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,N_diffuse_reps);
            
        case 'listening_room_vbap'
            speaker_dirs = load(data_path + 'simulated_responses/listening_room_speaker_directions_no_elevation.mat');
            speaker_dirs = speaker_dirs.dirs;
            speaker_dirs(:,1:2) = speaker_dirs(:,1:2)*pi/180;
            source_distance = 2;
            N_speakers = floor((9-1)/2);
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)
        
        case 'listening_room_vbap_noscat'
            speaker_dirs = load(data_path + 'simulated_responses/listening_room_speaker_directions_no_elevation.mat');
            speaker_dirs = speaker_dirs.dirs;
            speaker_dirs(:,1:2) = speaker_dirs(:,1:2)*pi/180;
            source_distance = 2;
            N_speakers = floor((9-1)/2);
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat)
            
        case 'ref_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes = indexes([1])
            savePressureFieldsReference(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,indexes,max_dist_p,0,delay);

        case 'dodeca_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes = indexes([1])
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            decoding_method = 'sad';
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,0,delay)

        case 'tdesign_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes = indexes([1])
            N_tdesign = 4;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,0,delay,N_tdesign)

        case 'tdesign_first_order_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes = indexes([1])
            N_tdesign = 1;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,0,delay)
 
        case 'dodeca_vbap_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes = indexes([1])
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            N_speakers = getLayoutAmbisonicOrder(speaker_dirs*180/pi);
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p,0,delay)

        case 'tdesign_vbap_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes = indexes([1])
            N_tdesign = 4;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p,0,delay)

        case 'tdesign_first_order_vbap_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes = indexes([1])
            N_tdesign = 1;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p,0,delay)
    
        case '2d_vbap_16_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes = indexes([1])
            hop_size = 2*pi/16;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((16-1)/2);
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p,0,delay)
        
        case '2d_vbap_32_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes = indexes([1])
            hop_size = 2*pi/32;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((32-1)/2);
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p,0,delay)
            
        case 'ref_noscat_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes_no_scat = indexes_no_scat([1])
            savePressureFieldsReferenceNoScat(file,X,Y,D,d,f_vec,doa_vec,indexes_no_scat,max_dist_p_no_scat,0,delay);
            
        case 'dodeca_noscat_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes_no_scat = indexes_no_scat([1])
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            decoding_method = 'sad';
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat,0,delay)

        case 'tdesign_noscat_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes_no_scat = indexes_no_scat([1])
            N_tdesign = 4;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat,0,delay,N_tdesign)

        case 'tdesign_first_order_noscat_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes_no_scat = indexes_no_scat([1])
            N_tdesign = 1;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat,0,delay)
 
        case 'dodeca_vbap_noscat_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes_no_scat = indexes_no_scat([1])
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            N_speakers = getLayoutAmbisonicOrder(speaker_dirs*180/pi);
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat,delay)

        case 'tdesign_vbap_noscat_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes_no_scat = indexes_no_scat([1])
            N_tdesign = 4;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat,delay)

        case 'tdesign_first_order_vbap_noscat_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes_no_scat = indexes_no_scat([1])
            N_tdesign = 1;
            N_speakers = N_tdesign;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat,delay)
    
        case '2d_vbap_16_noscat_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes_no_scat = indexes_no_scat([1])
            hop_size = 2*pi/16;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((16-1)/2);
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat,delay)
        
        case '2d_vbap_32_noscat_delay'
            % Delay the center point
            delay = delay_distance/c
            indexes_no_scat = indexes_no_scat([1])
            hop_size = 2*pi/32;
            speaker_dirs(:,1) = 0:hop_size:2*pi-hop_size;
            speaker_dirs(:,2) = 0;
            source_distance = 2;
            N_speakers = floor((32-1)/2);
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat,delay)
            
        case 'dodeca_eigenmike_encoding'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            decoding_method = 'sad';
            %savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,[],max_dist_p,0,0,0,encoding_regularization)
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,0,0,0,encoding_regularization)

        case 'tdesign_eigenmike_encoding'
            N_tdesign = 4;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            %savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,[],max_dist_p,0,0,N_tdesign,encoding_regularization)
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,0,0,N_tdesign,encoding_regularization)
            
        case 'tdesign_first_order_eigenmike_encoding'
            N_tdesign = 1;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            %savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,[],max_dist_p,0,0,0,encoding_regularization)
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,0,0,0,encoding_regularization)
            
        case 'dodeca_eigenmike_encoding_noscat'
            [~, speaker_dirs, ~] = platonicSolid('dodeca');
            source_distance = 2;
            decoding_method = 'sad';
            %savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,[],max_dist_p_no_scat,0,0,0,encoding_regularization)
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat,0,0,0,encoding_regularization)

        case 'tdesign_eigenmike_encoding_noscat'
            N_tdesign = 4;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            %savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,[],max_dist_p_no_scat,0,0,N_tdesign,encoding_regularization)
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat,0,0,N_tdesign,encoding_regularization)

        case 'tdesign_first_order_eigenmike_encoding_noscat'
            N_tdesign = 1;
            [~, speaker_dirs] = getTdesign(2*N_tdesign);
            source_distance = 2;
            decoding_method = 'sad';
            %savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,[],max_dist_p_no_scat,0,0,0,encoding_regularization)
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat,0,0,0,encoding_regularization)

        case 'dtu'
            dtu_setup = load(data_path + 'simulated_responses/DTU_ls_dirs_deg.mat');
            dtu_setup = dtu_setup.ls_dirs_deg;
            speaker_dirs = dtu_setup*pi/180;
            source_distance = 2;
            N_speakers = 5;
            decoding_method = 'allrad_normal';
            %savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,0,0,N_speakers)

        case 'dtu_noscat'
            dtu_setup = load(data_path + 'simulated_responses/DTU_ls_dirs_deg.mat');
            dtu_setup = dtu_setup.ls_dirs_deg;
            speaker_dirs = dtu_setup*pi/180;
            source_distance = 2;
            N_speakers = 5;
            decoding_method = 'allrad_normal';
            %savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat,0,0,N_speakers)
           
        case 'dtu_vbap'
            dtu_setup = load(data_path + 'simulated_responses/DTU_ls_dirs_deg.mat');
            dtu_setup = dtu_setup.ls_dirs_deg;
            speaker_dirs = dtu_setup*pi/180;
            source_distance = 2;
            N_speakers = getLayoutAmbisonicOrder(speaker_dirs);
            savePressureFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)
            
        case 'dtu_vbap_noscat'
            dtu_setup = load(data_path + 'simulated_responses/DTU_ls_dirs_deg.mat');
            dtu_setup = dtu_setup.ls_dirs_deg;
            speaker_dirs = dtu_setup*pi/180;
            source_distance = 2;
            N_speakers = getLayoutAmbisonicOrder(speaker_dirs);
            savePressureFieldsVBAPNoScat(file,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes_no_scat,max_dist_p_no_scat)
            
        case 'dtu_eigenmike_encoding'
            dtu_setup = load(data_path + 'simulated_responses/DTU_ls_dirs_deg.mat');
            dtu_setup = dtu_setup.ls_dirs_deg;
            speaker_dirs = dtu_setup*pi/180;
            source_distance = 2;
            N_speakers = 5;
            decoding_method = 'allrad';
            %savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,[],max_dist_p_no_scat,0,0,0,encoding_regularization)
            savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p,0,0,N_speakers,encoding_regularization)
 
        case 'dtu_eigenmike_encoding_noscat'
            dtu_setup = load(data_path + 'simulated_responses/DTU_ls_dirs_deg.mat');
            dtu_setup = dtu_setup.ls_dirs_deg;
            speaker_dirs = dtu_setup*pi/180;
            source_distance = 2;
            N_speakers = 5;
            decoding_method = 'allrad';
            %savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,[],max_dist_p_no_scat,0,0,0,encoding_regularization)
            savePressureFieldsReproductionNoScat(file,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes_no_scat,max_dist_p_no_scat,0,0,N_speakers,encoding_regularization)

        case 'dtu_velocity'
            dtu_setup = load(data_path + 'simulated_responses/DTU_ls_dirs_deg.mat');
            dtu_setup = dtu_setup.ls_dirs_deg;
            speaker_dirs = dtu_setup*pi/180;
            source_distance = 2;
            decoding_method = 'allrad';
            % Non-diffuse
            %saveVelocityFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance)
            saveVelocityFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p)

        case 'dtu_vbap_velocity'
            dtu_setup = load(data_path + 'simulated_responses/DTU_ls_dirs_deg.mat');
            dtu_setup = dtu_setup.ls_dirs_deg;
            speaker_dirs = dtu_setup*pi/180;
            source_distance = 2;
            N_speakers = getLayoutAmbisonicOrder(speaker_dirs);
            %saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers)
            saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_dist_p)
            
        otherwise
            fprintf("Couldn't find any case for setting: %s\n",lower(setting));
    end

catch e
    fprintf(1,'The identifier was:\n%s\n',e.identifier);
    fprintf(1,'There was an error! The message was:\n%s\n',e.message);
    for i = 1:max(size(e.stack))
        e.stack(i)
    end
    exit(0)
end
end
