function [] = saveVelocityFieldsVBAP(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,speaker_dirs,source_distance,N_speakers,indexes,max_distance,delay)
%SAVEVELOCITYFIELDSVBAP Calculate and save a number of velocity fields
%                       into a file.
%   Calculate velocity fields for plane waves with varying DOA for the 
%   reproduction settings using only Vector Base Amplitude Panning (VBAP).
%   Save the veloicty matrixes into a .mat file
%   for later analysis. The saved 5D matrix follows the form: 
%   M(X,Y,doa,f,3). For example, M(:,:,3,1,2) retrieves the complex
%   velocity field for y-component calculated with the first frequency in f_vec, and
%   the third DOA angle in doa_vec.
%
% ARGUMENTS:
% file - path to the file the pressure matrix is saved into
% head_diameter - diameter of the head in meters
% head_pos - position of the head from the center [azi ele r]
% X and Y - meshgrid for X and Y
% D - playback radius in meters
% d - grid resolution
% f_vec - frequency vector for all desired frequencies
% doa_vec - direction of arrival vector for all desired directions,
%           [azi1 azi2 ... azin; ele1 ele2 ... elen] in radians
% speaker_dirs - loudspeaker directions [azi ele] in radians
% source_distance - distance of the loudspeakers from the center
% N_speakers - the order of the playback setup
% indexes - for calculating only single points in the grid, OPTIONAL
% max_distance - maximum distance of the points from the center to make
%                calculations faster, OPTIONAL (must have 'indexes' set)
% delay - the amount of time delay to the simulated field, OPTIONAL

% NOTE: Here inclination angle is the angle from the +z-axis,
%       elevation angle is the angle from xy-plane towards z-axis

% Acoustics
c = 343;
rho_0 = 1.2;
Z_0 = c*rho_0;

if (head_pos(3) > D)
   error('Head/spherical scatterer is outside the reproduction field.'); 
end

%{
% Check if the file already exists
if (exist(file, 'file') == 2)
   prompt = sprintf("The file '%s' already exists, are you sure you want to overwrite it? Y/N: ",file);
   answer = input(prompt,'s');
   if (length(answer) > 1 || answer ~= 'Y' && answer ~= 'y')
       error('Script stopped by the user.');
   end
end
%}

% Helper function to convert directions from Matlab's azimuth-elevation to
% azimuth-inclination
aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

% Calculations for only the selected points
if (exist('indexes') && ~isempty(indexes))
    M = zeros(size(indexes,1), size(indexes,2), size(doa_vec,2), length(f_vec), 3);
% Calculations for all the points in the grid
else
    M = zeros(size(X,1), size(X,2), size(doa_vec,2), length(f_vec), 3);
end

% Speaker directions
speaker_dirs_incl = aziElev2aziIncl(speaker_dirs); % convert to azi-incl
Kspeakers = size(speaker_dirs_incl,1); % design size

% Shift the loudspeakers according to the head position
[shift_x,shift_y,shift_z] = sph2cart(head_pos(1),head_pos(2),head_pos(3));
shift_x = -shift_x;
shift_y = -shift_y;
shift_z = -shift_z;
shift_amnt_x = round(-shift_x/d);
shift_amnt_y = round(-shift_y/d);
[speakers_x,speakers_y,speakers_z] = sph2cart(speaker_dirs(:,1),speaker_dirs(:,2),ones(Kspeakers,1)*source_distance);
speakers_x = speakers_x + shift_x;
speakers_y = speakers_y + shift_y;
speakers_z = speakers_z + shift_z;
[new_azi,new_ele,new_r] = cart2sph(speakers_x,speakers_y,speakers_z);
speakers_dirs_ele_shifted = cat(2,new_azi,new_ele,new_r);

% Find valid loudspeaker pairs or triplets
if (all(speakers_dirs_ele_shifted(:,2) == 0))
    % 2D
    ls_dirs = speakers_dirs_ele_shifted(:,1)*180/pi';
    ls_groups = findLsPairs(ls_dirs);
else
    % 3D
    ls_dirs = speakers_dirs_ele_shifted(:,1:2)*180/pi';
    [ls_groups, ~] = findLsTriplets(ls_dirs);
end

% Compute inverse matrices for loudspeaker pairs or triplets
layoutInvMtx = invertLsMtx(ls_dirs, ls_groups);

% Compute vbap gains for the required source directions
g_speakers_vec = vbap_amplitude_norm(doa_vec'*180/pi, ls_groups, layoutInvMtx); % compute vbap gains

% Gain boost to set the source amplitude as m0 in the center of
% reproduction
%gain_boost = m0 / real(m0*exp(1i*mod((source_position(3)*f/c),1)*2*pi)*exp(1i*(w*t0*1 - k*sqrt(source_position(3)^2)))/(4*pi*sqrt(source_position(3)^2)));
%gain_boost = m0 / real(m0*exp(1i*mod((source_distance*f/c),1)*2*pi)*exp(1i*(w*t0*1 - k*sqrt(source_distance^2)))/(4*pi*sqrt(source_distance^2)));
g_speakers_vec = 4*pi*source_distance*g_speakers_vec;

% Start timer
tic

f_indx = 1;
for f = f_vec
    fprintf('Frequency %d/%d. ',f,f_vec(end))
    % Phase correction to make the max points of the waves to meet in the
    % center of the reproduction
    %initial_phases = [mod((source_position(3)*f/c),1)*2*pi];
    initial_phases = [mod((source_distance*f/c),1)*2*pi];
    
    toc
    doa_indx = 1;
    for doa = doa_vec
        t0 = 0;
        t = 0; % time point
        if (exist('delay'))
           t = t + delay;
        end
        w = 2*pi*f; % angular frequency
        k = w/c; % wavenumber
        m0 = 1; % magnitude
        r_scat = head_diameter/2; % radius of the scatterer
        N_max = ceil(exp(1)*k*(D)/2); % maximum order of expansion

        % Gains
        g_speakers = g_speakers_vec(doa_indx,:);

        % Apply the gains to K point sources (reproduced sound field)
        frequencies = ones(1,1)*f;

        % Initial velocity matrix
        % Calculations for only the selected points
        if (exist('indexes') && ~isempty(indexes))
            velocity = zeros(1,size(indexes,2));
        % Calculations for all the points in the grid
        else
            velocity = zeros(size(X,1));
        end
        
        % Calculate the pressure field
        for i=1:max(size(g_speakers))
            source_position = speakers_dirs_ele_shifted(i,:);
            m_speakers = g_speakers(i)*m0;

            % Calculations for only the selected points
            if (exist('indexes') && ~isempty(indexes))
                if (exist('max_distance') && max_distance ~= 0)
                    N_max = ceil(exp(1)*k*(max_distance)/2); % maximum order of expansion
                end
                [temp_velocity, ~] = point_sources_scatterer_velocity(t,frequencies,source_position,[m_speakers],initial_phases,X,Y,N_max,D,r_scat,indexes);
                %[temp_velocity, ~] = point_sources_scatterer_velocity(t,frequencies,source_position,[m_speakers],initial_phases,X,Y,N_speakers,D,r_scat,indexes);
            % Calculations for all the points in the grid
            else
                [temp_velocity, ~] = point_sources_scatterer_velocity(t,frequencies,source_position,[m_speakers],initial_phases,X,Y,N_max,D,r_scat);
                %[temp_velocity, ~] = point_sources_scatterer_velocity(t,frequencies,source_position,[m_speakers],initial_phases,X,Y,N_speakers,D,r_scat);
            end
            velocity = velocity + temp_velocity;
        end

        % Shift the velocity field according to the head position
        if (~exist('indexes') || isempty(indexes))
            velocity = circshift(velocity,shift_amnt_x,2);
            velocity = circshift(velocity,shift_amnt_y,1);
        end

        % Store the velocity field
        M(:,:,doa_indx,f_indx,:) = velocity;
        doa_indx = doa_indx + 1;
    end
    f_indx = f_indx + 1;
end
% Write the data to a file
save(file, 'M');
end

