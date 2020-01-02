function [] = saveVelocityFieldsReference(file,head_diameter,head_pos,X,Y,D,d,f_vec,doa_vec,indexes,max_distance,delay)
%SAVEVELOCITYFIELDSREFERENCE Calculate and save a number of velocity fields
%                            into a file.
%   Calculate velocity fields for plane waves with varying DOA. Save the
%   velocity matrixes into a .mat file for later analysis. The saved 5D
%   matrix follows the form: M(X,Y,doa,f,xyz). For example, M(:,:,3,1,2)
%   retrieves the y-values of the complex velocity field calculated with
%   the first frequency in f_vec, and the third DOA angle in doa_vec.
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

% Shift values
[shift_x,shift_y,shift_z] = sph2cart(head_pos(1),head_pos(2),head_pos(3));
shift_x = -shift_x;
shift_y = -shift_y;
shift_z = -shift_z;
shift_amnt_x = round(-shift_x/d);
shift_amnt_y = round(-shift_y/d);

% Values for calculations
f_indx = 1;

% Calculations for only the selected points
if (exist('indexes') && ~isempty(indexes))
    M = zeros(size(indexes,1), size(indexes,2), size(doa_vec,2), length(f_vec), 3);
% Calculations for all the points in the grid
else
    M = zeros(size(X,1), size(X,2), size(doa_vec,2), length(f_vec), 3);
end

% Start timer
tic

for f = f_vec
    fprintf('Frequency %d/%d. ',f,f_vec(end))
    toc
    doa_indx = 1;
    for doa = doa_vec
        t = 0; % time point
        if (exist('delay'))
           t = t + delay;
        end
        w = 2*pi*f; % angular frequency
        k = w/c; % wavenumber
        m0 = 1; % magnitude
        r_scat = head_diameter/2; % radius of the scatterer
        N_max = ceil(exp(1)*k*(D)/2); % maximum order of expansion
        
        % Phase correction to make the max point of the plane wave to meet in the
        % center after the shifting of the field
        [doa_x,doa_y,doa_z] = sph2cart(doa(1),doa(2),1);
        [head_x,head_y,head_z]  = sph2cart(head_pos(1),head_pos(2),head_pos(3));
        doa_cart = [doa_x doa_y doa_z];
        head_cart = [head_x head_y head_z];
        wave_distance = dot(head_cart,doa_cart)/norm(doa_cart);
        phi0 = -mod((wave_distance*f/c),1)*2*pi; % initial phase

        % Velocity
        % Calculations for only the selected points
        if (exist('indexes') && ~isempty(indexes))
            if (exist('max_distance') && max_distance ~= 0)
                N_max = ceil(exp(1)*k*(max_distance)/2); % maximum order of expansion
            end
            [v_cmplx] = plane_wave_scatterer_velocity(t,f,m0,phi0,doa,X,Y,N_max,D,r_scat,indexes);
        % Calculations for all the points in the grid
        else
            [v_cmplx] = plane_wave_scatterer_velocity(t,f,m0,phi0,doa,X,Y,N_max,D,r_scat);
            % Shift the velocity field according to the head position
            v_cmplx = circshift(v_cmplx,shift_amnt_x,2);
            v_cmplx = circshift(v_cmplx,shift_amnt_y,1);
        end
        % Store the velocity field
        M(:,:,doa_indx,f_indx,:) = v_cmplx;
        doa_indx = doa_indx + 1;
    end
    f_indx = f_indx + 1;
end
% Write the data to a file
save(file, 'M');
end

