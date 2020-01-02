function [plane_waves] = tilted_plane_waves(X,Y,frequencies,amplitudes,directions,attenuation)
% TILTED_PLANE_WAVES Compute tilted plane waves based on the angle from the x-axis.

%   PARAMETERS:
%   X = meshgrid for X (100x100 matrix, for example)
%   Y = meshgrid for Y (100x100 matrix, for example)
%   frequencies = a vector of frequencies for the plane waves
%   amplitudes = a vector of amplitudes for the plane waves
%   directions = a vector propagation directions from x-axis (in degrees),
%                   for example: from left to right 0
%                                from down to up 90, etc.
%   attenuation = boolean to switch distance attenuation on/off
%   RETURNS:
%   plane_waves = a vector of complex pressure matrixes

c = 343;
k = 2*pi*frequencies./c;        % wave numbers
ks = zeros(size(directions,2),2); % scaled direction cosines
% Tilts are used to weight x and y directions (to achieve the true
% direction)
tilts = zeros(size(directions,2),2);
plane_waves = zeros(size(X,1),size(X,2),size(directions,2));

for i=1:max(size(tilts))
    % Set the tilts appropriately
    if (abs(directions(i)) > 45)
        tilts(i,1) = directions(i);
    elseif (abs(directions(i)) < 45)
        tilts(i,2) = 90 - abs(directions(i));
    end
    % Scaled direction cosine
    ks(i,1) = k(i)*cos(deg2rad(tilts(i,1)));
    ks(i,2) = k(i)*cos(deg2rad(tilts(i,2)));
    % Plane wave
    plane_waves(:,:,i) = amplitudes(i)*exp(-1i*ks(i,1)*X).*exp(-1i*ks(i,2)*Y);
    % Distance attenuation
    att = sqrt((X*cos(deg2rad(tilts(i,1)))).^2+(Y*cos(deg2rad(tilts(i,2)))).^2) + 1;
    % Flip Y, if the angle of propagation is negative
    if (directions(i) < 0)
        plane_waves(:,:,i) = flip(plane_waves(:,:,i));
        att = flip(att);
    end
    % Add distance attenuation
    if (attenuation)
        plane_waves(:,:,i) = plane_waves(:,:,i)./att;
    end
end