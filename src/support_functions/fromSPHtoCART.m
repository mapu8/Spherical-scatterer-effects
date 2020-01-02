function [M] = fromSPHtoCART(X, micDistances, H_single_freq, mic_azi, doa_azi)
%FROMSPHTOCART Transform the measured points of sphericalScatterer() from
%              spherical coordinates to cartesian coordinates.
%   PARAMETERS:
%   X = meshgrid for X (100x100 matrix, for example)
%   micDistances = a vector of mic distance array
%   H_single_freq = simulated impulse response matrix for single frequency,
%                   obtained from sphericalScatterer()
%   mic_azi = azimuth angle of the microphone
%   doa_azi = a vector of azimuth angles of the DOA of the plane waves
%
%   RETURNS:
%   M = a matrix of pressure values in cartesian coordinates

M = zeros(size(X,1));
midx = size(M,1)/2;
midy = midx;
trueLength = X(1,2);  % The lenght between matrix elements in meters.
micIndexes = round((micDistances)./trueLength);   % Indexes of the mic positions in the matrix array.
% Loop over directions around the sphere
for dir=1:size(H_single_freq,2)
   % Spherical coordinates
   phi = mic_azi + doa_azi(dir);
   theta = 90;
   xmult = sin(theta)*cos(phi);
   ymult = sin(theta)*sin(phi);
   % Transform to cartesian coordinates
   x_indexes = round(xmult*micIndexes) + midx;
   y_indexes = round(ymult*micIndexes) + midy;
   for i=1:size(x_indexes,1)
      M(x_indexes(i),y_indexes(i)) = H_single_freq(i,dir);
   end
end
end

