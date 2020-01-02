function [a_nm] = encodeWithDTUMicrophone(doa_vec,regularization,fs,Lfilt,maxG_dB,MOA_response)
%ENCODEWITHDTUMICROPHONE Capture plane waves with a virtual DTU microphone, and encode
%                        the recorded signals to spherical harmonics up to 5th order in 3D
%                        and 7th order in Mixed Order Ambisonics (MOA) that utilizes the
%                        horizontal microphones for higher order response.
%
%   Simulate 52-channel spherical microphone developed in DTU [ref2] recordings with simulateSphArray(), and encode the
%   recording to 5th (3D) or 7th order (MOA) spherical harmonics (SH). The encoding to SH weights
%   a_nm is done as in [ref1,ref5]:
%
%       a_nm = diag(H_filt)*pseudo_inverse(Y_mics)*mic_signals.
%
%   Y_mics holds the SH functions for the microphone array, and mic_signals
%   are the simulated pressure signals obtained with simulateSphArray().
%   H_filt is the regularization filter, which can be chosen as one of the
%   following listed below.
%
%   Regularization methods:
%   - Regularized radial inverse (RADINV) [ref1], default
%   - Regularized least-squares (REGLS) [ref4]
%   - Soft-limiter (SOFTLIM) [ref3]
%   - No regularization (NONE)
%
%   Necessary external libraries can be found in [ref6-9]

% ARGUMENTS:
% doa_vec - direction of arrival vector for all desired directions,
%           [azi1 ele1; azi2 ele2; ...; azin elen] in radians
% regularization - {'RADINV_DTU','RADINV_DTU_MOA'} for a regularization filter
% fs - sampling frequency
% Lfilt - number of FFT points up to fs
% maxG_dB - maximum allowed amplification for the regularization filters
% MOA_response - set to 0, to output up to 5th order 3D SH 
%                set to 1, to output up to 5th 3D and 7th order horizontal SH weights
%                set to 2, to output up to 7th order horizontal only SH weights
%                set to 3, to output up to 7th order CH weights

% OUTPUTS:
% a_nm - spherical harmonic weights up to order N=5 (3D) or N=5/7 (MOA) and 
%        frequency fs/2 for each DOA


% Helper function to convert directions from Matlab's azimuth-elevation to
% azimuth-inclination
aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

% Microphone position angles in radians [azimuth, elevation]
mic_dirs_rad = ...
    [0.0                   1.40194351305251;
    -3.14159265358979      1.40194351305251;
    0.523598775598299      0.964876234040517;
    1.57079632679490       0.964876234040517;
    2.61799387799149       0.964876234040517;
    -2.61799387799150      0.964876234040517;
    -1.57079632679490      0.964876234040517;
    -0.523598775598299      0.964876234040517;
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
mic_dirs_incl = aziElev2aziIncl(mic_dirs_rad); % radians [azi incl]
nMics = size(mic_dirs_rad,1);

% Variables
c = 343;
R = 0.05; % radius of mic
nBins = Lfilt/2 + 1;
f_max = fs/2;
kR_max = 2*pi*f_max*R/c;
array_order = ceil(2*kR_max);
f = (0:Lfilt/2)'*fs/Lfilt;
kR = 2*pi*f*R/c;
arrayType = 'rigid';
Ndoas = size(doa_vec,1);

% Spherical harmonics for the mic array
% MOA response with both 3D and 2D spherical harmonic functions
if (MOA_response == 1 || MOA_response == 2)
    sht_order = 7;
    Y_mics = getSH(sht_order, mic_dirs_incl, 'real'); % real SH matrix for microphones
% 2D response with cylindrical harmonics (CH)
elseif (MOA_response == 3)
    
    warning("2D cylindrical harmonics are not implemented in this version of the code, aborting.")
    return
    
    sht_order = 7;
    mic_dirs_rad = mic_dirs_rad(mic_dirs_rad(:,2) == 0,:);
    mic_dirs_incl = aziElev2aziIncl(mic_dirs_rad); % radians [azi incl]
    nMics = size(mic_dirs_rad,1);
    Y_mics = getRCH_mod(sht_order, mic_dirs_incl(:,1)); % real CH matrix for microphones
    % Expand to Qx(sht_order+1)^2 matrix for calculations
    Y_temp = zeros(size(Y_mics,1), (sht_order+1)^2);
    Y_temp(:,1) = Y_mics(:,1);
    index = 2;
    for n = 0:sht_order-1
        Y_temp(:,(n+1)^2 + 1) = Y_mics(:,index);
        Y_temp(:,(n+2)^2) = Y_mics(:,index+1);
        index = index + 2;
    end
    Y_mics = Y_temp;
% 3D response
else
    sht_order = 5;
    Y_mics = getSH(sht_order, mic_dirs_incl, 'real'); % real SH matrix for microphones
end


% Pseudo-inversion of the microphone SH functions
pseudo_Y_mics = pinv(Y_mics);

% Simulate the response from the DTU microphone
[~, H_array_sim] = simulateSphArray(Lfilt, mic_dirs_rad, doa_vec, arrayType, R, array_order, fs);

% Sound field coefficients on the surface of a rigid sphere
b_N = sphModalCoeffs(sht_order, kR, arrayType, []);

% Regularization filter
if (exist('regularization'))
    switch (lower(regularization))
        case {'radinv','radinv_dtu'}
            [H_filt, ~] = arraySHTfiltersTheory_radInverse(R, nMics, sht_order, Lfilt, fs, maxG_dB);
            
        case {'radinv_dtu_moa','radinv_dtu_moa_2d_sh','radinv_dtu_moa_2d_ch'}
            [H_filt, ~] = arraySHTfiltersTheory_radInverse2D(R, nMics, sht_order, Lfilt, fs, maxG_dB);
        
        case 'regls'
            [H_filt, ~] = arraySHTfiltersTheory_regLS(R, mic_dirs_rad, sht_order, Lfilt, fs, maxG_dB);

        case 'softlim'
            [H_filt, ~] = arraySHTfiltersTheory_softLim(R, nMics, sht_order, Lfilt, fs, maxG_dB);
        
        case 'none'
            H_filt = 1./b_N;

        otherwise
            warning("Couldn't find regularization method: %s.\n Using radial inverse as default.",regularization);
            [H_filt, ~] = arraySHTfiltersTheory_radInverse(R, nMics, sht_order, Lfilt, fs, maxG_dB);
    end
else
    warning("Regularization method was not set, using radial inverse as default.")
    [H_filt, ~] = arraySHTfiltersTheory_radInverse(R, nMics, sht_order, Lfilt, fs, maxG_dB);
end

% Calculate the SH weights: a_nm = diag(H_filt)*pseudo_inverse(Y_mics)*mic_signals
a_nm = zeros(size(Y_mics,2),size(H_array_sim,1),size(H_array_sim,3));
if (exist('regularization') && strcmpi(regularization,'regls'))
    for doa_indx = 1:Ndoas
        for f_indx = 1:nBins
            a_nm(:,f_indx,doa_indx) = H_filt(:,:,f_indx)*H_array_sim(f_indx,:,doa_indx).';
        end
    end
else
    for doa_indx = 1:Ndoas
        for f_indx = 1:nBins
            a_nm(:,f_indx,doa_indx) = diag(replicatePerOrder(H_filt(f_indx,:),2))*pseudo_Y_mics*H_array_sim(f_indx,:,doa_indx).';
        end
    end
end

% Select all the SH weights up to 5th order. Above that, select only the
% weights with degree = +-order (m = -6,6,-7,7).
if (MOA_response == 1)
    MOA_lower_order = 5;
    a_nm_temp = zeros((MOA_lower_order+1)^2 + 2*(sht_order-MOA_lower_order),size(a_nm,2),size(a_nm,3));
    a_nm_temp(1:(MOA_lower_order+1)^2,:,:) = a_nm(1:(MOA_lower_order+1)^2,:,:);
    index = (MOA_lower_order+1)^2 + 1;
    for n = MOA_lower_order:sht_order-1
        a_nm_temp(index,:,:) = a_nm((n+1)^2 + 1,:,:);% /(4*pi) *(2*pi);
        a_nm_temp(index+1,:,:) = a_nm((n+2)^2,:,:);% /(4*pi) *(2*pi);
        index = index + 2;
    end
    a_nm = a_nm_temp;
% Only 2D SH functions or CH functions
elseif (MOA_response == 2 || MOA_response == 3)
    a_nm_temp = zeros(2*sht_order+1,size(a_nm,2),size(a_nm,3));
    a_nm_temp(1,:,:) = a_nm(1,:,:);
    index = 2;
    for n = 0:sht_order-1
        a_nm_temp(index,:,:) = a_nm((n+1)^2 + 1,:,:);% /(4*pi) *(2*pi);
        a_nm_temp(index+1,:,:) = a_nm((n+2)^2,:,:);% /(4*pi) *(2*pi);
        index = index + 2;
    end
    a_nm = a_nm_temp;
end

% REFERENCES
%
%   1.  Moreau, S., Daniel, J., Bertet, S., 2006, 
%       3D sound field recording with higher order ambisonics-objective measurements and validation of spherical microphone. 
%       In Audio Engineering Society Convention 120.
%
%   2.  Marcschall, M., Favrot, S., Buchholz, J., 2012,
%       Robustness of a mixed-order Ambisonics microphone array for sound field reproduction.
%       In Audio Engineering Society Convention 132.
%
%   3.  Bernsch?tz, B., P?rschmann, C., Spors, S., Weinzierl, S., Verst?rkung, B., 2011. 
%       Soft-limiting der modalen amplitudenverst?rkung bei sph?rischen mikrofonarrays im plane wave decomposition verfahren. 
%       Proceedings of the 37. Deutsche Jahrestagung f?r Akustik (DAGA 2011)
%
%   4.  Jin, C.T., Epain, N. and Parthy, A., 2014. 
%       Design, optimization and evaluation of a dual-radius spherical microphone array. 
%       IEEE/ACM Transactions on Audio, Speech, and Language Processing, 22(1), pp.193-204.
%
%   5.  Rafaely, B., 2015.
%       Fundamentals of Spherical Array Processing,
%       volume 8 of Springer Topics in Signal Processing. Springer-Verlag
%       Berlin Heidelberg.
%       193 p. ISBN 978-3-662-45663-7.
%
%   6.  Politis, A., 2014.
%       Spherical Array Processing. Matlab library. https://github.com/polarch/Spherical-Array-Processing
%   
%   7.  Politis, A, 2015.
%       Spherical Harmonic Transform. Matlab library. https://github.com/polarch/Spherical-Harmonic-Transform
%   
%   8.  Politis, A., 2015.
%       Higher Order Ambisonics. Matlab library. https://github.com/polarch/Higher-Order-Ambisonics
%
%   9.  Politis, A., 2015.
%       Array Response Simulator. Matlab library. https://github.com/polarch/Array-Response-Simulator
%

end

