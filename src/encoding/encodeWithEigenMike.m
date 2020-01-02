function [a_nm] = encodeWithEigenMike(doa_vec,regularization,fs,Lfilt,maxG_dB)
%ENCODEWITHEIGENMIKE Capture plane waves with a virtual EigenMike, and encode
%                    the recorded signals to spherical harmonics up to 4th order.
%
%   Simulate EigenMike [ref2] recordings with simulateSphArray(), and encode the
%   recording to 4th order spherical harmonics (SH). The encoding to SH weights
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
% regularization - {'RADINV','REGLS','SOFTLIM','NONE'} for a
%                  regularization filter
% fs - sampling frequency
% Lfilt - number of FFT points up to fs
% maxG_dB - maximum allowed amplification for the regularization filters

% OUTPUTS:
% a_nm - spherical harmonic weights up to order N=4 and frequency fs/2
%        for each DOA


% Helper function to convert directions from Matlab's azimuth-elevation to
% azimuth-inclination
aziElev2aziIncl = @(dirs) [dirs(:,1) pi/2-dirs(:,2)];

% Eigenmike angles in degrees [azimuth, elevation]
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
mic_dirs_rad = mic_dirs*pi/180; % radians [azi elev]
mic_dirs_incl = aziElev2aziIncl(mic_dirs_rad); % radians [azi incl]
nMics = size(mic_dirs,1);

% Variables
c = 343;
R = 0.042; % radius of EigenMike
nBins = Lfilt/2 + 1;
f_max = fs/2;
kR_max = 2*pi*f_max*R/c;
array_order = ceil(2*kR_max);
f = (0:Lfilt/2)'*fs/Lfilt;
kR = 2*pi*f*R/c;
arrayType = 'rigid';
Ndoas = size(doa_vec,1);

% Spherical harmonics for the mic array
sht_order = floor(sqrt(nMics)-1); % approximate for uniformly arranged mics
Y_mics = getSH(sht_order, mic_dirs_incl, 'real'); % real SH matrix for microphones

% Pseudo-inversion of the microphone SH functions
pseudo_Y_mics = pinv(Y_mics);

% Simulate the response from EigenMike
[~, H_array_sim] = simulateSphArray(Lfilt, mic_dirs_rad, doa_vec, arrayType, R, array_order, fs);

% Sound field coefficients on the surface of a rigid sphere
b_N = sphModalCoeffs(sht_order, kR, arrayType, []);

% Regularization filter
if (exist('regularization'))
    switch (lower(regularization))
        case 'radinv'
            [H_filt, ~] = arraySHTfiltersTheory_radInverse(R, nMics, sht_order, Lfilt, fs, maxG_dB);
        
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

% REFERENCES
%
%   1.  Moreau, S., Daniel, J., Bertet, S., 2006, 
%       3D sound field recording with higher order ambisonics-objective measurements and validation of spherical microphone. 
%       In Audio Engineering Society Convention 120.
%
%   2.  Mh Acoustics Eigenmike, https://mhacoustics.com/products#eigenmike1
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

