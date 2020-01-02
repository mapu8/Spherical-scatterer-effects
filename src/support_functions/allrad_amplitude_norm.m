function D_allrad = allrad_amplitude_norm(ls_dirs, order, MOA_2D_order, horizontal_response)
%ALLRAD_AMPLITUDE_NORM Implements the All-round Ambisonic Decoding of Zotter & Frank.
%NOTE: this script uses amplitude normalization instead of energy
%normalization for vbap
% ALLRAD implements the all-round ambisonic decoding published by Zotter &
% Frank in
%
%   Zotter, F., Frank, M. (2012). All-Round Ambisonic Panning and Decoding. 
%   Journal of the Audio Engineering Society, 60(10), 807:820.
%
% It is an ambisonic decoder which combines rendering to an ideal uniform
% virtual loudspeaker layout with energy preserving properties, then
% rendered to an arbitrary loudspeaker setup by means of vector-base
% amplitude panning (VBAP). The code requires the VBAP library and the
% Spherical Harmonic transform library (for VBAP and t-Designs accordingly)
% found in:
%
% <https://github.com/polarch/Vector-Base-Amplitude-Panning>
% <https://github.com/polarch/Spherical-Harmonic-Transform>
%
% Inputs:   
%   ls_dirs: speaker directions in [azi1 elev1; azi2 elev2;... ; aziL elevL]
%            convention, in degrees
%   order:   order of the HOA decoding matrix. For an irregular speaker
%            layout, a well behaved decoding order can be found by the
%            getLayoutAmbisonicOrder() function.
%   MOA_2D_order: upper limit of Mixed Order Ambisonics order, calculate
%                 only the degree=+-order SH functions for
%                 MOA_2D_order > order
%   horizontal_response: {'SH','CH'}, spherical or cylindrical horizontal
%                        functions
%
% Outputs:
%   D_allrad: [L x (order+1)^2] decoding matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%   Modified by Lauros Pajunen to support Mixed Order Ambisonics, 20/11/2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % t-value for the t-design
        %t = 2*order + 1;
        t = 20;
        % vbap gains for selected t-design
        [~, t_dirs_rad] = getTdesign(t);
        t_dirs = t_dirs_rad*180/pi;
        %G_td = vbap(t_dirs , findLsTriplets(ls_dirs), invertLsMtx(ls_dirs, findLsTriplets(ls_dirs))).';
        G_td = vbap_amplitude_norm(t_dirs , findLsTriplets(ls_dirs), invertLsMtx(ls_dirs, findLsTriplets(ls_dirs))).';
        
        % spherical harmonic matrix for t-design
        % convert to [azimuth zenith] for SH convention
        Y_td = getRSH(order, t_dirs).';
        if (exist('MOA_2D_order'))
            if (exist('horizontal_response'))
                switch lower(horizontal_response)
                    
                    % Spherical harmonics in 2D only
                    case 'sh'
                        Y_td = getRSH(MOA_2D_order, t_dirs).';
                        Y_temp = zeros(size(Y_td,1), 2*MOA_2D_order + 1);
                        Y_temp(:,1) = Y_td(:,1);
                        index = 2;
                        for n = 0:MOA_2D_order-1
                            Y_temp(:,index) = Y_td(:,(n+1)^2 + 1);
                            Y_temp(:,index+1) = Y_td(:,(n+2)^2);
                            index = index + 2;
                        end
                        Y_td = Y_temp;
                        
                    % Cylindrical harmonics only (NOT working yet)
                    case 'ch'
                        Y_td = getRCH_mod(MOA_2D_order, t_dirs_rad(:,1));
                        %Y_temp = zeros(size(Y_td,1),2*MOA_2D_order + 1);
                        %Y_temp(:,1) = Y_td(:,1);
                        %index = (order+1)^2 + 1;
                        %for n = order:MOA_2D_order-1
                        %    Y_temp(:,index) = Y_td(:,(n+1)^2 + 1);
                        %    Y_temp(:,index+1) = Y_td(:,(n+2)^2);
                        %    index = index + 2;
                        %end
                        %Y_td = Y_temp;
                        ls_dirs_2D = ls_dirs(ls_dirs(:,2) == 0,:);
                        G_td = vbap_amplitude_norm(t_dirs(:,1) , findLsPairs(ls_dirs_2D(:,1)), invertLsMtx(ls_dirs_2D(:,1), findLsPairs(ls_dirs_2D(:,1)))).';

                    otherwise
                        fprintf("Couldn't find any case for setting: %s\n",lower(horizontal_response));
                end
            % Full 3D Spherical harmonics up to 'order', and 2D SH for
            % the orders above until MOA_2D_order
            else
                Y_td = getRSH(MOA_2D_order, t_dirs).';
                Y_temp = zeros(size(Y_td,1),(order+1)^2 + 2*(MOA_2D_order-order));
                Y_temp(:,1:(order+1)^2) = Y_td(:,1:(order+1)^2);
                index = (order+1)^2 + 1;
                for n = order:MOA_2D_order-1
                    Y_temp(:,index) = Y_td(:,(n+1)^2 + 1);% /(4*pi) *(2*pi);
                    Y_temp(:,index+1) = Y_td(:,(n+2)^2);% /(4*pi) *(2*pi);
                    index = index + 2;
                end
                Y_td = Y_temp;
            end
        end
        
        % allrad decoder
        Ntd = size(t_dirs_rad,1);
        D_allrad = 4*pi/Ntd * G_td * Y_td;
       
end

