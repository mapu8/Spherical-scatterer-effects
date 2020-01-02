function [D, order] = ambiDecoder_amplitude_norm(ls_dirs, method, rE_WEIGHT, order, MOA_2D_order, horizontal_response)
%AMBIDECODER_AMPLITUDE_NORM Returns a HOA decoding matrix for a loudspeaker setup.
%NOTE: this script uses amplitude normalization instead of energy
%normalization for vbap
% AMBIDECODER computed a HOA decoding matrix for a loudspeaker setup, based
% on one of five implemented decoder designs. The methods are:
%   - Sampling decoder (SAD)
%   - Mode-matching decoder (MMD)
%   - Energy-preserving decoder (EPAD)
%   - All-round ambisonic panning (ALLRAD)
%   - Constant Angular Spread Decoder (CSAD)
%
%   The two first ones are traditional design approaches. The last three
%   ones are recently proposed design methods that are more flexible and
%   more psychoacoustically motivated. For references, check the
%   TEST_AMBI_SCRIPT.m. Additionally, the so-called max-rE weighting can be
%   enabled for all of the above decoders apart from CSAD. The max-rE
%   weighting maximises the norm of the energy vector for all decoding
%   directions, which in ambisonic literature is considered to reduce 
%   localization blur.
%
% Inputs:   
%   ls_dirs: speaker directions in [azi1 elev1; azi2 elev2;... ; aziL elevL]
%            convention, in degrees, for L loudspeakers in the layout
%   method:  {'SAD','MMD','EPAD','ALLRAD','CSAD'} for any one of the 
%            respective decoding methods
%   rE_WEIGHT:  Enable max-rE weighting. If CSAD is chosen, then this does
%               nothing.
%   order:      Ambisonic order of the array. If it is left blank, then an
%               equivalent ambisonic order is computed by
%               getLayoutAmbisonicOrder().
%   MOA_2D_order: upper limit of Mixed Order Ambisonics order, calculate
%                 only the degree=+-order SH functions for
%                 MOA_2D_order>MOA_3D_order
%   horizontal_response: {'SH','CH'}, spherical or cylindrical horizontal
%                        functions
%
% Outputs:
%   D:      The [L x (order+1)^2] ambisonic decoding matrix
%   order:  The order of the layout, useful to see what it is if it is not 
%           defined and it is computed internally.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Archontis Politis, 15/11/2015
%   archontis.politis@aalto.fi
%
%   Modified by Lauros Pajunen to support Mixed Order Ambisonics, 20/11/2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find speaker triplets
ls_num = size(ls_dirs,1);

if ~exist('order','var')
    % compute ambisonic equivalent order (from Zotter & Frank)
    order = getLayoutAmbisonicOrder(ls_dirs);
end

if ~exist('rE_WEIGHT','var')
    rE_WEIGHT = 1;
end

if rE_WEIGHT
    a_n = getMaxREweights(order);
    if (exist('MOA_2D_order'))
        if (exist('horizontal_response'))
            switch lower(horizontal_response)

                % 2D only
                case {'sh','ch'}
                    a_n = getMaxREweights(MOA_2D_order);
                    a_temp = zeros(2*MOA_2D_order + 1,1);
                    a_temp(1,:) = a_n(1,:);
                    index = 2;
                    for n = 0:MOA_2D_order-1
                        a_temp(index,:) = a_n((n+1)^2 + 1,:);
                        a_temp(index+1,:) = a_n((n+2)^2,:);
                        index = index + 2;
                    end
                    a_n = a_temp;
            end
        % 3D and 2D
        else
            a_n = getMaxREweights(MOA_2D_order);
            a_temp = zeros((order+1)^2 + 2*(MOA_2D_order-order),1);
            a_temp(1:(order+1)^2,1) = a_n(1:(order+1)^2,1);
            index = (order+1)^2 + 1;
            for n = order:MOA_2D_order-1
                a_temp(index,1) = a_n((n+1)^2 + 1);
                a_temp(index+1,1) = a_n((n+1)^2 + 1 + 2*(n+1));
                index = index + 2;
            end
            a_n = a_temp;
        end
    end
end

% compute real SH matrix for layout
Y_ls = getRSH(order, ls_dirs);
if (exist('MOA_2D_order'))
    if (exist('horizontal_response'))
        switch lower(horizontal_response)

            % Spherical harmonics in 2D only
            case 'sh'
                Y_ls = getRSH(MOA_2D_order, ls_dirs).';
                Y_temp = zeros(size(Y_ls,1), 2*MOA_2D_order + 1);
                Y_temp(:,1) = Y_ls(:,1);
                index = 2;
                for n = 0:MOA_2D_order-1
                    Y_temp(:,index) = Y_ls(:,(n+1)^2 + 1) /(4*pi) *(2*pi);
                    Y_temp(:,index+1) = Y_ls(:,(n+2)^2) /(4*pi) *(2*pi);
                    index = index + 2;
                end
                Y_ls = Y_temp.';
                
            % Cylindrical harmonics
            case 'ch'
                ls_dirs_2D = ls_dirs(ls_dirs(:,2) == 0,:);
                ls_num = size(ls_dirs_2D,1);
                Y_ls = getRCH_mod(MOA_2D_order, ls_dirs_2D(:,1)*pi/180).';
                %Y_temp = zeros(size(Y_ls,1), 2*MOA_2D_order + 1);
                %Y_temp(:,1) = Y_ls(:,1);
                %index = 2;
                %for n = 0:MOA_2D_order-1
                %    Y_temp(:,index) = Y_ls(:,(n+1)^2 + 1);
                %    Y_temp(:,index+1) = Y_ls(:,(n+2)^2);
                %    index = index + 2;
                %end
                %Y_ls = Y_temp;
        end
    % SH functions in 3D up to 'order' and 2D SH functions (m=+-n) for
    % 'order' < order < MOA_2D_order
    else
        Y_ls = getRSH(MOA_2D_order, ls_dirs).';
        Y_temp = zeros(size(Y_ls,1),(order+1)^2 + 2*(MOA_2D_order-order));
        Y_temp(:,1:(order+1)^2) = Y_ls(:,1:(order+1)^2);
        index = (order+1)^2 + 1;
        for n = order:MOA_2D_order-1
            Y_temp(:,index) = Y_ls(:,(n+1)^2 + 1) /(4*pi) *(2*pi);
            Y_temp(:,index+1) = Y_ls(:,(n+2)^2) /(4*pi) *(2*pi);
            index = index + 2;
        end
        Y_ls = Y_temp.';
    end
end

switch lower(method)
    case 'sad'
        D = (4*pi)/ls_num * Y_ls.';
    case 'mmd'
        D = pinv(Y_ls);
    case 'epad'
        [U,S,V] = svd(Y_ls);
        S_trunc = S(1:(order+1)^2,1:(order+1)^2);
        V_trunc = V(:,1:(order+1)^2);
        D = (4*pi)/ls_num * V_trunc*U.';
    case 'allrad'
        %D = allrad(ls_dirs, order);       
        if (exist('MOA_2D_order'))
            if (exist('horizontal_response'))
                D = allrad_amplitude_norm(ls_dirs, order, MOA_2D_order, horizontal_response); 
            else
                D = allrad_amplitude_norm(ls_dirs, order, MOA_2D_order);  
            end
        else
            D = allrad_amplitude_norm(ls_dirs, order);    
        end
    case 'allrad_normal'
        D = allrad(ls_dirs, order);   
    case 'csad'
        D = csad(ls_dirs, order);
end

% apply rE weights if defined
if rE_WEIGHT && ~isequal(method,'csad')
    D = D*diag(a_n);
end

end
