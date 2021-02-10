function [ ds, kf, sf, asD, asA, esD, esA, xpr ] =...
    generate_lsf( tx_pos, rx_pos, lsp_vals, lsp_xcorr, lsf_sos, ES_D_mu_A, ES_D_mu_min )
%GENERATE_LSF Calculates the spatially correlated values of LSPs
%
% Input:
%   tx_pos          Fixed Tx (Single mobility): Tx-position (metric) [ 3 x 1 ]
%                   Dual mobility: Tx-position (metric) [ 3 x N ]
%   rx_pos          Rx-positions (metric) [ 3 x N ]
%   lsp_vals        The distribution values of the LSPs (see below) [ 8 x 9 x F ]
%   lsp_xcorr       The Cross-correlation matrix for the LSPs [ 8 x 8 ]
%   lsf_sos         Pre-initialized SOS generators for the LSF generation [ 1 x 8 ]
%                   If left empty, SOS generator are initialized from the "lsp_vals"
%   ES_D_mu_A       Linear distance-dependence of the ES_D
%   ES_D_mu_min     Minimum value of ES_D_mu if using distance-dependence.
%
% Derived from input:
%   N  = number of users (from rx positions)
%   F  = number of frequencies (from lsp_vals)
%
% Content of variable lsp_vals (first dimension) and the order of the SOS generators:
%   DS  (log10 s)
%   KF  (dB)
%   SF  (dB)
%   ASD (log10 deg)
%   ASA (log10 deg)
%   ESD (log10 deg)
%   ESA (log10 deg)
%   XPR (dB)
%
% Content of variable lsp_vals (second dimension):           Variable    Unit (v = value)
%  1. reference value  with applied freq.-dep. scaling       (mu)        v
%  2. standard deviation with applied freq.-dep. scaling     (sigma)     v
%  3. decorrelation distance in meters                       (lambda)    meters
%  4. TX-RX 2D dist.-dep. of the reference value             (epsilon)   v / log10(m)
%  5. TX height-dep. of the reference value                  (zeta)      v / log10(m)
%  6. elevation-dep. of the reference value                  (alpha)     v / log10(rad)
%  7. TX-RX 2D dist.-dep. of the STD                         (kappa)     v / log10(m)
%  8. TX height-dep. of the STD                              (tau)       v / log10(m)
%  9. elevation-dep. of the STD                              (beta)      v / log10(rad)
%
% Output:
%   ds              Delay spread (s) [ F x N ]
%   kf              K-Factor (linear) [ F x N ]
%   sf              Shadow fading (linear) [ F x N]
%   asD             Azimuth angle spread of departure (deg) [ F x N ]
%   asA             Azimuth angle spread of arrival (deg) [ F x N ]
%   esD             Elevation angle spread of departure (deg) [ F x N ]
%   esA             Elevation angle spread of arrival (deg) [ F x N ]
%   xpr             Per-user XPR (linear) [ F x N ]
%
% 
% QuaDRiGa Copyright (C) 2011-2019
% Fraunhofer-Gesellschaft zur Foerderung der angewandten Forschung e.V. acting on behalf of its
% Fraunhofer Heinrich Hertz Institute, Einsteinufer 37, 10587 Berlin, Germany
% All rights reserved.
%
% e-mail: quadriga@hhi.fraunhofer.de
%
% This file is part of QuaDRiGa.
%
% The Quadriga software is provided by Fraunhofer on behalf of the copyright holders and
% contributors "AS IS" and WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES, including but not limited to
% the implied warranties of merchantability and fitness for a particular purpose.
%
% You can redistribute it and/or modify QuaDRiGa under the terms of the Software License for 
% The QuaDRiGa Channel Model. You should have received a copy of the Software License for The
% QuaDRiGa Channel Model along with QuaDRiGa. If not, see <http://quadriga-channel-model.de/>. 

% Parse input variables
if ~exist('rx_pos','var') || size( rx_pos,1) ~= 3
    error('QuaDRiGa:qd_builder:generate_lsf','Rx position is not given or has invalid format.')
end

% Number of users
N  = size( rx_pos, 2 );
oN = ones(1,N);
o8 = ones(8,1);         % Number of LSPs

if ~exist('tx_pos','var') || size( tx_pos,1) ~= 3
    error('QuaDRiGa:qd_builder:generate_lsf','Tx position is not given or has invalid format.')
end

if size( tx_pos,2 ) == 1
    tx_pos = tx_pos(:,oN);
    tx_pos_SOS = [];
elseif size( tx_pos,2 ) == N
    tx_pos_SOS = tx_pos;
else
    error('QuaDRiGa:qd_builder:generate_lsf','Invalid number of Tx-positions.')
end

if ~exist('lsp_vals','var') || size( lsp_vals,1) ~= 8 || size( lsp_vals,2) ~= 9
    error('QuaDRiGa:qd_builder:generate_lsf','LSP values are not given or have invalid format.')
end

% Throw warning if there is elevation and distance/height dependence at the same time
if ( any( lsp_vals(:,6) ~= 0 ) && ( any( lsp_vals(:,4) ~= 0 ) || any( lsp_vals(:,5) ~= 0 ) ) ) || ...
        ( any( lsp_vals(:,9) ~= 0 ) && ( any( lsp_vals(:,7) ~= 0 ) || any( lsp_vals(:,8) ~= 0 ) ) )
    warning('QuaDRiGa:qd_builder:generate_lsf',...
        'Simultaneous height/distance AND elevation-dependence creates wrong results.');
end

% Number of frequencies
F = size( lsp_vals, 3 );
oF = ones(1,F);

if ~exist('lsp_xcorr','var') || isempty( lsp_xcorr )
    lsp_xcorr = eye( 8 );
elseif size( lsp_xcorr,1) ~= 8 || size( lsp_xcorr,2) ~= 8
    error('QuaDRiGa:qd_builder:generate_lsf','LSP xcorr values are not given or have invalid format.')
end

if ~exist('lsf_sos','var') || isempty( lsf_sos )
    % Initialize SOS generators for the LSF model
    acf = 'Comb300';
    lambda = lsp_vals(:,3,1);
    lsf_sos = qd_sos([]);
    
    % Delays (are identical if Tx and Rx positions are swapped)
    lsf_sos(1,1) = qd_sos( acf, 'Normal', lambda(1) );
    lsf_sos(1,1).sos_phase(:,2) = lsf_sos(1,1).sos_phase(:,1);
    
    % K-Factor (is identical if Tx and Rx positions are swapped)
    lsf_sos(1,2) = qd_sos( acf, 'Normal', lambda(2) );
    lsf_sos(1,2).sos_phase(:,2) = lsf_sos(1,2).sos_phase(:,1);
    
    % Shadow-Fading (is identical if Tx and Rx positions are swapped)
    lsf_sos(1,3) = qd_sos( acf, 'Normal', lambda(3) );
    lsf_sos(1,3).sos_phase(:,2) = lsf_sos(1,3).sos_phase(:,1);
    
    % Azimuth angles (depature become arrival angles if positions are swapped)
    lsf_sos(1,4) = qd_sos( acf, 'Normal', lambda(4) ); % ASD
    lsf_sos(1,5) = qd_sos( acf, 'Normal', lambda(5) ); % ASA
    lsf_sos(1,4).sos_phase(:,1) = lsf_sos(1,5).sos_phase(:,2);
    lsf_sos(1,5).sos_phase(:,1) = lsf_sos(1,4).sos_phase(:,2);
    
    % Elevation angles (depature become arrival angles if positions are swapped)
    lsf_sos(1,6) = qd_sos( acf, 'Normal', lambda(6) ); % ASD
    lsf_sos(1,7) = qd_sos( acf, 'Normal', lambda(7) ); % ASA
    lsf_sos(1,6).sos_phase(:,1) = lsf_sos(1,7).sos_phase(:,2);
    lsf_sos(1,7).sos_phase(:,1) = lsf_sos(1,6).sos_phase(:,2);
    
    % XPR (is identical if Tx and Rx positions are swapped)
    lsf_sos(1,8) = qd_sos( acf, 'Normal', lambda(8) );
    lsf_sos(1,8).sos_phase(:,2) = lsf_sos(1,8).sos_phase(:,1);
    
elseif isa( lsf_sos, 'qd_sos' ) && size(lsf_sos,1) == 1 && size(lsf_sos,2) == 8
    % Test if all generators are set to uniform distribution
    for m = 1 : size( lsf_sos, 2)
        if ~strcmp( lsf_sos(1,m).distribution , 'Normal' )
            warning('QuaDRiGa:qd_builder:generate_lsf:Normal',...
                'All random generators should be set to Normal distribution.');
            warning ('off','QuaDRiGa:qd_builder:generate_lsf:Normal');
        end
    end
    warning ('on','QuaDRiGa:qd_builder:generate_lsf:Normal');
    
    % Test if DS, SF, KF, XPR are reciprocal
    r = [1,2,3,8];
    for n = r
        tmp = cat(3,lsf_sos(:,n).sos_phase);
        tmp = tmp(:,1,:) - tmp(:,2,:);
        if any(abs(tmp(:))>1e-6)
            warning('QuaDRiGa:qd_builder:generate_lsf','DS, SF, KF and XPR must be reciprocal.');
        end
    end
    
    % Test if azimuth spreads are reciprocal
    tmp1 = cat(3,lsf_sos(:,4).sos_phase);
    tmp2 = cat(3,lsf_sos(:,5).sos_phase);
    tmp3 = tmp1(:,1,:) - tmp2(:,2,:);
    tmp4 = tmp2(:,1,:) - tmp1(:,2,:);
    if any( abs( [tmp3(:);tmp4(:)] ) > 1e-6 )
        warning('QuaDRiGa:qd_builder:generate_lsf','Azimuth spreads are not reciprocal.');
    end
    
    % Test if elevation spreads are reciprocal
    tmp1 = cat(3,lsf_sos(:,6).sos_phase);
    tmp2 = cat(3,lsf_sos(:,7).sos_phase);
    tmp3 = tmp1(:,1,:) - tmp2(:,2,:);
    tmp4 = tmp2(:,1,:) - tmp1(:,2,:);
    if any( abs( [tmp3(:);tmp4(:)] ) > 1e-6 )
        warning('QuaDRiGa:qd_builder:generate_lsf','Elevation spreads are not reciprocal.');
    end
else
    error('QuaDRiGa:qd_builder:generate_lsf','lsf_sos has invalid format.');
end

if ~exist('ES_D_mu_A','var') || isempty( ES_D_mu_A )
    ES_D_mu_A = 0;
end

if ~exist('ES_D_mu_min','var') || isempty( ES_D_mu_min )
    ES_D_mu_min = -Inf;
end

% Generate LSF parameters
% Generate 3-D spatially correlated normal distributed random variables
% Combine the 3-D variables to 6-D variables
par = zeros(8,N);
for n = 1 : 8
    par(n,:) = val( lsf_sos(1,n), rx_pos, tx_pos_SOS );
end

% In a dual-mobility scenario, ASD and ASA are correlated when the Tx and Rx are close. This
% correlation hast to be included in the values.
if ~isempty( tx_pos_SOS ) % Dual-mobility
    % The distance between Tx and Rx
    d = sqrt( sum(abs( rx_pos - tx_pos ).^2,1) );
    
    % The distance-dependent correlations vor the four angles
    corr_ASD_ASA = [ lsf_sos(1,4).acfi( d ) ; lsf_sos(1,5).acfi( d ) ];
    corr_ASD_ASA( corr_ASD_ASA > 1-1e-5 ) =  1-1e-5;
    corr_ESD_ESA = [ lsf_sos(1,6).acfi( d ) ; lsf_sos(1,7).acfi( d ) ];
    corr_ESD_ESA( corr_ESD_ESA > 1-1e-5 ) =  1-1e-5;
    
    % Apply the correlation to each user
    for n = 1 : N
        R_sqrt = sqrtm( [ 1 corr_ASD_ASA(1,n) ; corr_ASD_ASA(2,n) 1 ] );
        par(4:5,n)    = R_sqrt * par(4:5,n);
        
        R_sqrt = sqrtm( [ 1 corr_ESD_ESA(1,n) ; corr_ESD_ESA(2,n) 1 ] );
        par(6:7,n)    = R_sqrt * par(6:7,n);
    end
end

% Apply the inter-parameter correlation model
R_sqrt = sqrtm( lsp_xcorr );
par = R_sqrt * par;

% Calculate the variables d2D, hBS and alpha from the TX and RX positions
hBS   = tx_pos(3,:);
hBS( hBS < 1e-4 ) = 1e-4;   % hBS cannot be 0 or smaller than 0
d2D   = sqrt( sum( ( tx_pos(1:2,:) - rx_pos(1:2,:) ).^2,1 ) );
d2D( d2D < 1e-4 ) = 1e-4;   % d2D cannot be 0 or smaller than 0
alpha = atan( hBS ./ d2D );

% Calculate mean values
mu = lsp_vals(:,oN,:) + lsp_vals(:,4*oN,:) .* log10( d2D(o8,:,oF) ) + ...
    lsp_vals(:,5*oN,:) .* log10( hBS(o8,:,oF) ) + lsp_vals(:,6*oN,:) .* log10( alpha(o8,:,oF) );

% Apply distant-dependent ESD scaling with linear distance
if ES_D_mu_A ~= 0 
    mu(6,:,:) = mu(6,:,:) + ES_D_mu_A .* d2D(1,:,oF)/1000;
end 

% Apply mimimum value for the ESD
mu( 6, mu(6,:) < ES_D_mu_min ) = ES_D_mu_min;

% Calculate STD
sigma = lsp_vals(:,2*oN,:) + lsp_vals(:,7*oN,:) .* log10( d2D(o8,:,oF) ) + ...
    lsp_vals(:,8*oN,:) .* log10( hBS(o8,:,oF) ) + lsp_vals(:,9*oN,:) .* log10( alpha(o8,:,oF) );

% Apply mu and sigma from the parameter table for each frequency
par = par( :,:,oF ) .* sigma + mu;

% Transform to linear values
par( [1,4:7],:,: ) = 10.^( par( [1,4:7],:,: ) );
par( [2,3,8],:,: ) = 10.^( 0.1 * par( [2,3,8],:,: ) );

% Generate output values
ds   = permute( par(1,:,:) , [3,2,1] );
kf   = permute( par(2,:,:) , [3,2,1] );
sf   = permute( par(3,:,:) , [3,2,1] );
asD  = permute( par(4,:,:) , [3,2,1] );
asA  = permute( par(5,:,:) , [3,2,1] );
esD  = permute( par(6,:,:) , [3,2,1] );
esA  = permute( par(7,:,:) , [3,2,1] );
xpr  = permute( par(8,:,:) , [3,2,1] );

end
