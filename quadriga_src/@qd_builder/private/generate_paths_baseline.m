function [ pow, taus, AoD, AoA, EoD, EoA ] = generate_paths_baseline( L, tx_pos, rx_pos, spreads, r_DS, LNS_ksi, alt_AS_mapping )
%GENERATE_PATHS Generate multipath components according to the 3GPP 38.901 baseline model
%
% Input:
%   L               Number of paths (including LOS) [ 1 x 1 ]
%   tx_pos          Fixed Tx (Single mobility): Tx-position (metric) [ 3 x 1 ]
%   rx_pos          Rx-positions (metric) [ 3 x N ]
%   spreads         Vector of delay and angular spreads (see below) [ 6 x N ]
%   r_DS            Delay distribution proportionality factor (must be > 1); [ 1 x 1 ]
%   LNS_ksi         Per cluster shadowing std / [dB]
%   alt_AS_mapping  Option 1: Random angles scaled by angular spread (correct calibration)
%                   Option 0: Mapping as described in the 3GPP TR (incorrect calibration)
%
%   L  = number of paths
%   Ln = number of NLOS paths ( L-2 for ground-reflection, L-1 otherwise )
%   N  = number of users (from rx positions)
%
% Content of variable spreads:
%   DS (s)
%   ASD (deg)
%   ASA (deg)
%   ESD (deg)
%   ESA (deg)
%   KF (linear)
%
% Output:
%   pow             Path powers (linear, normalized to sum 1) [ N x L ]
%   taus            Path delays (seconds) [ N x L ]
%   AoD             Azimuth angles of departure (rad) [ N x L ]
%   AoA             Azimuth angles of arrival (rad) [ N x L ]
%   EoD             Elevation angles of departure (rad) [ N x L ]
%   EoA             Elevation angles of arrival (rad) [ N x L ]
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
speed_of_light = qd_simulation_parameters.speed_of_light;

if ~exist('L','var') || isempty( L )
    error('QuaDRiGa:qd_builder:generate_paths','Numer of paths must be given.')
end

if ~exist('rx_pos','var') || size( rx_pos,1) ~= 3
    error('QuaDRiGa:qd_builder:generate_paths','Rx position is not given or has invalid format.')
end

% Number of users
N = size( rx_pos, 2 );
oN = ones(1,N);

if ~exist('tx_pos','var') || size( tx_pos,1) ~= 3
    error('QuaDRiGa:qd_builder:generate_paths','Tx position is not given or has invalid format.')
end

if size( tx_pos,2 ) == 1
    tx_pos = tx_pos(:,oN);
    tx_pos_SOS = [];
elseif size( tx_pos,2 ) == N
    tx_pos_SOS = tx_pos;
else
    error('QuaDRiGa:qd_builder:generate_paths','Invalid number of Tx-positions.')
end

% distances between BS and MT
d_2d = hypot( tx_pos(1,:) - rx_pos(1,:), tx_pos(2,:) - rx_pos(2,:) );
d_3d = sqrt( sum((rx_pos - tx_pos).^2,1) );
d_2d( d_2d<1e-5 ) = 1e-5;
d_3d( d_3d<1e-5 ) = 1e-5;

% Calculate angles between BS and MT
angles = zeros( 5,N );
angles(1,:) = atan2( rx_pos(2,:) - tx_pos(2,:) , rx_pos(1,:) - tx_pos(1,:) );   % Azimuth at BS
angles(2,:) = pi + angles(1,:);                                                 % Azimuth at MT
angles(3,:) = atan2( ( rx_pos(3,:) - tx_pos(3,:) ), d_2d );                     % Elevation at BS
angles(4,:) = -angles(3,:);                                                     % Elevation at MT
angles(5,:) = -atan2( ( rx_pos(3,:) + tx_pos(3,:) ), d_2d );                    % Ground Reflection Elevation at BS and MT
angles = angles.';

if ~exist('spreads','var') || size( spreads,1) ~= 6 || size( spreads,2) ~= N
    error('QuaDRiGa:qd_builder:generate_paths','Tx position is not given or has invalid format.')
end

% Number of frequencies
F = size( spreads, 3 );
if F > 1
    error('QuaDRiGa:qd_builder:generate_paths','Multi-frequency simulations are not supported by the 3GPP baseline SSF model.')
end

% Split spread variables into individual variables and convert from [deg] to [rad]
spreadsP = permute( spreads , [2,3,1] );
DS  = spreadsP(:,:,1);
ASD = spreadsP(:,:,2) * pi/180;
ASA = spreadsP(:,:,3) * pi/180;
ESD = spreadsP(:,:,4) * pi/180;
ESA = spreadsP(:,:,5) * pi/180;
KF  = spreadsP(:,:,6);

% Get the number of deterministic and non-deterministic paths
Lf = 1;
Ln = L-1;
oLn = ones(1,Ln);
oLf = ones(1,Lf);
oL  = ones(1,L);

if ~exist('r_DS','var') || isempty( r_DS )
    r_DS = 2.2;
end

if ~exist('LNS_ksi','var') || isempty( LNS_ksi )
    LNS_ksi = 3;
end

% Generate paths
if L == 1   % Only LOS component is present
    
    pow  = ones(N,1,F);
    taus = zeros(N,1);
    AoD  = angles(:,1);
    AoA  = angles(:,2);
    EoD  = angles(:,3);
    EoA  = angles(:,4);
      
else % Additional NLOS componenets are present
    
    % Step 5: Generate cluster delays
    % Generate initial delays following an exponenetial distribution; eq. (7.5-1)
    X = rand( N,L );
    taus = -r_DS * DS(:,oL) .* log( X );
    
    % Subtract the minimum delay and sort the normalised delays to ascending order; eq. (7.5-2)
    taus = sort( taus - min(taus,[],2)*oL , 2 );
    
    % aAditional scaling of delays is required to compensate for the effect of LOS peak; eq. (7.5-3)
    KF_dB = 10*log10( KF );
    C_tau = 0.7705 - 0.0433*KF_dB + 0.0002*KF_dB.^2 + 0.000017*KF_dB.^3;
    C_tau( C_tau > 1 ) = 1;
    
    % Scaled delays (not to be used in cluster power generation); eq. (7.5-4)
    taus_LOS = taus ./ C_tau(:,oL);
    
    % Step 6: Generate cluster powers
    % Cluster powers are calculated assuming a single slope exponential power delay profile; eq. (7.5-5)
    Z = randn( N,L )*LNS_ksi;
    pow = exp( -taus .* (((r_DS-1)./( r_DS*DS ))*oL) ) .* 10.^(Z/10);
    
    % Power of the single LOS ray; eq. (7.5-7)
    pow(:,1) = KF./(KF+1);
    
    % Normalize the cluster powers so that the sum of all cluster powers is equal to one; eq. (7.5-8)
    pow(:,2:end) = (1./(KF+1))*oLn .* pow(:,2:end) ./ (sum( pow(:,2:end),2 )*oLn);
     
    % Step 7: Generate arrival angles and departure angles for both azimuth and elevation.
    % C_phi_NLOS is defined as a scaling factor related to the total number of clusters.
    % The "-1" accounts for the LOS component which is always present in QuaDRiGa.
    if L <= 5
        C_phi_NLOS = 0.779;
    elseif L >= 21
        C_phi_NLOS = 1.289;
    else
        L_phi_NLOS = [     4,     5,     8,    10,    11,    12,    14,    15,    16,    19,    20 ];
        C_phi_NLOS = [ 0.779, 0.860, 1.018, 1.090, 1.123, 1.146, 1.190, 1.211, 1.226, 1.273, 1.289 ];
        C_phi_NLOS = qf.interp( L_phi_NLOS, 0, C_phi_NLOS, L-1 );
    end
    
    % C_phi_LOS depends on the K-Factor, eq. (7.5-10)
    C_phi_LOS = 1.1035 - 0.028*KF_dB - 0.002*KF_dB.^2 + 0.0001*KF_dB.^3;
    C_phi_LOS( KF_dB < 0 & C_phi_LOS < 1 ) = 1;
    C_phi = C_phi_NLOS * C_phi_LOS;
    
    if alt_AS_mapping % Alternative method leading to correct calibration results
        AoA = 2*rand(size(pow))-1;                          % Random angles in range -57 ... 57 deg
        AoA = AoA .* ceil( -pow+(max( pow,[],2)*oL) );      % Set angle of strongest path to 0 deg
        as  = qf.calc_angular_spreads( AoA , pow );         % Calculate AS
        AoA = AoA .* (ASA./as * oL);
        
        AoA_LOS = angle( exp( 1j*angles(:,2) ) );           % Apply LOS angle
        AoA = AoA - ( AoA(:,1) - AoA_LOS )*oL;
        
    else % Strict method as defined in the TR
        % The AOAs are determined by applying the inverse Gaussian function with input parameters pow
        % and RMS angle spread ASA, eq. (7.5-9)
        AoA = 2*(ASA/1.4)*oL .* sqrt( -log( pow ./ (max( pow,[],2)*oL) ) ) ./ (C_phi*oL);
        
        % Assign positive or negative sign and introduce random variation, enforce the first cluster to
        % the LOS direction, eq. (7.5-12)
        X = 2*(randi(2,N,L)-1.5);
        Y = randn( N,L ).*((ASA/7)*oL);
        AoA_LOS = angle( exp( 1j*angles(:,2) ) );
        AoA = ( X.*AoA + Y ) - ( X(:,1).*AoA(:,1) + Y(:,1) - AoA_LOS )*oL;
    end
    
    % Wrap AoAs around the unit circle
    AoA = angle( exp( 1j*AoA ) );
    
    if alt_AS_mapping % Alternative method leading to correct calibration results
        AoD = 2*rand(size(pow))-1;                          % Random angles in range -57 ... 57 deg
        AoD = AoD .* ceil( -pow+(max( pow,[],2)*oL) );      % Set angle of strongest path to 0 deg
        as  = qf.calc_angular_spreads( AoD , pow );         % Calculate AS
        AoD = AoD .* (ASD./as * oL);
        
        AoD_LOS = angle( exp( 1j*angles(:,1) ) );           % Apply LOS angle
        AoD = AoD - ( AoD(:,1) - AoD_LOS )*oL;
        
    else % Strict method as defined in the TR
        % The generation of AOD follows a procedure similar to AOA
        AoD = 2*(ASD/1.4)*oL .* sqrt( -log( pow ./ (max( pow,[],2)*oL) ) ) ./ (C_phi*oL);
        X = 2*(randi(2,N,L)-1.5);
        Y = randn( N,L ).*((ASD/7)*oL);
        AoD_LOS = angle( exp( 1j*angles(:,1) ) );
        AoD = ( X.*AoD + Y ) - ( X(:,1).*AoD(:,1) + Y(:,1) - AoD_LOS )*oL;
    end
    
    % Wrap AoDs around the unit circle
    AoD = angle( exp( 1j*AoD ) );
    
    % C_phi_NLOS is defined as a scaling factor related to the total number of clusters.
    % The "-1" accounts for the LOS component which is always present in QuaDRiGa.
    if L <= 9
        C_phi_NLOS = 0.889;
    elseif L >= 21
        C_phi_NLOS = 1.178;
    else
        L_phi_NLOS = [     8,    10,    11,    12,     15,    19,    20 ];
        C_phi_NLOS = [ 0.889, 0.957, 1.031, 1.104, 1.1088, 1.184, 1.178 ];
        C_phi_NLOS = qf.interp( L_phi_NLOS, 0, C_phi_NLOS, L-1 );
    end
    
    % C_phi_LOS depends on the K-Factor, eq. (7.5-15)
    C_phi_LOS = 1.3086 + 0.0339*KF_dB - 0.0077*KF_dB.^2 + 0.0002*KF_dB.^3;
    C_phi_LOS( KF_dB < 0 & C_phi_LOS < 1 ) = 1;
    C_phi = C_phi_NLOS * C_phi_LOS;
    
    if alt_AS_mapping % Alternative method leading to correct calibration results
        EoA = 2*rand(size(pow))-1;                          % Random angles in range -57 ... 57 deg
        EoA = EoA .* ceil( -pow+(max( pow,[],2)*oL) );      % Set angle of strongest path to 0 deg
        as  = qf.calc_angular_spreads( EoA , pow );         % Calculate AS
        EoA = EoA .* (ESA./as * oL);
        
        EoA_LOS = angle( exp( 1j*angles(:,4) ) );
        EoA = EoA - ( EoA(:,1) - EoA_LOS )*oL;
        
    else % Strict method as defined in the TR
        % The generation of ZOA assumes that the composite PAS in the zenith dimension of all clusters
        % is Laplacian. The ZOAs are determined by applying the inverse Laplacian function, eq. (7.5-14)
        EoA = -(ESA*oL) .* log( pow ./ (max( pow,[],2)*oL) ) ./ (C_phi*oL);
        
        % Assign positive or negative sign and introduce random variation, enforce the first cluster to
        % the LOS direction, eq. (7.5-12)
        X = 2*(randi(2,N,L)-1.5);
        Y = randn( N,L ).*((ESA/7)*oL);
        EoA_LOS = angle( exp( 1j*angles(:,4) ) );
        EoA = ( X.*EoA + Y ) - ( X(:,1).*EoA(:,1) + Y(:,1) - EoA_LOS )*oL;
    end
    
    % If ZOA is is wrapped within [0, 360°], and ZOA is in [ 180° , 360° ], then set to ( 360° - ZOA )
    % QuaDRiGa uses EOA. The wrapping is equivalent to inverting the real part.
    C = exp( 1j*EoA );
    ii = real(C) < 0;
    C(ii) = -real(C(ii)) + 1j*imag(C(ii));
    EoA = angle( C );
    
    % The generation of EOD follows a procedure similar to EOA
    if alt_AS_mapping % Alternative method leading to correct calibration results
        EoD = 2*rand(size(pow))-1;                          % Random angles in range -57 ... 57 deg
        EoD = EoD .* ceil( -pow+(max( pow,[],2)*oL) );      % Set angle of strongest path to 0 deg
        as  = qf.calc_angular_spreads( EoD , pow );         % Calculate AS
        EoD = EoD .* (ESD./as * oL);
        
        EoD_LOS = angle( exp( 1j*angles(:,3) ) );
        EoD = EoD - ( EoD(:,1) - EoD_LOS )*oL;
        
    else % Strict method as defined in the TR
        EoD = -(ESD*oL) .* log( pow ./ (max( pow,[],2)*oL) ) ./ (C_phi*oL);
        X = 2*(randi(2,N,L)-1.5);
        Y = randn( N,L ).*((ESD/7)*oL);
        EoD_LOS = angle( exp( 1j*angles(:,3) ) );
        EoD = ( X.*EoD + Y ) - ( X(:,1).*EoD(:,1) + Y(:,1) - EoD_LOS )*oL;
    end

    C = exp( 1j*EoD );
    ii = real(C) < 0;
    C(ii) = -real(C(ii)) + 1j*imag(C(ii));
    EoD = angle( C );
    
%     % Remove clusters with less than -25 dB power compared to the maximum cluster power.
%     ind = pow ./ (max( pow,[],2 ) * oL) < 10.^(0.1*-25);
%     pow(ind) = 0;
    
end


end
