function [ fbs_pos, lbs_pos, AoD_c, AoA_c, EoD_c, EoA_c ] = generate_fbs_lbs( tx_pos, rx_pos, taus, ...
    AoD, AoA, EoD, EoA, NumSubPaths, SubPathCPL, PerClusterAS, method )
%GENERATE_FBS_LBS Generates the positions of the FBS and the LBS
%
% Input:
%   tx_pos          Tx-positions (metric)                               [ 3 x N ]
%   rx_pos          Rx-positions (metric)                               [ 3 x N ]
%   taus            Cluster delays (seconds)                            [ N x L ]
%   AoD             Azimuth angles of departure (rad)                   [ N x L ]
%   AoA             Azimuth angles of arrival (rad)                     [ N x L ]
%   EoD             Elevation angles of departure (rad)                 [ N x L ]
%   EoA             Elevation angles of arrival (rad)                   [ N x L ]
%   NumSubPaths     The number of sub-paths per cluster                 [ 1 x L ] - optional
%   SubPathCPL      Uniform random numbers for the subpath copling      [ 4 x ML x F ] - optional
%   PerClusterAS    Per-cluster angular spreads in (deg)                [ 4 x F ] - optional
%   method          String selecting the scatterer generation method
%
% If NumSubPaths is not given, then there is only one sub-path per cluster. The variables
% SubPathCPL and PerClusterAS are ignored in this case.
%
% Derived from input:
%   N   = number of users (from rx positions)
%   L   = number of paths (from taus)
%   Ln  = NLOS clusters (from angles)
%   Lf  = Deterministic clusters (LOS + ground reflection)
%   ML  = total number of sub-paths = sum(NumSubPaths)
%   MLn = total number of NLOS sub-paths
%   F   = number of frequencies (from SubPathCPL)
%
% Supported scatterer generation methods specified by "method":
%   bounce2         Splits the path in two bounces and mimizes the distance between FBS and LBS. If
%                   this is not possible, a single-bounce model is used where FBS = LBS. The
%                   departure angle change in this case.
%
% Content of variable PerClusterAS:
%   PerClusterAS_D (deg)
%   PerClusterAS_A (deg)
%   PerClusterES_D (deg)
%   PerClusterES_A (deg)
%
% Output:
%   fbs_pos     	Positions of the first-bounce-scatterers            [ 3 x ML x N x F ]
%   lbs_pos     	Positions of the last-bounce-scatterers             [ 3 x ML x N x F ]
%   AoD_c           Corrected azimuth angles of departure (rad)         [ N x L ]
%   AoA_c           Corrected azimuth angles of arrival (rad)           [ N x L ]
%   EoD_c           Corrected elevation angles of departure (rad)       [ N x L ]
%   EoA_c           Corrected elevation angles of arrival (rad)         [ N x L ]
%
% The function relies on solving an optimization problem that might have no solution. In this case,
% thee might be changes in the angles that make the problem solvable. The corrected angles contain
% these chenges.
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

if ~exist('rx_pos','var') || size( rx_pos,1) ~= 3
    error('QuaDRiGa:qd_builder:generate_fbs_lbs','Rx position is not given or has invalid format.')
end

% Number of users
N = size( rx_pos, 2 );
oN = ones(1,N);

if ~exist('tx_pos','var') || size( tx_pos,1) ~= 3 || size( tx_pos,2) ~= N
    error('QuaDRiGa:qd_builder:generate_fbs_lbs','Tx position is not given or has invalid format.')
end

if ~exist('taus','var') || size( taus,1) ~= N || ~isreal(taus)
    error('QuaDRiGa:qd_builder:generate_fbs_lbs','"taus" is not given or has invalid format.')
end

% Number of clusters
L = size( taus, 2 );
oL = ones(1,L);

if ~exist('AoD','var') || size( AoD,1) ~= N || size( AoD,2) ~= L || ~isreal(AoD)
    error('QuaDRiGa:qd_builder:generate_fbs_lbs','"AoD" is not given or has invalid format.')
end
if ~exist('AoA','var') || size( AoA,1) ~= N || size( AoA,2) ~= L || ~isreal(AoA)
    error('QuaDRiGa:qd_builder:generate_fbs_lbs','"AoA" is not given or has invalid format.')
end
if ~exist('EoD','var') || size( EoD,1) ~= N || size( EoD,2) ~= L || ~isreal(EoD)
    error('QuaDRiGa:qd_builder:generate_fbs_lbs','"EoD" is not given or has invalid format.')
end
if ~exist('EoA','var') || size( EoA,1) ~= N || size( EoA,2) ~= L || ~isreal(EoA)
    error('QuaDRiGa:qd_builder:generate_fbs_lbs','"EoA" is not given or has invalid format.')
end

if ~exist('NumSubPaths','var') || isempty( NumSubPaths ) 
    NumSubPaths = oL;
elseif size( NumSubPaths,1) ~= 1 || size( NumSubPaths,2) ~= L
    error('QuaDRiGa:qd_builder:generate_fbs_lbs','"NumSubPaths" is not given or has invalid format.')
end

% Total number of sub-paths
ML = sum( NumSubPaths );
oML = ones(1,ML);

if ~exist('SubPathCPL','var') || isempty( SubPathCPL ) 
    SubPathCPL = ones(4,ML);
elseif size( SubPathCPL,1) ~= 4 || size( SubPathCPL,2) ~= ML
    error('QuaDRiGa:qd_builder:generate_fbs_lbs','"SubPathCPL" is not given or has invalid format.')
end

% Number of Frequencies
F = size( SubPathCPL,3 );
oF = ones(1,F);

if ~exist('PerClusterAS','var') || isempty( PerClusterAS ) 
    PerClusterAS = zeros(4,F);
elseif size( PerClusterAS,1) ~= 4 || ~( size( PerClusterAS,2) == 1 || size( PerClusterAS,2) == F )
    error('QuaDRiGa:qd_builder:generate_fbs_lbs','"PerClusterAS" is not given or has invalid format.')
end
if size( PerClusterAS,2) == 1
    PerClusterAS = PerClusterAS(:,oF);
end

% Distances between BS and MT
d_2d = hypot( tx_pos(1,:) - rx_pos(1,:), tx_pos(2,:) - rx_pos(2,:) );
d_2d( d_2d<1e-5 ) = 1e-5;

% Calculate angles between BS and MT
angles = zeros( 5,N );
angles(1,:) = atan2( rx_pos(2,:) - tx_pos(2,:) , rx_pos(1,:) - tx_pos(1,:) );   % Azimuth at BS
angles(2,:) = pi + angles(1,:);                                                 % Azimuth at MT
angles(3,:) = atan2( ( rx_pos(3,:) - tx_pos(3,:) ), d_2d );                     % Elevation at BS
angles(4,:) = -angles(3,:);                                                     % Elevation at MT
angles = angles.';

% Deterimne the number if direct paths
if NumSubPaths(1) == 1 && ...
        all( abs( exp(1j*AoD(:,1)) - exp(1j*angles(:,1)) ) < 1e-12 ) && ...
        all( abs( EoD(:,1) - angles(:,3) ) < 1e-12 ) && ...
        all( abs( exp(1j*AoA(:,1)) - exp(1j*angles(:,2)) ) < 1e-12 ) && ...
        all( abs( EoA(:,1) - angles(:,4) ) < 1e-12 ) && ...
        all( abs( taus(:,1) ) < 1e-12 )
    if L>1 && NumSubPaths(2) == 1 && ... 
            all( abs( AoD(:,2) - AoD(:,1) ) < 1e-12 ) && ...    % GR has same azimuth angles as LOS path has.
            all( abs( AoA(:,2) - AoA(:,1) ) < 1e-12 )
        Lf = 2;     % First path is LOS path and second path is GR path
    else
        Lf = 1;     % First path is LOS path
    end
else
    Lf = 0;         % Only NLOS paths are given
end

if Lf == 0 && any( abs(taus( :,1 )) < 1e-12 )
    error('QuaDRiGa:qd_builder:generate_fbs_lbs','No LOS paths are defined, but delays have a value of 0.')
end

% Generate NLOS sub-paths and apply sub-path coupling

% The optimized order of the 20 sub-paths 
offset = [-0.8844 1.1481 0.6797 -1.1481 -0.6797 1.5195 -0.0447 0.3715 -1.5195,...
    -0.2492 -0.5129 0.5129 0.8844 -2.1551 -0.1413 0.1413 0.0447 2.1551 -0.3715 0.2492 ];

% Combine cluster angles [ N x L x 4 ]
Ac = cat(3,AoD,AoA,EoD,EoA);

% Placeholder for the sub-path angles  [ N x ML x 4 x F ]
As = zeros( N,ML,4,F );

% The start and end-indices for the sub-paths [ 1 x L ]
Ms = [0,cumsum(NumSubPaths(1:end-1))]+1;
Me = cumsum(NumSubPaths);

for iF = 1 : F              % Frequency
    for iA = 1 : 4          % AoD, AoA, EoD, EoA

        % Apply LOS angles
        As(:,1:Lf,iA,iF) = Ac(:,1:Lf,iA);
        
        % Create sub-path angles
        for iL = Lf+1 : L       % Cluster

            % Get the offset angles
            if NumSubPaths(iL) == 1
                of = 0;
            else
                of = offset( 1:NumSubPaths(iL) );
                if NumSubPaths(iL) < 20
                    of = (of-mean(of));
                    of = of./sqrt(mean(of.^2));
                end
                of = of .* 0.017453292519943;  % deg to rad
            end
            
            % Add subpath-angles to cluster angles
            ax = oN' * PerClusterAS(iA,iF) * of;% + Ac(:,iL,iA) * ones(1,NumSubPaths(iL));
            
            % Change order to account for random sub-path coupling
            [~,ind] = sort( SubPathCPL(iA,Ms(iL):Me(iL),iF) );
            
            % Store subpath angles
            As(:,Ms(iL):Me(iL),iA,iF) = ax(:,ind);
        end
    end
    
    % Apply cluster angles
    for iL = Lf+1 : L
        [As(:,Ms(iL):Me(iL),1,iF), As(:,Ms(iL):Me(iL),3,iF)] =...
            rotate_sphere( As(:,Ms(iL):Me(iL),1,iF), As(:,Ms(iL):Me(iL),3,iF), Ac(:,iL,1), Ac(:,iL,3) );
        [As(:,Ms(iL):Me(iL),2,iF), As(:,Ms(iL):Me(iL),4,iF)] =...
            rotate_sphere( As(:,Ms(iL):Me(iL),2,iF), As(:,Ms(iL):Me(iL),4,iF), Ac(:,iL,2), Ac(:,iL,4) );
    end
end

% Create clusters

% The distance vector from the Tx to the initial position
r = rx_pos - tx_pos;
norm_r = sqrt(sum(r.^2)).';

% Get the total path length of the NLOS component
dist = clst_expand( taus, NumSubPaths );  % Delay for each sub-path
dist = dist * speed_of_light + norm_r(:,oML);

fbs_pos = zeros(N,ML,3,F);
lbs_pos = zeros(N,ML,3,F);

for iF = 1 : F
    
    phi_d_lm    = As(:,:,1,iF); % AoD [ N x ML ]
    theta_d_lm  = As(:,:,3,iF); % EoD [ N x ML ]
    phi_a_lm    = As(:,:,2,iF); % AoA [ N x ML ]
    theta_a_lm  = As(:,:,4,iF); % EoA [ N x ML ]
    
    % Get the direction of the last bounce scatterer (LBS) seen from the receivers position.
    [ ahat_lm_x, ahat_lm_y, ahat_lm_z ] = sph2cart(phi_a_lm, theta_a_lm, 1);
    ahat_lm = cat(3, ahat_lm_x, ahat_lm_y, ahat_lm_z );
    ahat_lm = permute(ahat_lm,[3,1,2]);  % [ 3 x N x ML  ]
    
    % Get the direction of the first bounce scatterer (FBS) seen from the transmitter center position.
    [ bhat_lm_x, bhat_lm_y, bhat_lm_z ] = sph2cart(phi_d_lm, theta_d_lm, 1);
    bhat_lm = cat(3, bhat_lm_x, bhat_lm_y, bhat_lm_z );
    bhat_lm = permute(bhat_lm,[3,1,2]);  % [ 3 x N x ML ]
    
    % This implements the "bounce2" model
    [ norm_a_lm, norm_b_lm,~,valid ] = solve_multi_bounce_opti( ahat_lm, bhat_lm, r, dist, 2 );
    norm_a_lm = permute( norm_a_lm, [2,3,1] );
    norm_b_lm = permute( norm_b_lm, [2,3,1] );
    valid = permute( valid, [2,3,1] );

    if Lf > 0.5     % LOS
        norm_a_lm(:,1) = 0.5*norm_r;
        norm_b_lm(:,1) = 0.5*norm_r;
        valid(:,1) = true;
    end
    if Lf > 1.5     % GR
       valid(:,2) = false;  % Use single-bounce for ground reflection
    end
        
    % For the invalid paths, use single-bounce-model
    iv = ~valid;
    if any( iv(:) )
        ahat_iv = [ reshape( ahat_lm_x(iv),1,[]) ; reshape( ahat_lm_y(iv),1,[]) ; reshape( ahat_lm_z(iv),1,[]) ];
        rx = r(oML,:)'; ry = r(oML*2,:)'; rz = r(oML*3,:)';
        rx = [ reshape(rx(iv),1,[]) ; reshape(ry(iv),1,[]) ; reshape(rz(iv),1,[]) ];
        
        [ b, norm_b_lm(iv), norm_a_lm(iv) ] = solve_cos_theorem( ahat_iv , rx , reshape(dist(iv),1,[])  );
        
        tmp = reshape(norm_b_lm(iv),1,[]);
        bhat_lm_x(iv) = b(1,:)./tmp;
        bhat_lm_y(iv) = b(2,:)./tmp;
        bhat_lm_z(iv) = b(3,:)./tmp;
    end
       
    % Calculate the FBS position (relative to initial Rx-pos)
    fbs_pos(:,:,1,iF) = norm_b_lm .* bhat_lm_x - r(oML,:).';
    fbs_pos(:,:,2,iF) = norm_b_lm .* bhat_lm_y - r(oML*2,:).';
    fbs_pos(:,:,3,iF) = norm_b_lm .* bhat_lm_z - r(oML*3,:).';

    % Calculate the LBS position
    lbs_pos(:,:,1,iF) = norm_a_lm .* ahat_lm_x;
    lbs_pos(:,:,2,iF) = norm_a_lm .* ahat_lm_y;
    lbs_pos(:,:,3,iF) = norm_a_lm .* ahat_lm_z;
end

% Calculate the scatterer positions in global coordinates.
fbs_pos = fbs_pos + permute( rx_pos(:,:,oML,oF) , [2,3,1,4] );
lbs_pos = lbs_pos + permute( rx_pos(:,:,oML,oF) , [2,3,1,4] );

% Reshape output to [ 3 x ML x N x F ]
fbs_pos = permute( fbs_pos,[3,2,1,4] );
lbs_pos = permute( lbs_pos,[3,2,1,4] );

% Update the angles from the scatterer positions
AoD_c = zeros(N,L);
AoA_c = zeros(N,L);
EoD_c = zeros(N,L);
EoA_c = zeros(N,L);

for iL = 1 : L
    % Calculate average departure angles
    A = fbs_pos(:,Ms(iL):Me(iL),:,:);
    A = reshape( permute( A,[3,1,2,4] ),N,3,[] );
    if size( A,3 ) > 1
        A = mean( A,3 );
    end
    A = A-tx_pos';
    AoD_c(:,iL) = atan2(A(:,2),A(:,1) );
    EoD_c(:,iL) = asin( A(:,3)./sqrt(sum(A.^2,2)) );
    
    % Calculate average arrival angles
    A = lbs_pos(:,Ms(iL):Me(iL),:,:);
    A = reshape( permute( A,[3,1,2,4] ),N,3,[] );
    if size( A,3 ) > 1
        A = mean( A,3 );
    end
    A = A-rx_pos';
    AoA_c(:,iL) = atan2(A(:,2),A(:,1) );
    EoA_c(:,iL) = asin( A(:,3)./sqrt(sum(A.^2,2)) );
    
end

end
