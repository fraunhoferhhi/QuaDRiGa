function [ pow, taus, AoD, AoA, EoD, EoA ] = generate_paths( L, tx_pos, rx_pos, spreads, gr_epsilon_r, path_sos, r_DS )
%GENERATE_PATHS Generate multipath components
%
% Input:
%   L               Number of paths (including LOS and GR) [ 1 x 1 ]
%   tx_pos          Fixed Tx (Single mobility): Tx-position (metric) [ 3 x 1 ]
%                   Dual mobility: Tx-position (metric) [ 3 x N ]
%   rx_pos          Rx-positions (metric) [ 3 x N ]
%   spreads         Vector of delay and angular spreads (see below) [ 6 x N x F ]
%   gr_epsilon_r    Relative permeability for ground reflection [ N x F ]
%                   if set to empty [], ground reflection is disabled (default)
%   path_sos        Pre-initialized SOS generators for the parameter generation [ Ln x 5 ]
%                   Column 1 : Delays (Uniform)
%                   Column 2 : AoA (Normal)
%                   Column 3 : AoD (Normal)
%                   Column 4 : EoA (Normal)
%                   Column 5 : EoD (Normal)
%                   Alternative : Correlation distance for spatial consistency in [m]; [ 1 x 1 ]
%   r_DS            Delay distribution proportionality factor (must be > 1); [ 1 x 1 ]
%                   Is only used in single-Frequency simulations
%
% Derived from input:
%   L  = number of paths
%   Ln = number of NLOS paths ( L-2 for ground-reflection, L-1 otherwise )
%   N  = number of users (from rx positions)
%   F  = number of frequencies (from spreads)
%
% Content of variable spreads:
%   DS  (s)
%   ASD (deg)
%   ASA (deg)
%   ESD (deg)
%   ESA (deg)
%   KF  (linear)
%
% Output:
%   pow             Path powers (linear, normalized to sum 1) [ N x L x F ]
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
angles(3,:) = atan2( ( rx_pos(3,:) - tx_pos(3,:) ), d_2d );                    % Elevation at BS
angles(4,:) = -angles(3,:);                                                     % Elevation at MT
angles(5,:) = -atan2( ( rx_pos(3,:) + tx_pos(3,:) ), d_2d );                   % Ground Reflection Elevation at BS and MT
angles = angles.';

if ~exist('spreads','var') || size( spreads,1) ~= 6 || size( spreads,2) ~= N
    error('QuaDRiGa:qd_builder:generate_paths','Spreads are not given or have invalid format.')
end

% Number of frequencies
F = size( spreads, 3 );
oF = ones(1,F);

% Split spread variables into individual variables and convert from [deg] to [rad]
spreadsP = permute( spreads , [2,3,1] );
DS  = spreadsP(:,:,1);
ASD = spreadsP(:,:,2) * pi/180;
ASA = spreadsP(:,:,3) * pi/180;
ESD = spreadsP(:,:,4) * pi/180;
ESA = spreadsP(:,:,5) * pi/180;
KF  = spreadsP(:,:,6);

if ~exist('gr_epsilon_r','var') || isempty( gr_epsilon_r )
    gr_enabled = false;
elseif size( gr_epsilon_r,1) ~= F || size( gr_epsilon_r,2) ~= N
    error('QuaDRiGa:qd_builder:generate_paths','Ground reflection coefficient has invalid format.')
else
    gr_enabled = true;
    gr_epsilon_r = gr_epsilon_r.';
end

% Process ground reflection parameters
if gr_enabled
    % GR path length
    d_gf = sqrt( sum(([rx_pos(1:2,:);-rx_pos(3,:)] - tx_pos).^2,1) );
    
    % Incident angle of the GR
    theta_r = -angles(:,5);
    
    % Delay of the ground reflection relative to the LOS component
    tau_gr = ( d_gf-d_3d ).' / speed_of_light;
    
    % The reflection coefficient
    Z         = sqrt( gr_epsilon_r - (cos(theta_r)).^2 * oF );
    R_par     = (gr_epsilon_r .* (sin(theta_r)*oF) - Z) ./ (gr_epsilon_r .* (sin(theta_r)*oF) + Z);
    R_per     = ( sin(theta_r)*oF - Z) ./ ( sin(theta_r)*oF + Z);
    Rsq       = ( 0.5*(abs(R_par).^2 + abs(R_per).^2) );
end

% Get the number of deterministic and non-deterministic paths
if gr_enabled
    Lf = 2;
    Ln = L-2;
else
    Lf = 1;
    Ln = L-1;
end
oLn = ones(1,Ln);
oLf = ones(1,Lf);
oL  = ones(1,L);

if ~exist('path_sos','var') || isempty( path_sos )
    path_sos = [];
    
elseif isnumeric( path_sos ) && numel( path_sos ) == 1
    
    if path_sos == 0
        path_sos = [];
    else % Create sos_model
        SC_lambda = path_sos;
        acf = 'Comb300';
        path_sos = qd_sos([]);
        for i_cluster = 1 : Ln
            % Delays (are identical if Tx and Rx positions are swapped)
            path_sos(i_cluster,1) = qd_sos( acf, 'Uniform', SC_lambda ); % Delays
            path_sos(i_cluster,1).sos_phase(:,2) = path_sos(i_cluster,1).sos_phase(:,1);
            
            % Azimuth angles (depature become arrival angles if positions are swapped)
            path_sos(i_cluster,2) = qd_sos( acf, 'Uniform', SC_lambda ); % ASD
            path_sos(i_cluster,3) = qd_sos( acf, 'Uniform', SC_lambda ); % ASA
            path_sos(i_cluster,2).sos_phase(:,1) = path_sos(i_cluster,3).sos_phase(:,2);
            path_sos(i_cluster,3).sos_phase(:,1) = path_sos(i_cluster,2).sos_phase(:,2);
            
            % Elevation angles (depature become arrival angles if positions are swapped)
            path_sos(i_cluster,4) = qd_sos( acf, 'Uniform', SC_lambda ); % ASD
            path_sos(i_cluster,5) = qd_sos( acf, 'Uniform', SC_lambda ); % ASA
            path_sos(i_cluster,4).sos_phase(:,1) = path_sos(i_cluster,5).sos_phase(:,2);
            path_sos(i_cluster,5).sos_phase(:,1) = path_sos(i_cluster,4).sos_phase(:,2);
        end
    end
    
elseif isa( path_sos, 'qd_sos' ) && size(path_sos,1) == Ln && size(path_sos,2) == 5
    % Test if all generators are set to uniform distribution
    for n = 1 : size( path_sos, 1)
        for m = 1 : size( path_sos, 2)
            if ~strcmp( path_sos(n,m).distribution , 'Uniform' )
                warning('QuaDRiGa:qd_builder:generate_paths:Uniform',...
                    'All random generators should be set to uniform distribution.');
                warning ('off','QuaDRiGa:qd_builder:generate_paths:Uniform');
            end
        end
    end
    warning ('on','QuaDRiGa:qd_builder:generate_paths:Uniform');
    
    % Test if delays are reciprocal
    tmp = cat(3,path_sos(:,1).sos_phase);
    tmp = tmp(:,1,:) - tmp(:,2,:);
    if any(abs(tmp(:))>1e-6)
        warning('QuaDRiGa:qd_builder:generate_paths','Cluster delays are not reciprocal.');
    end
    
    % Test if azimuth angles are reciprocal
    tmp1 = cat(3,path_sos(:,2).sos_phase);
    tmp2 = cat(3,path_sos(:,3).sos_phase);
    tmp3 = tmp1(:,1,:) - tmp2(:,2,:);
    tmp4 = tmp2(:,1,:) - tmp1(:,2,:);
    if any( abs( [tmp3(:);tmp4(:)] ) >1e-6 )
        warning('QuaDRiGa:qd_builder:generate_paths','Azimuth angles are not reciprocal.');
    end
    
    % Test if elevation angles are reciprocal
    tmp1 = cat(3,path_sos(:,4).sos_phase);
    tmp2 = cat(3,path_sos(:,5).sos_phase);
    tmp3 = tmp1(:,1,:) - tmp2(:,2,:);
    tmp4 = tmp2(:,1,:) - tmp1(:,2,:);
    if any( abs( [tmp3(:);tmp4(:)] ) >1e-6 )
        warning('QuaDRiGa:qd_builder:generate_paths','Elevation angles are not reciprocal.');
    end
else
    error('QuaDRiGa:qd_builder:generate_paths','path_sos has invalid format.');
end

if ~exist('r_DS','var') || isempty( r_DS )
    r_DS = 2.2;
end

% Generate paths
if L == 1   % Only LOS component is present
    
    pow  = ones(N,1,F);
    taus = zeros(N,1);
    AoD  = angles(:,1);
    AoA  = angles(:,2);
    EoD  = angles(:,3);
    EoA  = angles(:,4);
    
elseif gr_enabled && L == 2   % Only LOS and GR componenet
    
    pow  = zeros( N, L, F );
    pow( :,1,: ) = (1-Rsq);
    pow( :,2,: ) = Rsq;
    taus = [ zeros(N,1), tau_gr ];
    AoD  = angles(:,1) * [1 1];
    AoA  = angles(:,2) * [1 1];
    EoD  = [ angles(:,3), angles(:,5) ];
    EoA  = [ angles(:,4), angles(:,5) ];
    
else % Additional NLOS componenets are present
    
    % Placeholder for the output variables
    pow  = zeros( N, L, F );
    taus = zeros( N, L );
    
    % Generate initial delays following an exponenetial distribution
    if isempty( path_sos )
        randC = rand( N,Ln );
    else
        randC = val( path_sos( 1:Ln,1 ), rx_pos,tx_pos_SOS ).';  % Uniform distribution
    end
    taus(:,Lf+1:end) = -log( randC );
    
    % Add ground reflection delay and calculate split the LOS power in a LOS and GR part
    if gr_enabled
        taus(:,2) = tau_gr;
        Pf = [ (1-mean(Rsq,2)) mean(Rsq,2) ];
    else % only LOS
        Pf = oN';
    end
    
    % Generate correlated random variables for the angles
    if Ln > 2
        if isempty( path_sos )
            randI = rand( N,Ln,4 );
        else
            % Generate Normal-distribute spatially correlated random variables
            randI = zeros( N,Ln,4 );
            for i_ang = 1:4                                             % Uniform distribution
                randI(:,:,i_ang) = val( path_sos( 1:Ln,i_ang+1 ), rx_pos,tx_pos_SOS ).';
            end
        end
    end
    
    % Generate angles for the average angular spreads
    mu_fixed = zeros(N,4);  % The mean fixed angle for LOS and GR
    path_angles = zeros( N,L,4 );
    for i_ang = 1 : 4
        
        % The LOS angle and the angle for the GR are deterministic and know from the BS and MT
        % positions. They are assembled before generating the NLOS angles. NLOS angles are
        % distributed around the average fixed angle.
        
        if gr_enabled
            if i_ang == 3 % EoD
                ang_fixed = [ angles(:,3) angles(:,5) ];
            elseif i_ang == 4 % EoA
                ang_fixed = [ angles(:,4) angles(:,5) ];
            else % AoD and AoA
                ang_fixed = angles(:,i_ang) * [1 1];
            end
        else % LOS only
            ang_fixed = angles( :,i_ang );
        end
        [ ~, mu_fixed(:,i_ang) ] = qf.calc_angular_spreads( ang_fixed, Pf, 0 );
        
        % Scale angles within +/- pi
        if Ln == 1
            randC = ones(N,1) / sqrt(2);    % 40.5 degree
        elseif Ln == 2
            randC = ones(N,1) * [-1 1] / sqrt(2);       % +/- 40.5 degree
        else % Ln > 2
            % Read variables from previously generated random numbers
            randC = (randI(:,:,i_ang)-0.5)*pi;
        end
        
        % Set the fixed angles (LOS and GR)
        path_angles(:,:,i_ang) = [ ang_fixed - mu_fixed(:,oLf*i_ang), randC ];
    end
    
    % Generate initial path powers
    % Calculate DS scaling coefficient
    if F > 1
        gDS = DS ./ (( max(DS,[],2) + min(DS,[],2) ) *oF);
        gDS( gDS < 0.15 ) = 0.15;
        gDS( gDS > 0.85 ) = 0.85;
        gDS = -1.5 * log( 1.2 * gDS - 0.15 );
        
        % If we have identical DS values for multiple frequencies, we use the scaling factor
        ind = abs( max(DS,[],2) ./ min(DS,[],2) - 1 ) < 1e-6;
        gDS(ind,:) = r_DS-1;
    else
        gDS = r_DS(oN,1) - 1;
    end
        
    % Calculate ASD scaling coefficient
    if F > 1
        gASD = 0.75 * ASD./(max(ASD,[],2)*oF); 
        gASD( gASD < 0.25 ) = 0.25;
        gASD = -2.2 * log( 1.5 * gASD - 0.35 );
    else
        gASD = 0.56*oN';
    end
    
    % Calculate ASA scaling coefficient
    if F > 1
        gASA = 0.75 * ASA./(max(ASA,[],2)*oF); 
        gASA( gASA < 0.25 ) = 0.25;
        gASA = -2.2 * log( 1.5 * gASA - 0.35 );
    else
        gASA = 0.56*oN';
    end
    
    % Calculate ESD scaling coefficient
    if F > 1
        gESD = 0.75 * ESD./(max(ESD,[],2)*oF); 
        gESD( gESD < 0.25 ) = 0.25;
        gESD = -3.4 * log( 1.2 * gESD - 0.1 );
    else
        gESD = 0.76*oN';
    end
    
    % Calculate ESA scaling coefficient
    if F > 1
        gESA = 0.75 * ESA./(max(ESA,[],2)*oF); 
        gESA( gESA < 0.25 ) = 0.25;
        gESA = -3.4 * log( 1.2 * gESA - 0.1 );
    else
        gESA = 0.76*oN';
    end
    
    for i_freq = 1 : F
        % NLOS powers
        pow(:,:,i_freq) = exp( -taus .* gDS(:,oL*i_freq) ) .* ...
            exp( -path_angles(:,:,1).^2 .* gASD(:,oL*i_freq) ) .* ...
            exp( -path_angles(:,:,2).^2 .* gASA(:,oL*i_freq) ) .* ...
            exp( -abs(path_angles(:,:,3)) .* gESD(:,oL*i_freq) ) .* ...
            exp( -abs(path_angles(:,:,4)) .* gESA(:,oL*i_freq) );
        
        % Split LOS power into LOS and GR power
        if gr_enabled
            Pf = [ 1-Rsq(:,i_freq), Rsq(:,i_freq) ];
        else
            Pf = oN';
        end
        
        % Apply K-Factor
        Pf = Pf .* ((KF(:,i_freq) ./ ( 1 + KF(:,i_freq) ))*oLf);
        Pn = 1-sum(Pf,2);
        
        pow(:,1:Lf,i_freq) = 0;
        pow(:,:,i_freq)    = pow(:,:,i_freq) .* (Pn ./ sum(pow(:,:,i_freq),2) * oL);
        pow(:,1:Lf,i_freq) = Pf;
    end
    
    % Applying the Delay Spread
    % Calculate the delay scaling factor (has analytic solution)
    s = zeros( N,F );   % Scaling factor
    for i_freq = 1 : F
        if gr_enabled
            a = DS(:,i_freq).^2;
            b = pow(:,2,i_freq) .* tau_gr.^2;
            c = sum( pow(:,3:end,i_freq) .* taus(:,3:end).^2 , 2 );
            d = pow(:,2,i_freq) .* tau_gr;
            e = sum( pow(:,3:end,i_freq) .* taus(:,3:end) , 2 );
            tmp = a.*c - a.*e.^2 - b.*c + b.*e.^2 + c.*d.^2;
            tmp( tmp<0 ) = 0;   % Minimum DS is given by GR, cannot get smaller
            s(:,i_freq) = ( sqrt( tmp ) + d.*e ) ./ (c-e.^2);
        else
            a = DS(:,i_freq).^2;
            b = sum( pow(:,2:end,i_freq) .* taus(:,2:end) , 2 );
            c = sum( pow(:,2:end,i_freq) .* taus(:,2:end).^2 , 2 );
            s(:,i_freq) = sqrt( a ./ (c-b.^2) );
        end
    end
    if F > 1 % Average scaling coefficient for multiple frequencies
        s = mean(s,2);
    end
    taus(:,Lf+1:end) = taus(:,Lf+1:end) .* ( s*oLn );  % Scale the delays
    
    
    % Applying the Angular Spread
    % Combine all angular spread values to a single variable
    AS = cat( 3, ASD,ASA,ESD,ESA );
    
    for i_ang = 1 : 4
        
        s = zeros( N,F ); % Scaling factor
        for i_freq = 1 : F
            
            % Current angular spread
            ang  = path_angles(:,:,i_ang);
            as   = qf.calc_angular_spreads( ang , pow(:,:,i_freq), 1 );
            
            % The following itertive optimization tries to find the optimal STD of the NLOS angles
            % such that the given angular spread is reached. This iteration is needed to correctly 
            % uncorporate the ground reflection.
            
            % Calculate the scaling coefficient as start values for an iterative optimization
            sL = 10*log10( AS(:,i_freq,i_ang) ./ as  );

            lp = 1;
            upd  = true( N,1 );         % The values that need to be updated
            while any( upd ) && lp < 10
                % Apply the scaling coefficients to the angles and update the angles and the angular spread
                ang( upd, Lf+1:end ) = (10.^(0.1*sL(upd))*oLn) .* path_angles( upd, Lf+1:end, i_ang );
                as = qf.calc_angular_spreads( ang(upd,:), pow(upd,:,i_freq),0 );
                
                % Calculate the new scaling coefficient
                sN = zeros( size( sL ) );
                sN(upd)  = 10*log10( AS(upd,i_freq,i_ang) ./ as  );
                
                upd = abs( sN ) > 1e-4;
                sL(upd) = sL(upd) + sN(upd);
                lp = lp + 1;
            end
            
            s(:,i_freq) = 10.^( 0.1*sL );
        end
        
        if F > 1 % Average scaling coefficient for multiple frequencies
            s = mean(s,2);
        end
        
        % The initial angles range from -pi/2 to pi/2, the initial angular spread is 30°
        % Due to the wrapping of the angles, thy cannot exceed
        if i_ang < 2.5 % Azimuth angles
            s( s > 3 ) = 3;
        else % Elevation angles
            s( s > 1.5 ) = 1.5;
        end
        
        % Scale NLOS angles
        path_angles(:,Lf+1:end,i_ang) = ( s*oLn ) .* path_angles(:,Lf+1:end,i_ang);
    end
    
    % Elevation angles outside the +/- pi/2 range would change the azimuth angles as well
    % We therefore fix mirror the elevation angles at the poles
    tmp = path_angles(:,:,[3,4]);
    ind = tmp > pi/2;
    tmp(ind) = pi-tmp(ind);
    ind = tmp < -pi/2;
    tmp(ind) = -pi-tmp(ind);
    path_angles(:,:,[3,4]) = tmp;
    
    % Apply the LOS angles
    % Wrap angles around the unit circle (does not change AS)
    path_angles  = mod( real(path_angles) + pi, 2*pi) - pi;
    
    % Transform departure angles to Cartesian coordinates
    path_angles = permute( path_angles, [3,2,1] );
    C = zeros( 3,L,N );
    C(1,:,:) = cos( path_angles(1,:,:) ) .* cos( path_angles(3,:,:) );
    C(2,:,:) = sin( path_angles(1,:,:) ) .* cos( path_angles(3,:,:) );
    C(3,:,:) = sin( path_angles(3,:,:) );
    
    % Build rotation matrix
    R = zeros( 3,3,N );
    ce = cos( mu_fixed(:,3) );
    se = sin( mu_fixed(:,3) );
    ca = cos( mu_fixed(:,1) );
    sa = sin( mu_fixed(:,1) );
    R(1,1,:) = ce.*ca;
    R(1,2,:) = -sa;
    R(1,3,:) = -se.*ca;
    R(2,1,:) = ce.*sa;
    R(2,2,:) = ca;
    R(2,3,:) = -se.*sa;
    R(3,1,:) = se;
    R(3,2,:) = 0;
    R(3,3,:) = ce;
    
    % Apply rotation in Cartesian coordinates
    for n = 1 : N
        C(:,:,n) = R(:,:,n) * C(:,:,n);
    end
    
    % Calculate rotated departure angles
    hypotxy = hypot( C(1,:,:),C(2,:,:) );
    EoD = atan2(C(3,:,:),hypotxy);
    AoD = atan2(C(2,:,:),C(1,:,:));
    EoD = permute( EoD, [3,2,1] );
    AoD = permute( AoD, [3,2,1] );
    
    % Transform arrival angles to Cartesian coordinates
    C = zeros( 3,L,N );
    C(1,:,:) = cos( path_angles(2,:,:) ) .* cos( path_angles(4,:,:) );
    C(2,:,:) = sin( path_angles(2,:,:) ) .* cos( path_angles(4,:,:) );
    C(3,:,:) = sin( path_angles(4,:,:) );
    
    % Build rotation matrix
    R = zeros( 3,3,N );
    ce = cos( mu_fixed(:,4) );
    se = sin( mu_fixed(:,4) );
    ca = cos( mu_fixed(:,2) );
    sa = sin( mu_fixed(:,2) );
    R(1,1,:) = ce.*ca;
    R(1,2,:) = -sa;
    R(1,3,:) = -se.*ca;
    R(2,1,:) = ce.*sa;
    R(2,2,:) = ca;
    R(2,3,:) = -se.*sa;
    R(3,1,:) = se;
    R(3,2,:) = 0;
    R(3,3,:) = ce;
    
    % Apply rotation in Cartesian coordinates
    for n = 1 : N
        C(:,:,n) = R(:,:,n) * C(:,:,n);
    end
    
    % Calculate rotated departure angles
    hypotxy = hypot( C(1,:,:),C(2,:,:) );
    EoA = atan2(C(3,:,:),hypotxy);
    AoA = atan2(C(2,:,:),C(1,:,:));
    EoA = permute( EoA, [3,2,1] );
    AoA = permute( AoA, [3,2,1] );
end

end
