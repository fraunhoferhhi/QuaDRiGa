function set_fbs_lbs( h_builder, lbs_pos, fbs_pos, pow, freespaceLOS )
%SET_FBS_LBS Set predefined LBS and FBS positions
%
% Calling object:
%   Object array
%
% Description:
%   This function assigns LBS and FBS positions to the builder object. It calculated the delays and
%   angles from those positions and updates the path powers. When using this function, sub-paths
%   are disabled.
%
% Input:
%   lbs_pos
%   The positions of the last-bounce-scatterers (NLOS only), dimensions: [3 x L-1]
%
%   fbs_pos
%   The positions of the first-bounce-scatterers (NLOS only), dimensions: [3 x L-1]. If this
%   variable is not given, LBS and FBS are identical.
%
%   pow
%   The path powers (optional). If this variable is not given, path powers are calculated from PL,
%   SF and KF distributions in the builder. The power-differences of the NLOS paths follow an
%   exponential decay based on the path-length differences. All LSPs are updated based on the power
%   values.
%
%   freespaceLOS
%   If set to 1, the LOS path power is calculated with the equivalent power of a freespace path.
%   SF, and KF are updated accordingly. Default: 0
%
%
% QuaDRiGa Copyright (C) 2011-2020
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


if ~exist('lbs_pos','var') || isempty(lbs_pos)
    error('QuaDRiGa:qd_builder:set_fbs_lbs','"lbs_pos" must be given.');
end
nc = size(lbs_pos,2); % Number of NLOS clusters

if ~exist('fbs_pos','var') || isempty(fbs_pos)
    fbs_pos = lbs_pos;
elseif any( size(lbs_pos) ~= size(fbs_pos) )
    error('QuaDRiGa:qd_builder:set_fbs_lbs','"fbs_pos" and "" must have the same size.');
end

if ~exist('pow','var')
    pow = [];
end

if ~exist('freespaceLOS','var') || isempty(freespaceLOS)
    freespaceLOS = false;
else
    freespaceLOS = logical(freespaceLOS(1));
end

if numel(h_builder) > 1
    sic = size( h_builder );
    for i_cb = 1 : numel(h_builder)
        [ i1,i2,i3,i4 ] = qf.qind2sub( sic, i_cb );
        if h_builder( i1,i2,i3,i4 ).no_rx_positions > 0
            set_fbs_lbs( h_builder( i1,i2,i3,i4 ), lbs_pos, fbs_pos, pow, freespaceLOS );
        end
    end
else
    h_builder = h_builder(1,1);
    show_progress_bars = h_builder.simpar.show_progress_bars;
    h_builder.simpar.show_progress_bars = 0;
    
    if h_builder.dual_mobility == -1
        h_builder.check_dual_mobility;
    end
    
    Nrx = h_builder.no_rx_positions;
    Orx = ones(1,Nrx);
    Onc = ones(1,nc);
    
    if Nrx > 1 && size(lbs_pos,3)==1
        lbs_pos = lbs_pos(:,:,Orx);
        fbs_pos = fbs_pos(:,:,Orx);
    elseif Nrx ~= size(lbs_pos,3)
        error('QuaDRiGa:qd_builder:set_fbs_lbs','"lbs_pos" must match the number of MTs.');
    end
    
    h_builder.NumClusters = nc+1;                               % Set number of clusters
    h_builder.NumSubPaths = ones(1,nc+1);                       % One path per cluster
    
    % Update scenpar
    scenpar = h_builder.scenpar;
    scenpar.SC_lambda = 0;
    scenpar.NumClusters = h_builder.NumClusters;
    scenpar.NumSubPaths = 1;
    scenpar.PerClusterDS = 0;
    scenpar.PerClusterDS_gamma = 0;
    h_builder.scenpar = scenpar;
    %h_builder.lsp_xcorr = eye(8);
    
    % Generate values for SF and XPR
    h_builder.gen_parameters(0);
    h_builder.gen_parameters;
    
    % Get some common parameters
    ptx = permute( h_builder.tx_position, [1,3,2] );            % TX-position
    prx = permute( h_builder.rx_positions, [1,3,2] );           % RX-position
    d3d = permute( sqrt(sum(abs(ptx-prx).^2,1)) , [1,3,2] );    % 3D LOS distance
    ang = h_builder.get_angles;                                 % LOS angles
    
    % Set angles and delays from scatterer positions
    a = lbs_pos - prx(:,Onc,:);
    b = fbs_pos - ptx(:,Onc,:);
    da = sqrt( sum(abs( a ).^2,1) );
    db = sqrt( sum(abs( b ).^2,1) );
    dc = sqrt( sum(abs( lbs_pos - fbs_pos ).^2,1) );
    dpath = permute( (da+db+dc) ,[3,2,1] ) - d3d(Onc,:)';
    h_builder.AoA = [ang(2,:)'*pi/180,permute( atan2( a(2,:,:),a(1,:,:) ) , [3,2,1] )];
    h_builder.AoD = [ang(1,:)'*pi/180,permute( atan2( b(2,:,:),b(1,:,:) ) , [3,2,1] )];
    h_builder.EoA = [ang(4,:)'*pi/180,permute( atan2( a(3,:,:),sqrt(a(1,:,:).^2+a(2,:,:).^2) ) , [3,2,1] )];
    h_builder.EoD = [ang(3,:)'*pi/180,permute( atan2( b(3,:,:),sqrt(b(1,:,:).^2+b(2,:,:).^2) ) , [3,2,1] )];
    h_builder.taus = [ zeros(Nrx,1), dpath ] ./ h_builder.simpar.speed_of_light;
    
    pl = 10.^(-0.1*h_builder.get_pl);                           % Path Loss (linear)
    sf = h_builder.sf;                                      % Shadow Fading (linear)
    p_total = (pl.*sf)';                                    % Total received power
    
    % Freespace path powers derived from path lengths, no losses from scattering
    dpath = [d3d',permute( (da+db+dc) ,[3,2,1] )];
    p_free = 20*log10(dpath) + 32.45 + 20*log10(h_builder.simpar.center_frequency/1e9);
    p_free = 10.^(-p_free/10);
    
    isempty_pow = false;
    if isempty( pow )
        isempty_pow = true;
        pow = p_free;
    elseif size(pow,2) == nc && size(pow,1) == Nrx
        pow = [p_free(:,1),pow];
    elseif size(pow,2) == nc && size(pow,1) == 1
        pow = [p_free(:,1),pow(Orx,:)];
    elseif size(pow,2) == nc+1 && size(pow,1) == Nrx
        % pow = pow; :-)
    elseif size(pow,2) == nc+1 && size(pow,1) == 1
        pow = pow(Orx,:);
    else
        error('QuaDRiGa:qd_builder:set_fbs_lbs','Size of "pow" must match the number of MTs and paths.');
    end
    
    p_nlos = sum(pow(:,2:end),2);                           % Sum-power of the NLOS paths
    p_los  = h_builder.kf.'.*p_nlos;                        % LOS power calculated from KF
    scale  = p_total./(p_los+p_nlos);                       % Scaling coefficient
    
    if freespaceLOS                     % Use PL and SF from model as NLOS referecne, add deterministic LOS
        pow(:,2:end) = pow(:,2:end).*scale(:,Onc);          % Scale only NLOS paths, keep freespace LOS
    elseif isempty_pow                  % Use KF from builder to scale power values
        pow(:,1) = p_los;                                   % Assign LOS powers based on KF
        pow = pow.*scale(:,ones(1,nc+1));                   % Scale all paths to match PL and SF
    end
    
    h_builder.sos = [];
    
    h_builder.sf = sum(pow,2)'./pl;
    h_builder.pow = pow./(sum(pow,2)*ones(1,nc+1));
    
    h_builder.ds = qf.calc_delay_spread( h_builder.taus, h_builder.pow )';
    h_builder.kf = (pow(:,1)./sum(pow(:,2:end),2))';
    
    [ as, es ] = qf.calc_angular_spreads_sphere( h_builder.AoA, h_builder.EoA, h_builder.pow, 1 );
    h_builder.asA = as'*180/pi;
    h_builder.esA = es'*180/pi;
    
    [ as, es ] = qf.calc_angular_spreads_sphere( h_builder.AoD, h_builder.EoD, h_builder.pow, 1 );
    h_builder.asD = as'*180/pi;
    h_builder.esD = es'*180/pi;
    
    [ fbs_posB, lbs_posB ] = generate_fbs_lbs( h_builder.tx_position, h_builder.rx_positions,...
        h_builder.taus, h_builder.AoD, h_builder.AoA, h_builder.EoD, h_builder.EoA,...
        h_builder.NumSubPaths, h_builder.subpath_coupling, [0;0;0;0], 'legacy' );
    
    h_builder.lbs_pos = lbs_posB;
    h_builder.fbs_pos = fbs_posB;
    
    h_builder.simpar.show_progress_bars = show_progress_bars;
    
end

end

