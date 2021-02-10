function h_channel = get_channels( h_builder, vb_dots )
%GET_CHANNELS Calculate the channel coefficients
%
% Calling object:
%   Object array
%
% Output:
%   h_channel
%   A vector of qd_channel objects
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

% Array indexing is needed for Octave
verbose = h_builder(1,1,1,1).simpar.show_progress_bars;
if verbose && nargin == 1
    fprintf('Channels     [');
    vb_dots = 50;
    tStart = clock;
end
m0=0;

if numel(h_builder) > 1
    
    % Equally distribute the dots in the progress bar
    sic = size( h_builder );
    vb_dots = zeros( 1,numel(h_builder) );
    for i_cb = 1 : numel(h_builder)
        [ i1,i2,i3,i4 ] = qf.qind2sub( sic, i_cb );
        if verbose
            vb_dots(i_cb) = h_builder(i1,i2,i3,i4).no_rx_positions;
        else
            % Workaround for Octave 4
            if numel( sic ) == 4
                h_builder(i1,i2,i3,i4).simpar.show_progress_bars = false;
            elseif numel( sic ) == 3
                h_builder(i1,i2,i3).simpar.show_progress_bars = false;
            else % 2 and 1
                h_builder(i1,i2).simpar.show_progress_bars = false;
            end
        end
    end
    if verbose
        vb_dots = init_progress_dots(vb_dots);
    end
    
    % Call each builder in the builder array and concatinate the output channels
    cnt = 1;
    h_channel = qd_channel;
    for i_cb = 1 : numel(h_builder)
        [ i1,i2,i3,i4 ] = qf.qind2sub( sic, i_cb );
        if h_builder( i1,i2,i3,i4 ).no_rx_positions > 0
            tmp = h_builder( i1,i2,i3,i4 ).get_channels( vb_dots(i_cb) );
            h_channel( 1, cnt : cnt+size(tmp,2)-1 ) = tmp;
            cnt = cnt + size(tmp,2);
        end
    end
    
else
    % Fix for octave 4.0 (conversion from object-array to single object)
    h_builder = h_builder(1,1);
    
    % Check if we have a single-frequency builder
    if numel( h_builder.simpar.center_frequency ) > 1
        error('QuaDRiGa:qd_builder:get_channels','get_channels only works for single-freqeuncy simulations.');
    end
    
    % Check if SSF parameters have been generated already
    if isempty( h_builder.taus )
        error('QuaDRiGa:qd_builder:get_channels','Small-scale fading parametes have not been generated yet.');
    end
    
    % Check if the builder is a dual-mobility builder and that the inputs are correctly formatted
    dual_mobility = h_builder.dual_mobility;
    if dual_mobility == -1
        h_builder.check_dual_mobility;
    end
    
    % Set initial parameters
    n_clusters      = h_builder.NumClusters;
    o_clusters      = ones(1,n_clusters,'uint8');
    n_paths         = sum(h_builder.NumSubPaths);
    n_subpaths      = h_builder.NumSubPaths;
    n_mobiles       = h_builder.no_rx_positions;
    o4              = ones(1,4);
    t2              = 2*ones(1,2);
    
    % These variables are often needed. Pre-computing them saves a lot of time
    use_3GPP_baseline = h_builder.simpar.use_3GPP_baseline; % logical
    use_ground_reflection = logical( h_builder.scenpar.GR_enabled ); % logical
    if use_3GPP_baseline && use_ground_reflection
        warning('QuaDRiGa:qd_builder:use_baseline_ground_reflection',...
            'Ground reflection is not supported by the 3GPP baseline model.');
        use_ground_reflection = false;
        warning('off','QuaDRiGa:qd_builder:use_baseline_ground_reflection');
    end
    wave_no = 2*pi/h_builder.simpar.wavelength;
    
    % Create new channel object
    h_channel = qd_channel;
    
    % The loop for each user position
    for i_mobile = 1 : n_mobiles
        if verbose; m1=ceil(i_mobile/n_mobiles*vb_dots); if m1>m0;
                for m2=1:m1-m0; fprintf('o'); end; m0=m1; end;
        end;
        
        % Check if the positions are correct
        if any( abs( h_builder.rx_track(1,i_mobile).initial_position - h_builder.rx_positions(:,i_mobile) ) > 1e-5 )
            if i_mobile > 1 && ~qf.eqo( h_builder.rx_track(1,i_mobile), h_builder.rx_track(1,1) )
                warning('QuaDRiGa:qd_builder:get_channels',...
                    'Rx position in track does not match the initial position, using initial position.');
            end
            h_builder.rx_track(1,i_mobile).initial_position = h_builder.rx_positions(:,i_mobile);
        end
        if any( abs( h_builder.tx_track(1,i_mobile).initial_position - h_builder.tx_position(:,i_mobile) ) > 1e-5 )
            if i_mobile > 1 && ~qf.eqo( h_builder.tx_track(1,i_mobile), h_builder.tx_track(1,1) )
                warning('QuaDRiGa:qd_builder:get_channels',...
                    'Tx position in track does not match the initial position, using initial position.');
            end
            h_builder.tx_track(1,i_mobile).initial_position = h_builder.tx_position(:,i_mobile);
        end
        
        % The array antenna interpolation is the most time intense operation in the channel builder.
        % We save some computing time by reading the array antennas here and passing them as
        % variables to the interpolate function later on.
        
        % Update Tx-array data
        if i_mobile == 1 || ~qf.eqo( h_builder.tx_array(1,i_mobile), h_builder.tx_array(1,i_mobile-1) )
            n_txant             = h_builder.tx_array(1,i_mobile).no_elements;
            [ tx_patV, tx_patH, tx_azimuth_grid, tx_elevation_grid, tx_element_pos ] = ...
                wrap_grid( h_builder.tx_array(1,i_mobile) );
            tx_patV = reshape( tx_patV, [], n_txant );
            tx_patH = reshape( tx_patH, [], n_txant );
        end
        
        % Update Rx-array data
        if i_mobile == 1 || ~qf.eqo( h_builder.rx_array(1,i_mobile), h_builder.tx_array(1,i_mobile-1) )
            n_rxant             = h_builder.rx_array(1,i_mobile).no_elements;
            [ rx_patV, rx_patH, rx_azimuth_grid, rx_elevation_grid, rx_element_pos ] = ...
                wrap_grid( h_builder.rx_array(1,i_mobile) );
            rx_patV = reshape( rx_patV, [], n_rxant );
            rx_patH = reshape( rx_patH, [], n_rxant );
        end
        
        % Read some commonly needed variables in order to save time.
        n_links     = n_rxant*n_txant;
        o_links     = ones(1,n_links,'uint8');
        n_snapshots = h_builder.rx_track(1,i_mobile).no_snapshots;
        o_snapshots = ones(1,n_snapshots,'uint8');
        initial_pos = h_builder.rx_track(1,i_mobile).segment_index( min( [h_builder.rx_track(1,i_mobile).no_segments,2] ));
        
        % Extract the random initial phases
        pin = h_builder.pin(i_mobile,:); % double
        
        if use_3GPP_baseline
            % If we don't use drifting and have a linear track, then the Doppler component is only
            % dependent on the rotating phases of the taps. So, we don't recalculate the antenna
            % response for each snapshot.
            
            % Get the angles of the 20 subpaths and perform random coupling.
            [ aod,eod,aoa,eoa,delay ] = get_subpath_angles( h_builder, i_mobile );
            
            % Calculate the distance-dependent phases
            lambda  = h_builder.simpar.wavelength;
            if h_builder.simpar.use_absolute_delays
                d_lms   = h_builder.simpar.speed_of_light * delay;
            else
                r       = h_builder.rx_positions(:,i_mobile) - h_builder.tx_position(:,i_mobile);
                norm_r  = sqrt(sum(r.^2)).';
                d_lms   = norm_r + h_builder.simpar.speed_of_light * delay;
            end
            phase   = 2*pi/lambda * mod(d_lms, lambda);
            phase   = clst_expand( phase, n_subpaths );
            
            % Doppler component
            % Without drifting, the Doppler component is calculated by plane wave approximation
            % using the distance from the initial position.
            tmp = h_builder.rx_track(1,i_mobile).positions;
            dist_rx = sqrt( sum([ tmp(1,:) - tmp(1,1) ; tmp(2,:) - tmp(2,1) ; tmp(3,:) - tmp(3,1)   ].^2 ) );
            no_snap_process = 1;       % Only process the first snapshot, everyting else is approximated
            update_tx_ant = true;      % Interpolate Tx array response
            
        else
            % Calculate the scatterer positions
            lbs_pos = h_builder.lbs_pos(:,:,i_mobile);
            fbs_pos = h_builder.fbs_pos(:,:,i_mobile);
            
            % Initialize the path delays
            delay = zeros( n_snapshots, n_clusters, n_rxant, n_txant );
            no_snap_process = n_snapshots;
        end
        
        % Travel directions
        rx_orientation = h_builder.rx_track(1,i_mobile).orientation;
        tx_orientation = h_builder.tx_track(1,i_mobile).orientation;
        if size( tx_orientation,2 ) == 1
            tx_orientation = tx_orientation(:,o_snapshots);
        end
        
        xprmat = zeros(4,n_paths);                         % Initialization of the pol. rotation
        if use_3GPP_baseline
            % This is the polarization coupling from WINNER / 3GPP. The polarization is initialized
            % by random phases which are scaled by the XPR. Polarization drifting is not supported.
            
            if size( h_builder.kappa,3 ) == 1
                % Parameters were generated using the advances non-baseline model
                % We need to update the polarization here to match the baseline model
                
                % Generate XPR, 3GPP TR 38.901, Step 9
                gamma = randn( 1, n_paths  ) * h_builder.lsp_vals(8,2) + h_builder.lsp_vals(8,1); % dB
                gamma = 10.^(0.1*gamma);        % Linear
                gamma(:,1) = Inf;               % LOS

                % Draw initial random phases, 3GPP TR 38.901, Step 10
                kappa = 2*pi*rand( 1, n_paths, 4 ) - pi;
                kappa(:,1,:) = 0;               % LOS
            else
                gamma = h_builder.gamma(i_mobile,:);
                kappa = h_builder.kappa(i_mobile,:,:);
            end
            
            % sqrt(1/XPR)-Values
            xpr = sqrt( reshape( 1./gamma ,1,n_paths ) );
            
            % Random initial phases
            xprmat = exp( 1j*permute( kappa , [3,2,1] ));
            
            % Global XPR-Matrix
            xprmat = xprmat .* [ones(1,n_paths);xpr;xpr;ones(1,n_paths)];
            
            % Identity Matrix for the LOS-Path
            xprmat( 1,1 ) = 1;
            xprmat( 2,1 ) = 0;
            xprmat( 3,1 ) = 0;
            xprmat( 4,1 ) = -1;
            
            % Set rotation angles to 0
            gamma = zeros( 1,n_paths);

        else
            % Conversion of the XPR into rotation angle for linear and circular polarization
            gamma = h_builder.gamma(i_mobile,:);
            kappa = exp(1j*h_builder.kappa(i_mobile,:));
        end
        
        % Placeholder for the coefficient calculation
        cn    = zeros( n_links , n_clusters , n_snapshots );
        
        % Placeholder for the radiated power
        ppat  = zeros( n_links , n_clusters , n_snapshots );
        
        % Do for each snapshot
        for i_snapshot = 1 : no_snap_process          % Track positions
            
            c  = zeros( n_links , n_paths );          % The pattern coefficient matrix
            cp = zeros( n_links , n_paths );          % The phase coefficient matrix
            
            % Update the drifting angles, phases and delays.
            if ~use_3GPP_baseline
                if i_snapshot == 1 % Initialize
                    [ aod, eod, aoa, eoa, phase, delay(i_snapshot,:,:,:),...
                        aod_los, eod_los, aoa_los, eoa_los, theta_r, update_tx_ant ] =...
                        update_drifting( h_builder(1,1), i_snapshot, i_mobile, fbs_pos, lbs_pos );
                else % Update
                    [ aod, eod, aoa, eoa, phase, delay(i_snapshot,:,:,:),...
                        aod_los, eod_los, aoa_los, eoa_los, theta_r, update_tx_ant ] =...
                        update_drifting( h_builder(1,1), i_snapshot );
                end
            end
            
            % Include the direction on travel in the angles
            [ ~, aoa_c, eoa_c, deg_a ] = qf.calc_ant_rotation( rx_orientation(3,i_snapshot),...
                -rx_orientation(2,i_snapshot), rx_orientation(1,i_snapshot), aoa, eoa );
            if update_tx_ant
                [ ~, aod_c, eod_c, deg_d ] = qf.calc_ant_rotation( tx_orientation(3,i_snapshot),...
                    -tx_orientation(2,i_snapshot), tx_orientation(1,i_snapshot), aod, eod );
            end
            
            % Apply the additional polarization rotation for the NLOS paths
            deg = deg_a - gamma( :,:, ones(1,size(aoa,3)) );
            
            % Calculate the additional polarization scaling factors for the ground reflection
            if use_ground_reflection
                epsilon_r = h_builder.gr_epsilon_r( i_mobile );   % Relative permittivity
                
                Z      = sqrt( epsilon_r - (cos(theta_r)).^2 );
                R_par  = (epsilon_r * sin(theta_r) - Z) ./ (epsilon_r * sin(theta_r) + Z);
                R_per  = ( sin(theta_r) - Z) ./ ( sin(theta_r) + Z);
                
                % Read the path power scaling that was used in "generate_initial_paths.m"
                P_LOS = h_builder.pow(i_mobile,1);
                P_GR  = h_builder.pow(i_mobile,2);
                Rsq   = 2 * P_GR / (P_LOS + P_GR);
                
                % Compensate for the power scaling in "generate_initial_paths.m"
                Sl = 1 / sqrt( 1-Rsq/2 );       % LOS path
                Sv = sqrt(2/Rsq) * R_par;       % GR path vertical pol.
                Sh = sqrt(2/Rsq) * R_per;       % GR path horizontal pol.
            end
            
            if use_3GPP_baseline % Planar waves
                % Interpolate the receive antenna patterns
                [Vr,Hr,Pr] = h_builder.rx_array(1,i_mobile).interpolate( aoa_c ,...
                    eoa_c, 1:n_rxant, rx_azimuth_grid, rx_elevation_grid ,...
                    rx_patV, rx_patH, rx_element_pos );
                Pr  = reshape( Pr,1,n_paths,n_rxant );
                
                % Interpolate the transmit antenna patterns
                [Vt,Ht,Pt] = h_builder.tx_array(1,1).interpolate( ...
                    aod_c, eod_c, 1:n_txant, tx_azimuth_grid, tx_elevation_grid, ...
                    tx_patV, tx_patH, tx_element_pos );
                Pt = reshape( Pt, 1 ,n_paths, n_txant );
                
                % Calculate the Doppler profile.
                doppler = reshape( cos(aoa_c+pi).*cos(eoa) ,1,n_paths );
                
                % Apply scaling for the GR path
                if use_ground_reflection
                    xprmat(:,1)     = Sl .* xprmat(:,1);
                    xprmat([1,2],2) = Sv .* xprmat([1,2],2);
                    xprmat([3,4],2) = Sh .* xprmat([3,4],2);
                end
                
            else % Spherical waves
                % Interpolate the receive antenna patterns (NLOS)
                Vr  = zeros(1,n_paths,n_rxant);
                Hr  = zeros(1,n_paths,n_rxant);
                for i_rx = 1 : n_rxant
                    [Vr(1,:,i_rx),Hr(1,:,i_rx)] =...
                        h_builder.rx_array(1,i_mobile).interpolate( ...
                        aoa_c(:,:,i_rx) , eoa_c(:,:,i_rx) , i_rx , ...
                        rx_azimuth_grid , rx_elevation_grid , rx_patV , ...
                        rx_patH, rx_element_pos);
                end
                
                % Interpolate the receive antenna patterns (NLOS)
                if update_tx_ant
                    Vt  = zeros(1,n_paths,n_txant);
                    Ht  = zeros(1,n_paths,n_txant);
                    for i_tx = 1 : n_txant
                        [Vt(1,:,i_tx),Ht(1,:,i_tx)] =...
                            h_builder.tx_array(1,i_mobile).interpolate( ...
                            aod_c(:,:,1,i_tx) , eod_c(:,:,1,i_tx) , i_tx , ...
                            tx_azimuth_grid , tx_elevation_grid , tx_patV , ...
                            tx_patH, tx_element_pos);
                    end
                end
                
                % Include the direction on travel in the LOS angles
                [ ~, aoa_los_c, eoa_los_c, deg_LOS_a ] = qf.calc_ant_rotation( rx_orientation(3,i_snapshot),...
                    -rx_orientation(2,i_snapshot), rx_orientation(1,i_snapshot), aoa_los, eoa_los );
                [ ~, aod_los_c, eod_los_c, deg_LOS_d ] = qf.calc_ant_rotation( tx_orientation(3,i_snapshot),...
                    -tx_orientation(2,i_snapshot), tx_orientation(1,i_snapshot), aod_los, eod_los );
                
                % Calculate the Rx-Array-Pattern LOS response for
                % each element separately.
                if use_ground_reflection
                    Vr_LOS  = zeros(1,2,n_rxant,n_txant);
                    Hr_LOS  = zeros(1,2,n_rxant,n_txant);
                    Vt_LOS  = zeros(1,2,n_rxant,n_txant);
                    Ht_LOS  = zeros(1,2,n_rxant,n_txant);
                else
                    Vr_LOS  = zeros(1,1,n_rxant,n_txant);
                    Hr_LOS  = zeros(1,1,n_rxant,n_txant);
                    Vt_LOS  = zeros(1,1,n_rxant,n_txant);
                    Ht_LOS  = zeros(1,1,n_rxant,n_txant);
                end
                
                for i_rx = 1 : n_rxant
                    [ Vr_LOS(1,:,i_rx,:), Hr_LOS(1,:,i_rx,:) ] =...
                        h_builder.rx_array(1,i_mobile).interpolate( ...
                        aoa_los_c(1,:,i_rx,:), eoa_los_c(1,:,i_rx,:), i_rx, ...
                        rx_azimuth_grid , rx_elevation_grid , rx_patV , ...
                        rx_patH, rx_element_pos);
                end
                
                for i_tx = 1 : n_txant
                    [ Vt_LOS(1,:,:,i_tx), Ht_LOS(1,:,:,i_tx) ] =...
                        h_builder.tx_array(1,1).interpolate( ...
                        aod_los_c(1,:,:,i_tx), eod_los_c(1,:,:,i_tx), i_tx, ...
                        tx_azimuth_grid, tx_elevation_grid, tx_patV,...
                        tx_patH, tx_element_pos );
                end
            end
            
            % Calculate the Tx polarization rotation matrix
            if update_tx_ant && ~use_3GPP_baseline
                xprmat_tx = zeros(2,n_paths,n_txant);
                for i_tx = 1 : n_txant
                    xprmat_tx(1,:,i_tx) =  cos( deg_d(1,:,1,i_tx) );
                    xprmat_tx(2,:,i_tx) =  sin( deg_d(1,:,1,i_tx) );
                end
            end
            
            % The main loop to calculate the channel coefficients
            for i_rx = 1 : n_rxant     % Rx elements
                
                % Rx uses spherical waves, we need to calculate one XPRMAT for each Rx antenna separately.
                if ~use_3GPP_baseline
                    % Calculate a common XPRMAT for all Tx and Rx antennas
                    xprmat_rx(1,:) =  cos( deg(1,:,i_rx) );
                    xprmat_rx(2,:) = -sin( deg(1,:,i_rx) );
                    xprmat_rx(3,:) =  xprmat_rx(2,:);
                    xprmat_rx(4,:) = -xprmat_rx(1,:);
                    
                    % Include circular phase offset
                    xprmat_rx([1,2],:) = xprmat_rx([1,2],:) .* conj(kappa([1,1],:));
                    xprmat_rx([3,4],:) = xprmat_rx([3,4],:) .* kappa([1,1],:);
                end
                
                for i_tx = 1 : n_txant                           % Transmit elements
                    ind = (i_tx-1)*n_rxant + i_rx;               % Index of element in c
                    
                    % Tx uses spherical waves, we need to update XPRMAT for each Tx antenna
                    if ~use_3GPP_baseline
                        % Calculate XPRMAT by combining TX and RX matrices
                        xprmat = xprmat_rx .* xprmat_tx(o4,:,i_tx);
                        xprmat([1,2],:) = xprmat([1,2],:) + xprmat_rx([3,4],:) .* xprmat_tx(t2,:,i_tx);
                        xprmat([3,4],:) = xprmat([3,4],:) - xprmat_rx([1,2],:) .* xprmat_tx(t2,:,i_tx);
                        
                        % Update XPRMAT for the LOS path
                        TC = -cos( deg_LOS_d(1,1,i_rx,i_tx) );
                        TS = -sin( deg_LOS_d(1,1,i_rx,i_tx) );
                        xprmat([1,2,4],1) = [ -TC,TS ; TS,TC ; TC,-TS ] *...
                            [ cos( deg_LOS_a(1,1,i_rx,i_tx) ) ; sin( deg_LOS_a(1,1,i_rx,i_tx) ) ];
                        xprmat(3,1) = xprmat(2,1);
                        
                        % Apply scaling for the GR path
                        if use_ground_reflection
                            TC = -cos( deg_LOS_d(1,2,i_rx,i_tx) );
                            TS = -sin( deg_LOS_d(1,2,i_rx,i_tx) );
                            SV = Sv(i_rx,i_tx);
                            SH = Sh(i_rx,i_tx);
                            xprmat(:,2) = [ -TC*SV,TS*SH ; TS*SH,TC*SV ; TS*SV,TC*SH ; TC*SH,-TS*SV ] *...
                                [ cos( deg_LOS_a(1,2,i_rx,i_tx) ) ; sin( deg_LOS_a(1,2,i_rx,i_tx) ) ];
                            xprmat(:,1) = Sl .* xprmat(:,1);
                        end
                    end
                    
                    % Get antenna responses
                    PatTx = [ Vt(1,:,i_tx) ; Ht(1,:,i_tx) ];
                    PatRx = [ Vr(1,:,i_rx) ; Hr(1,:,i_rx) ];
                    if ~use_3GPP_baseline
                        PatTx(1,1) = Vt_LOS(1,1,i_rx,i_tx);
                        PatTx(2,1) = Ht_LOS(1,1,i_rx,i_tx);
                        PatRx(1,1) = Vr_LOS(1,1,i_rx,i_tx);
                        PatRx(2,1) = Hr_LOS(1,1,i_rx,i_tx);
                        if use_ground_reflection
                            PatTx(1,2) = Vt_LOS(1,2,i_rx,i_tx);
                            PatTx(2,2) = Ht_LOS(1,2,i_rx,i_tx);
                            PatRx(1,2) = Vr_LOS(1,2,i_rx,i_tx);
                            PatRx(2,2) = Hr_LOS(1,2,i_rx,i_tx);
                        end
                    end
                    
                    % Get the channel coefficients without random phases
                    c(ind,:) = sum( [ sum( PatTx .* xprmat([1 3],:)) ;...
                        sum( PatTx .* xprmat([2 4],:))] .* PatRx );
                    
                    % The phases
                    if use_3GPP_baseline
                        % In drifting mode, we have to update the coefficient
                        % matrix with the time-variant Doppler profile.
                        cp(ind,:) = exp( -1j*( pin + wave_no*( Pt(1,:,i_tx) + Pr(1,:,i_rx) ) + phase( 1 , : ) ));
                    else
                        % The phases already contain the effect of the
                        % AoD. Hence, the parallel projection of the
                        % array antennas is not needed.
                        cp(ind,:) = exp( -1j*( pin + phase( 1 , : , i_rx , i_tx )));
                    end
                end
            end
            
            % There atr be random phases in the sub-paths of the antenna
            % patterns. This changes the power when summing up the coefficients.
            
            % Combine antenna patterns and phases
            ccp = c.*cp;
            ppat(:,:,i_snapshot) = clst_sum( abs(ccp).^2 , n_subpaths );
            cn(:,:,i_snapshot)   = clst_sum( ccp , n_subpaths );
        end
        
        if use_3GPP_baseline
            % Only one snapshot is calculated, the others are
            % emulated by phase rotation.
            
            % Combine pattern and phase for the first snapshopt
            c = c.*cp;
            
            for i_snapshot = 2 : n_snapshots
                % Generate rotating Dopplers for the sucessive snapshots
                cp = exp( -1j * wave_no * doppler * dist_rx(i_snapshot) );
                cp = cp( ones(1,n_links) , : );
                
                % Combine antenna patterns and phases
                ccp = c.*cp;
                ppat(:,:,i_snapshot) = clst_sum( abs(ccp).^2 , n_subpaths );
                cn(:,:,i_snapshot)   = clst_sum( ccp , n_subpaths );
            end
        end
        
        % The path powers
        p_cl = h_builder.pow(i_mobile*ones(1,n_links),: );
        
        % The powers in the current channel coefficients (complex sum)
        p_coeff = sum( abs(cn).^2, 3 ) ./ size(cn,3);
        
        % The powers of the antenna patterns at the given angles (power-sum)
        p_pat = sum( ppat,3 ) ./ size(ppat,3);
        
        % Correct the powers
        p_correct = sqrt( p_cl .* p_pat ./ p_coeff ./ n_subpaths(o_links,:) );
        p_correct( p_pat < 1e-30 ) = 0; % Fix NaN caused by 0/0
        cn = p_correct(:,:,ones(1,n_snapshots)) .* cn;
        
        % Now we apply the K-Factor and the shadowing profile
        if use_3GPP_baseline
            
            % Get the PL for the initial position only
            [ ~, ~, path_loss , scale_sf ] = h_builder.get_pl( h_builder.rx_track(1,i_mobile),...
                [],h_builder.tx_track(1,i_mobile) );
            rx_power = -path_loss + 10*log10( h_builder.sf(1,i_mobile) ) .* scale_sf;
            rx_power = sqrt( 10.^( 0.1 * rx_power ) );
            
            % The initial KF is already applied in path powers. Here,
            % we only need to apply the SF and the path loss.
            cn = cn * rx_power;
            
        else
            
            % Get shadowing profile along the track from the correlation map. The first vector is
            % the K-Factor and the second vector is the SF. The initial K-Factor is already applied
            % in the path powers. We thus need to correct for that factor.
            [sf,kf] = h_builder.get_sf_profile( h_builder.rx_track(1,i_mobile), h_builder.tx_track(1,i_mobile) );
            
            % Scaling factor for the KF
            kf  = kf./h_builder.kf(1,i_mobile);
            if use_ground_reflection
                p1 = h_builder.pow( i_mobile,1 ) + h_builder.pow( i_mobile,2 );
                p1 = p1 ./ sum( h_builder.pow(i_mobile,:) );
            else
                p1 = h_builder.pow( i_mobile,1 );
            end
            kf_power_scale = sqrt( 1+p1*(kf-1) );
            kf = sqrt(kf);
            
            % Calculate the path loss
            [ path_loss , scale_sf ] = h_builder.get_pl( h_builder.rx_track(1,i_mobile), [], h_builder.tx_track(1,i_mobile) );
            
            % Calculate the Rx power
            rx_power = -path_loss + 10*log10( sf ) .* scale_sf;
            rx_power = sqrt( 10.^( 0.1 * rx_power ) );
            
            o_tmp = ones(1,n_txant*n_rxant,'uint8');
            kf = permute( kf(o_tmp,:) , [1,3,2] );
            
            rx_power = rx_power ./ kf_power_scale;
            rx_power = permute( rx_power(o_tmp,:) , [1,3,2] );
            
            cn(:,1,:) = cn(:,1,:).*kf;
            if use_ground_reflection
                cn(:,2,:) = cn(:,2,:).*kf;
            end
            cn = cn.*rx_power(:,o_clusters,:);
        end
        
        % Apply antenna coupling
        Ct = h_builder.tx_array(1,1).coupling;
        Cr = h_builder.rx_array(1,i_mobile).coupling.';
        
        % Reshape objects
        cn = reshape( cn , n_rxant , n_txant , n_clusters , n_snapshots );
        
        % Apply the antenna coupling
        c = zeros( size(Cr,1) , size(Ct,2) , n_clusters , n_snapshots );
        for i_snapshot = 1:n_snapshots
            for i_cluster = 1:n_clusters
                c(:,:,i_cluster,i_snapshot) = Cr * cn(:,:,i_cluster,i_snapshot) * Ct;
            end
        end
        clear cn
        
        if use_3GPP_baseline
            h_channel(1,i_mobile) = qd_channel( c , delay' , initial_pos );
        else
            % When we use high precision, the delays on all elements are different. However, antenna
            % coupling will merge the coefficients of different antennas. This needs to be
            % considered by the delays too. The delays on different elements are weighted by the
            % powers in the coupling matrix.
            Cr_dl = zeros( size( Cr ));
            for i_rx = 1:size( Cr,1 )
                tmp = abs( Cr( i_rx , : ) ).^2;
                Cr_dl( i_rx , : ) = tmp./sum(tmp);
            end
            
            Ct_dl = zeros( size( Ct ));
            for i_tx = 1:size( Ct,2 )
                tmp = abs( Ct( : , i_tx ) );
                Ct_dl( : , i_tx ) = tmp./sum(tmp);
            end
            
            % Here, we scale the delays for each path by the coupling powers.
            delay = permute( delay, [3,4,2,1] );
            dl = zeros( size(Cr,1) , size(Ct,2) , n_clusters , n_snapshots );
            for i_snapshot = 1:n_snapshots
                for i_cluster = 1:n_clusters
                    dl(:,:,i_cluster,i_snapshot) = Cr_dl * delay(:,:,i_cluster,i_snapshot) * Ct_dl;
                end
            end
            h_channel(1,i_mobile) = qd_channel( c , dl , initial_pos );
        end
        clear c
        
        % Store the channel name
        % This is important because the merger uses the name string to connect the channels.
        channel_name = h_builder.name;
        if isempty( channel_name ) || isempty( regexp( channel_name , '_', 'once' ) )
            channel_name = 'Scen_*';
        end
        % The "*" is added when there are multiple tx positions in the builder
        tmp = regexp( channel_name , '\*' );
        if ~isempty( tmp )
            channel_name = [ channel_name(1:tmp-1), h_builder.tx_track(1,i_mobile).name ];
        end
        channel_name = [ channel_name ,'_', h_builder.rx_track(1,i_mobile).name ];
        h_channel(1,i_mobile).name = channel_name;
        h_channel(1,i_mobile).rx_position = h_builder.rx_track(1,i_mobile).positions_abs;
        h_channel(1,i_mobile).tx_position = h_builder.tx_track(1,i_mobile).positions_abs;
        h_channel(1,i_mobile).center_frequency = h_builder.simpar.center_frequency(1,1);
        
        % Save Additional LSF and SSF information
        clear par_struct
        if use_ground_reflection
            par_struct.has_ground_reflection = 1;
        end
        par_struct.ds_parset = h_builder.ds( i_mobile ); % [s]
        par_struct.kf_parset = 10*log10( h_builder.kf( i_mobile ) ); % [db]
        if use_3GPP_baseline
            par_struct.pg_parset = 10*log10( rx_power.^2 ); % [db]
        else
            par_struct.pg_parset = 10*log10( mean(rx_power(:)).^2 ); % [db]
            par_struct.pg = 10*log10(abs( reshape( mean(mean(rx_power,1),2) , 1,[] ) ).^2);
        end
        par_struct.asD_parset = h_builder.asD( i_mobile ); % [deg]
        par_struct.asA_parset = h_builder.asA( i_mobile ); % [deg]
        par_struct.esD_parset = h_builder.esD( i_mobile ); % [deg]
        par_struct.esA_parset = h_builder.esA( i_mobile ); % [deg]
        par_struct.XPR_parset = 10*log10( h_builder.xpr( i_mobile ) ); % [db]
        
        % Save the individual per-path values
        par_struct.AoD_cb = h_builder.AoD( i_mobile,: ) * 180/pi; % [deg]
        par_struct.AoA_cb = h_builder.AoA( i_mobile,: ) * 180/pi; % [deg]
        par_struct.EoD_cb = h_builder.EoD( i_mobile,: ) * 180/pi; % [deg]
        par_struct.EoA_cb = h_builder.EoA( i_mobile,: ) * 180/pi; % [deg]
        par_struct.pow_cb = h_builder.pow( i_mobile,: );          % [W]
        
        % Calculate the spreads at the output of the builder
        par_struct.ds_cb  = qf.calc_delay_spread( h_builder.taus( i_mobile,: ), h_builder.pow( i_mobile, : ) );
        par_struct.asD_cb = qf.calc_angular_spreads( h_builder.AoD( i_mobile,: ), h_builder.pow( i_mobile, : ) ) * 180/pi;
        par_struct.asA_cb = qf.calc_angular_spreads( h_builder.AoA( i_mobile,: ), h_builder.pow( i_mobile, : ) ) * 180/pi;
        par_struct.esD_cb = qf.calc_angular_spreads( h_builder.EoD( i_mobile,: ), h_builder.pow( i_mobile, : ) ) * 180/pi;
        par_struct.esA_cb = qf.calc_angular_spreads( h_builder.EoA( i_mobile,: ), h_builder.pow( i_mobile, : ) ) * 180/pi;
        
        if ~use_3GPP_baseline
            par_struct.NumSubPaths = h_builder.NumSubPaths;
            par_struct.fbs_pos = fbs_pos;
            par_struct.lbs_pos = lbs_pos;
        end
        
        h_channel(1,i_mobile).par = par_struct;
    end
end

% Fix for octave
if numel( h_channel ) == 1;
    h_channel = h_channel(1,1);
end

if verbose && nargin == 1
    fprintf('] %5.0f seconds\n',round( etime(clock, tStart) ));
end

end
