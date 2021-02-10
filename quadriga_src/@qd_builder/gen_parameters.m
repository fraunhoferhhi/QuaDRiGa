function gen_parameters( h_builder, usage, vb_dots )
%GEN_PARAMETERS Generates LSF parameters, SSF parameters and scatterer positions
%
% Calling object:
%   Object array
%
% Description:
%   This function generates all parameters that are needed for the channel coefficient generation.
%   The outputs of the method are stored in the class properties. This includes the following
%   steps:
%
%   * Initialize the random generators by calling "qd_builder.init_sos". If the random generators
%     are already initialized, use the existing initialization.
%
%   * Generate correlated large scale parameters for each user position. Those parameters are
%     needed by the channel builder to calculate initial SSF parameters for each track or segment
%     which are then evolved into time varying channels.
%
%   * Generates the small-scale-fading parameters for the channel builder. Already existing
%     parameters are overwritten. However, due to the spatial consistency of the model, identical
%     values will be obtained for the same rx positions. Spatial consistency can be disabled by
%     setting "qd_builder.scenpar.SC_lambda = 0".* Calculates the positions of the scatterers.
%
%
% Input:
%   usage
%   Controls the working mode of the method. The allowed options are:
%
%   * usage = 0
%     Clears all exisiting LSF and SSF parameters including the SOS random generators.
%
%   * usage = 1
%     Generates only the LSF parameters. Exisiting LSF parameters will be overwritten and all SSF
%     parameters will be cleared. If the SOS generators are not initialized, they are initialized
%     first. Existing SOS generators in "qd_builder.sos" are reused. This leads to identical results
%     when calling the method multiple times.
%
%   * usage = 2
%     Generates the SSF parameters. Exisiting SSF parameters will be overwritten. Existing SOS
%     generators and LSF parameters will be reused.
%
%   * usage = 3
%     Calculates the scattering cluster positions from the exisiting SSF parameters. In some cases,
%     this may lead to changes in the departure angles (AoD, EoD). If LSF or SSF parameters are not
%     initialized, they are initialized first.
%
%   * usage = 4 (default)
%     Clears exisiting LSF parameters, SSF parameters and cluster positions and calculates new ones.
%     Existing SOS generators are reused.
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

if ~exist('usage','var') || isempty(usage)
    usage = 4;
end

verbose = h_builder(1,1).simpar.show_progress_bars;
if exist('vb_dots','var') && vb_dots == 0
    verbose = 0;
end

if verbose && nargin < 3
    fprintf('SSF Corr.    [');
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
    
    for i_cb = 1 : numel(h_builder)
        [ i1,i2,i3,i4 ] = qf.qind2sub( sic, i_cb );
        if h_builder( i1,i2,i3,i4 ).no_rx_positions > 0
            gen_parameters( h_builder( i1,i2,i3,i4 ), usage, vb_dots(i_cb) );
        end
    end
    
else
    % Fix for octave 4.0 (conversion from object-array to single object)
    h_builder = h_builder(1,1);
    
    % Check the validity of the input variables
    if h_builder.dual_mobility == -1
        h_builder.check_dual_mobility;
    end
    
    if usage == 0
        
        % Delete all existing values
        h_builder.ds = [];
        h_builder.kf = [];
        h_builder.sf = [];
        h_builder.asD = [];
        h_builder.asA = [];
        h_builder.esD = [];
        h_builder.esA = [];
        h_builder.xpr = [];
        h_builder.gr_epsilon_r = [];
        h_builder.taus = [];
        h_builder.pow = [];
        h_builder.AoD = [];
        h_builder.AoA = [];
        h_builder.EoD = [];
        h_builder.EoA = [];
        h_builder.gamma = [];
        h_builder.kappa = [];
        h_builder.pin = [];
        h_builder.fbs_pos = [];
        h_builder.lbs_pos = [];
        h_builder.NumClusters = [];
        h_builder.NumSubPaths = [];
        
        % Clear all SOS genreators and the subpath coupling
        h_builder.sos = [];
        h_builder.gr_sos = [];
        h_builder.path_sos = [];
        h_builder.xpr_sos = [];
        h_builder.pin_sos = [];
        h_builder.clst_dl_sos = [];
        h_builder.subpath_coupling = [];
        
    else
        
        % Initialize the random generators.
        % If they are already initialized, the exisiting values will not be changed.
        init_sos( h_builder );
        
        % Delete exisiting values
        switch usage
            case 2
                names = { 'taus','pow','AoD','AoA','EoD','EoA','gamma','pin','kappa',...
                    'gr_epsilon_r','fbs_pos','lbs_pos'};
                
            case 3
                names = { 'pin','fbs_pos','lbs_pos' };
                
            otherwise % case 1 or 4
                names = {'ds','kf','sf','asD','asA','esD','esA','xpr',...
                    'taus','pow','AoD','AoA','EoD','EoA','gamma','pin','kappa',...
                    'gr_epsilon_r','fbs_pos','lbs_pos'};
        end
        for n = 1 : numel( names )
            if ~isempty( h_builder.( names{n} ) )
                warning('QuaDRiGa:qd_builder:gen_ssf_parameters:exisitng',...
                    ['Exisitng parameters "',names{n},'" were deleted when calling "gen_parameters".']);
            end
            h_builder.(names{n}) = [];
        end
        
        % If there LSF parameters missing, generate them
        if usage == 2 && ( isempty(h_builder.ds) || isempty(h_builder.kf) || isempty(h_builder.sf) || ...
                isempty(h_builder.asD) || isempty(h_builder.asA) || isempty(h_builder.esD) || ...
                isempty(h_builder.esA) || isempty(h_builder.xpr) )
            h_builder.gen_parameters(1);
        end
        
        % If there SSF parameters missing, generate them
        if usage == 3 && ( isempty(h_builder.taus) || isempty(h_builder.pow) || isempty(h_builder.AoD) || ...
                isempty(h_builder.AoA) || isempty(h_builder.EoD) || isempty(h_builder.EoA) || ...
                isempty(h_builder.gamma) || isempty(h_builder.kappa) )
            h_builder.gen_parameters(2);
        end
        
        % Set the number of clusters
        n_clusters = h_builder.scenpar.NumClusters;
        
        % Get the number of frequencies
        n_freq = numel( h_builder.simpar.center_frequency );
        
        % Get the number of mobiles
        n_mobiles       = h_builder.no_rx_positions;
        
        use_3GPP_baseline = h_builder.simpar.use_3GPP_baseline;
        use_ground_reflection = logical( h_builder.scenpar.GR_enabled );
        if use_ground_reflection && n_clusters == 1
            error('QuaDRiGa:qd_builder:gen_parameters:gr_2_clusters',...
                'You need at least 2 clusters to enable the ground reflection option.');
        end
        if use_3GPP_baseline && use_ground_reflection
            warning('QuaDRiGa:qd_builder:use_baseline_ground_reflection',...
                'Ground reflection is not supported by the 3GPP baseline model.');
            use_ground_reflection = false;
            warning('off','QuaDRiGa:qd_builder:use_baseline_ground_reflection');
        end
        
        % Get the total number of paths
        if use_ground_reflection
            n_los_clusters = 2;
            n_nlos_clusters = n_clusters-2;
            n_subpaths = [1,1,ones(1,n_nlos_clusters) * h_builder.scenpar.NumSubPaths];
        else
            n_los_clusters = 1;
            n_nlos_clusters = n_clusters-1;
            n_subpaths = [1,ones(1,n_nlos_clusters) * h_builder.scenpar.NumSubPaths];
        end
        n_paths         = sum( n_subpaths );                                % Total number of paths
        n_LnM           = n_nlos_clusters * h_builder.scenpar.NumSubPaths;  % Total number of NLOS paths
        
        if h_builder.dual_mobility
            tx_position = h_builder.tx_position;
        else
            tx_position = h_builder.tx_position(:,1);
        end
        
        if usage == 1 || usage == 4
            % Generate LSF parameters
            if ~h_builder.lsp_xcorr_chk
                error('QuaDRiGa:qd_builder:gen_parameters:lsp_xcorr_chk',...
                    ['LSP cross-correlation matix of "',h_builder.name,'" is not positive definite.']);
            end
            [ ds, kf, sf, asD, asA, esD, esA, xpr ] = ...
                generate_lsf( tx_position, h_builder.rx_positions, h_builder.lsp_vals,...
                h_builder.lsp_xcorr, h_builder.sos,...
                h_builder.scenpar.ES_D_mu_A, h_builder.scenpar.ES_D_mu_min);
            
            % Save LSF parameters to builder object
            h_builder.ds  = ds;
            h_builder.kf  = kf;
            h_builder.sf  = sf;
            h_builder.asD = asD;
            h_builder.asA = asA;
            h_builder.esD = esD;
            h_builder.esA = esA;
            h_builder.xpr = xpr;
        else
            % Read LSF parameters from the builder object
            ds = h_builder.ds;
            kf = h_builder.kf;
            asD = h_builder.asD;
            asA = h_builder.asA;
            esD = h_builder.esD;
            esA = h_builder.esA;
            xpr = h_builder.xpr;
        end
        
        if usage == 2 || usage == 4
            
            % Set the number of clusters
            h_builder.NumClusters = n_clusters;
            
            % Set the number of sub-paths per cluster
            h_builder.NumSubPaths = n_subpaths;
         
            % Create input variables for the path generation function
            spreads = [ ds; asD; asA; esD; esA; kf ];
            spreads = reshape( spreads , n_freq , 6 , [] );
            spreads = permute( spreads, [2,3,1] );
            
            % Ground Reflection Parameters
            if use_ground_reflection
                if h_builder.scenpar.GR_epsilon == 0
                    gr_epsilon_r = generate_GR_parameters( tx_position, h_builder.rx_positions, ...
                        h_builder.simpar.center_frequency/1e9, h_builder.gr_sos );
                else
                    gr_epsilon_r = ones( n_freq, n_mobiles )*h_builder.scenpar.GR_epsilon;
                end
                h_builder.gr_epsilon_r = gr_epsilon_r;
            else
                gr_epsilon_r = [];
            end
            
            % Generate delays and angles
            if h_builder.simpar.use_3GPP_baseline
                [ pow, taus, AoD, AoA, EoD, EoA ] = generate_paths_baseline( h_builder.scenpar.NumClusters,...
                    tx_position, h_builder.rx_positions, spreads, h_builder.scenpar.r_DS, h_builder.scenpar.LNS_ksi, 1 );
            else
                [ pow, taus, AoD, AoA, EoD, EoA ] = generate_paths( h_builder.scenpar.NumClusters,...
                    tx_position, h_builder.rx_positions, spreads, gr_epsilon_r, h_builder.path_sos, h_builder.scenpar.r_DS );
            end
            
            % Assign values to the builder
            h_builder.pow  = pow;
            h_builder.taus = taus;
            h_builder.AoD  = AoD;
            h_builder.AoA  = AoA;
            h_builder.EoD  = EoD;
            h_builder.EoA  = EoA;
            
            if h_builder.simpar.use_3GPP_baseline
                % Generate XPR, 3GPP TR 38.901, Step 9
                randC       = randn( n_mobiles, n_LnM  );
                xpr_sigma   = h_builder.lsp_vals(8,2);
                xpr_mu      = h_builder.lsp_vals(8,1);
                xpr_nlos    = randC * xpr_sigma + xpr_mu; % dB
                los_xpr     = Inf(n_mobiles,1);
                h_builder.gamma = cat( 2, los_xpr, 10.^(0.1*xpr_nlos) );
                
                % Draw initial random phases, 3GPP TR 38.901, Step 10
                randC       = rand( n_mobiles, n_LnM, 4 );
                los_phase   = zeros( n_mobiles,1,4 );
                h_builder.kappa = cat( 2, los_phase, 2*pi*randC - pi );
                
            else
                % Use spatially consistent polarization rotation
                xpr_sigma = 0.1 * reshape( h_builder.lsp_vals(8,2,:), n_freq, 1 ) * ones(1,n_mobiles);
                [ gamma, kappa ] = generate_pol_rot( n_nlos_clusters, h_builder.scenpar.NumSubPaths,...
                    tx_position, h_builder.rx_positions, 10*log10(xpr), xpr_sigma, h_builder.xpr_sos );
                los_angle = zeros(n_mobiles,n_los_clusters,n_freq);
                h_builder.gamma = cat( 2, los_angle, gamma );
                h_builder.kappa = cat( 2, los_angle, kappa );
            end
        end
        
        % Update progress bar (50% point)
        if verbose; m1=ceil(1/2*vb_dots); if m1>m0;
                for m2=1:m1-m0; fprintf('o'); end; m0=m1; end;
        end;
        
        if usage == 3 || usage == 4
            
            % Generate random initial phases
            pin = zeros( n_mobiles,n_paths,n_freq );
            if h_builder.simpar.use_random_initial_phase
                for iF = 1 : n_freq  % Phases are independent for each frequency in multi-frequency simulations
                    
                    % Generate spatially correlated random variables for the NLOS phases
                    if isempty( h_builder.pin_sos )
                        randC = rand( n_mobiles,n_LnM );
                    else
                        if h_builder.dual_mobility
                            randC = val( h_builder.pin_sos(:,:,iF), h_builder.rx_positions,tx_position ).';   % Uniform
                        else
                            randC = val( h_builder.pin_sos(:,:,iF), h_builder.rx_positions ).';   % Uniform
                        end
                        randC = reshape( randC, n_mobiles, n_nlos_clusters, h_builder.scenpar.NumSubPaths );
                        randC = reshape( permute( randC, [1,3,2] ) , n_mobiles,n_LnM );
                    end
                    
                    % Set phases
                    pin(:,n_los_clusters+1:end,iF) = 2*pi*randC-pi;
                end
            end
            h_builder.pin = pin;
            
            % Apply PerClusterDS and AS
            if h_builder.NumClusters == h_builder.scenpar.NumClusters
                generate_subpaths( h_builder );
            end
            
            % Calculate scatterer positions
            if ~h_builder.simpar.use_3GPP_baseline
                PerClusterAS = [h_builder.scenpar.PerClusterAS_D; h_builder.scenpar.PerClusterAS_A;...
                    h_builder.scenpar.PerClusterES_D; h_builder.scenpar.PerClusterES_A ] * ones(1,n_freq);
                
                [ h_builder.fbs_pos, h_builder.lbs_pos, h_builder.AoD, h_builder.AoA, h_builder.EoD, h_builder.EoA ] = ...
                    generate_fbs_lbs( h_builder.tx_position, h_builder.rx_positions,...
                    h_builder.taus, h_builder.AoD, h_builder.AoA, h_builder.EoD, h_builder.EoA, ...
                    h_builder.NumSubPaths, h_builder.subpath_coupling, PerClusterAS, 'bounce2' );
            end
        end
    end
    
    % Update progress bar( 100% point )
    if verbose; m1=ceil((2/2*vb_dots)); if m1>m0;
            for m2=1:m1-m0; fprintf('o'); end; m0=m1; end;
    end;
end

if verbose && nargin < 3
    fprintf('] %5.0f seconds\n',round( etime(clock, tStart) ));
end

end
