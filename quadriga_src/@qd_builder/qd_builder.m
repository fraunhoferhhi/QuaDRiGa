classdef qd_builder < handle
%QD_BUILDER Class for generating the channel coefficients and the model pareters
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

properties
    name = '';                  % Name of the parameter set object
end

properties(Dependent)
    scenario                    % The name of the scenario (text string)
    scenpar                     % The parameter table
end

properties
    plpar = [];                 % Parameters for the path loss (scenario-dependent)
    
    simpar = qd_simulation_parameters;    	% Object of class qd_simulation_parameters
    tx_array = [];              % Handles of qd_arrayant objects for each Tx
    rx_array = [];              % Handles of qd_arrayant objects for each Rx
    tx_track = [];              % Handles of Track objects for each Tx
    rx_track = [];              % Handles of Track objects for each Rx
    
    % The transmitter position obtained from the corresponding 'layout.tx_position'
    tx_position = [];
    
    % The list of initial positions for which LSPs are generated
    %   This variable is obtained from 'track.initial_position' and 'layout.rx_position'
    rx_positions = [];
    
    sos = [];                  	% The large-scale parameter SOS generators
    gr_sos = [];                % The SOS generator for the ground reflection coefficient
    path_sos = [];              % The SOS generators for the generation of MPCs
    xpr_sos = [];               % The SOS generators for the generation of the linear NLOS polarization
    pin_sos = [];               % The SOS generators for the generation of initial phases
    clst_dl_sos = [];           % The SOS generators for the generation per-cluster delay offsets
    
    ds = [];                    % The RMS delay spread in [s] for each receiver position
    kf = [];                    % The Rician K-Factor [linear scale] for each receiver position
    sf = [];                    % The shadow fading [linear scale] for each receiver position
    asD = [];                   % The azimuth spread of departure in [deg] for each receiver position
    asA = [];                   % The azimuth spread of arrival in [deg] for each receiver position
    esD = [];                   % The elevation spread of departure in [deg] for each receiver position
    esA = [];                   % The elevation spread of arrival in [deg] for each receiver position
    xpr  = [];                  % The cross polarization ratio [linear scale] for each receiver position
    gr_epsilon_r = [];          % The relative permittivity for the ground reflection
    
    NumClusters                 % The number of clusters.
    NumSubPaths                 % The number of sub-paths per cluster
    
    % The initial delays for each path in [s]. Rows correspond to the
    % MTs, columns to the paths.
    taus
    
    % The normalized initial power (squared average amplitude) for each path. Rows correspond to the
    % MT, columns to the paths. The sum over all columns must be 1.
    pow
    
    AoD                         % The initial azimuth of departure angles for each path in [rad].
    AoA                         % The initial azimuth of arrival angles for each path in [rad].
    EoD                         % The initial elevation of departure angles for each path in [rad].
    EoA                         % The initial elevation of departure angles for each path in [rad].
    
    % The polarization rotation angle for the lineaar XPR in [rad]. 
    % For 3GPP baseline simulations, this property stores the per-path XPR (linear units).
    gamma
    
    % The phase offset angle for the circular XPR in [rad]. 
    % For 3GPP baseline simulations, this property stores the initial random phases. The dimensions
    % correspond to polarization matrix index  '[ 1 3 ; 2 4 ]', the subpath number and the MT.
    kappa
    
    pin                         % The initial phases in [rad] for each sub-path.
    
    % A random index list for the mutual coupling of subpaths at the Tx
    % and Rx. The dimensions correspond to the subpath index (1-20),
    % the angle (AoD, AoA, EoD, EoA), the path number and the MT.
    subpath_coupling
    
    fbs_pos = [];               % The positions of the first-bounce scatterers
    lbs_pos = [];               % The positions of the last-bounce scatterers
    
    % Indicates if the builder is a dual-mobility builder
    dual_mobility = -1;
end

properties(Dependent,SetAccess=protected)
    % Number of receiver positions associated to this 'parameter_set' object
    %   Note that each segment in longer tracks is considered a new Rx
    %   position.
    no_rx_positions
    
    lsp_vals                    % The distribution values of the LSPs
    lsp_xcorr_chk               % Indicator if cross-correlation matrix is positive definite
end

properties(Dependent)
    lsp_xcorr                   % The cross-correlation matrix for the LSPs
end

properties(Dependent,Hidden)
    scenpar_nocheck
end

% Data storage
properties(Access=private)
    Pscenario               = '';
    Pscenpar                = [];
    pow_wo_kf
end

properties(Hidden)
    OctEq = false; % For qf.eq_octave
end

methods
    
    % Constructor
    function h_builder = qd_builder( scenario, check_parfiles, scenpar )
        if exist( 'scenpar' , 'var' ) && ~isempty( scenpar )
            % If scenpar is given, we skip all checks, assuming that the values are correct
            h_builder.Pscenpar = scenpar;
            
        elseif exist( 'scenario' , 'var' ) && ~isempty( scenario )
            if exist( 'check_parfiles' , 'var' )
                if ~( all(size(check_parfiles) == [1 1]) ...
                        && (isnumeric(check_parfiles) || islogical(check_parfiles)) ...
                        && any( check_parfiles == [0 1] ) )
                    error('QuaDRiGa:qd_builder:WrongInput','??? "check_parfiles" must be 0 or 1')
                end
            else
                check_parfiles = true;
            end
            set_scenario_table( h_builder, scenario, check_parfiles );
        end
        
        % Initialize simpar with default settings
        h_builder.simpar = qd_simulation_parameters;
    end
    
    % Get functions
    function out = get.scenario(h_builder)
        out = h_builder.Pscenario;
    end
    function out = get.scenpar(h_builder)
        out = h_builder.Pscenpar;
    end
    function out = get.scenpar_nocheck(h_builder)
        out = h_builder.Pscenpar;
    end
    function out = get.no_rx_positions(h_builder)
        out = size( h_builder.rx_positions,2 );
    end
    function out = get.lsp_vals(h_builder)
        if isempty( h_builder.Pscenpar )
            out = [];
        else
            % Carrier frequency in GHz
            f_GHz = h_builder.simpar.center_frequency / 1e9;
            oF = ones( 1,numel( f_GHz ));
            o8 = ones( 8,1 );
            
            scenpar = h_builder.Pscenpar;
            
            % Reference frequency offset (omega)
            omega = [ scenpar.DS_omega; scenpar.KF_omega; scenpar.SF_omega; scenpar.AS_D_omega; ...
                scenpar.AS_A_omega; scenpar.ES_D_omega; scenpar.ES_A_omega; scenpar.XPR_omega ];
            
            % Reference values (mu) including the frequency scaling
            mu = [ scenpar.DS_mu; scenpar.KF_mu; 0; scenpar.AS_D_mu; scenpar.AS_A_mu;...
                scenpar.ES_D_mu; scenpar.ES_A_mu; scenpar.XPR_mu ];
            gammap = [ scenpar.DS_gamma;scenpar.KF_gamma;0; scenpar.AS_D_gamma; ...
                scenpar.AS_A_gamma ;scenpar.ES_D_gamma ;scenpar.ES_A_gamma; scenpar.XPR_gamma];
            mu = mu( :,oF ) + gammap(:,oF) .* log10( omega(:,oF) + f_GHz(o8,:) );
            
            % STD of the LSPs (sigma) including the frequency scaling
            sigma = [ scenpar.DS_sigma;scenpar.KF_sigma;scenpar.SF_sigma; scenpar.AS_D_sigma;...
                scenpar.AS_A_sigma;scenpar.ES_D_sigma;scenpar.ES_A_sigma; scenpar.XPR_sigma ];
            delta = [ scenpar.DS_delta;scenpar.KF_delta;scenpar.SF_delta; scenpar.AS_D_delta;...
                scenpar.AS_A_delta;scenpar.ES_D_delta;scenpar.ES_A_delta; scenpar.XPR_delta ];
            sigma = sigma( :,oF ) + delta(:,oF) .* log10( omega(:,oF) + f_GHz(o8,:) );
            sigma( sigma<0 ) = 0;
            
            % Decorr dist (lambda)
            lambda = [ scenpar.DS_lambda;scenpar.KF_lambda;scenpar.SF_lambda;...
                scenpar.AS_D_lambda;scenpar.AS_A_lambda;scenpar.ES_D_lambda;...
                scenpar.ES_A_lambda; scenpar.XPR_lambda ];
            lambda = lambda( :,oF );
            
            % Distance-dependence (epsilon) of the reference value
            epsilon = [ scenpar.DS_epsilon; scenpar.KF_epsilon; 0;...
                scenpar.AS_D_epsilon; scenpar.AS_A_epsilon; scenpar.ES_D_epsilon;...
                scenpar.ES_A_epsilon; scenpar.XPR_epsilon ];
            epsilon = epsilon( :,oF );
            
            % Height-dependence (zeta) of the reference value
            zeta = [ scenpar.DS_zeta; scenpar.KF_zeta; 0;...
                scenpar.AS_D_zeta; scenpar.AS_A_zeta; scenpar.ES_D_zeta;...
                scenpar.ES_A_zeta; scenpar.XPR_zeta ];
            zeta = zeta( :,oF );
            
            % Elevation-dependence (alpha) of the reference value
            alpha = [ scenpar.DS_alpha; scenpar.KF_alpha; 0;...
                scenpar.AS_D_alpha; scenpar.AS_A_alpha; scenpar.ES_D_alpha;...
                scenpar.ES_A_alpha; scenpar.XPR_alpha ];
            alpha = alpha( :,oF );
            
            % Distance-dependence (kappa) of the reference STD
            kappap = [ scenpar.DS_kappa; scenpar.KF_kappa; scenpar.SF_kappa;...
                scenpar.AS_D_kappa; scenpar.AS_A_kappa; scenpar.ES_D_kappa;...
                scenpar.ES_A_kappa; scenpar.XPR_kappa ];
            kappap = kappap( :,oF );
            
            % Height-dependence (tau) of the reference STD
            tau = [ scenpar.DS_tau; scenpar.KF_tau; scenpar.SF_tau;...
                scenpar.AS_D_tau; scenpar.AS_A_tau; scenpar.ES_D_tau;...
                scenpar.ES_A_tau; scenpar.XPR_tau ];
            tau = tau( :,oF );
            
            % Elevation-dependence (beta) of the reference STD
            beta = [ scenpar.DS_beta; scenpar.KF_beta; scenpar.SF_beta;...
                scenpar.AS_D_beta; scenpar.AS_A_beta; scenpar.ES_D_beta;...
                scenpar.ES_A_beta; scenpar.XPR_beta ];
            beta = beta( :,oF );
            
            % Assemble output
            out = [mu(:), sigma(:), lambda(:), epsilon(:), zeta(:), alpha(:), kappap(:), tau(:), beta(:)];
            out = permute( reshape( out, [],numel( f_GHz ),9 ),[1,3,2] );
        end
    end
    function out = get.lsp_xcorr(h_builder)
        if isempty( h_builder.Pscenpar )
            out = [];
        else
            value = h_builder.Pscenpar;
            a = value.ds_kf;          % delay spread vs k-factor
            b = value.ds_sf;          % delay spread vs shadowing std
            c = value.asD_ds;         % departure AS vs delay spread
            d = value.asA_ds;         % arrival AS vs delay spread
            e = value.esD_ds;         % departure ES vs delay spread
            f = value.esA_ds;         % arrival ES vs delay spread
            g = value.sf_kf;          % shadowing std vs k-factor
            h = value.asD_kf;         % departure AS vs k-factor
            k = value.asA_kf;         % arrival AS vs k-factor
            l = value.esD_kf;         % departure DS vs k-factor
            m = value.esA_kf;         % arrival DS vs k-factor
            n = value.asD_sf;         % departure AS vs shadowing std
            o = value.asA_sf;         % arrival AS vs shadowing std
            p = value.esD_sf;         % departure ES vs shadowing std
            q = value.esA_sf;         % arrival ES vs shadowing std
            r = value.asD_asA;        % departure AS vs arrival AS
            s = value.esD_asD;        % departure ES vs departure AS
            t = value.esA_asD;        % arrival ES vs departure AS
            u = value.esD_asA;        % departure ES vs arrival AS
            v = value.esA_asA;        % arrival ES vs arrival AS
            w = value.esD_esA;        % departure ES vs arrival ES
            x1 = value.xpr_ds;        % xpr vs delay spread
            x2 = value.xpr_kf;        % xpr vs k-factor
            x3 = value.xpr_sf;        % xpr vs shadowing
            x4 = value.xpr_asd;       % xpr vs departure AS
            x5 = value.xpr_asa;       % xpr vs arrival AS
            x6 = value.xpr_esd;       % xpr vs departure ES
            x7 = value.xpr_esa;       % xpr vs arrival ES

            out = [ 1  a  b  c  d  e  f x1;...
                a  1  g  h  k  l  m x2;...
                b  g  1  n  o  p  q x3;...
                c  h  n  1  r  s  t x4;...
                d  k  o  r  1  u  v x5;...
                e  l  p  s  u  1  w x6;...
                f  m  q  t  v  w  1 x7;...
                x1  x2  x3  x4  x5  x6  x7 1];
        end
    end
    function out = get.lsp_xcorr_chk(h_builder)
        out = false;
        if ~isempty( h_builder.Pscenpar )
            [~,p] = chol( h_builder.lsp_xcorr, 'lower');
            if p <= 0
                out = true;
            end
        end
    end
    
    % Set functions
    function set.scenario(h_builder,value)
        if ~( ischar(value) )
            error('QuaDRiGa:qd_builder:WrongInput','??? "scenario" must be a string.')
        end
        set_scenario_table( h_builder, value);
    end
    
    function set.scenpar(h_builder,value)
        if ~( isstruct(value) )
            error('QuaDRiGa:qd_builder:WrongInput','??? "scenpar" must be a structure.')
        end
        set_scenario_table( h_builder, value );
    end
    
    function set.scenpar_nocheck(h_builder,value)
        % Faster when we know that "scenpar" is correct
        h_builder.Pscenpar = value;
    end
    
    function set.lsp_xcorr(h_builder,value)
        if any( size( value ) ~= 8 )
            error('QuaDRiGa:qd_builder:WrongInput','??? "lsp_xcorr" must be a 8x8 matrix.')
        end
        sp = h_builder.Pscenpar;
        sp.ds_kf = value( 1,2 ) ;
        sp.ds_sf = value( 1,3 ) ;
        sp.asD_ds = value( 1,4 ) ;
        sp.asA_ds = value( 1,5 ) ;
        sp.esD_ds = value( 1,6 ) ;
        sp.esA_ds = value( 1,7 ) ;
        sp.xpr_ds = value( 1,8 ) ;
        sp.sf_kf = value( 2,3 ) ;
        sp.asD_kf = value( 2,4 ) ;
        sp.asA_kf = value( 2,5 ) ;
        sp.esD_kf = value( 2,6 ) ;
        sp.esA_kf = value( 2,7 ) ;
        sp.xpr_kf = value( 2,8 ) ;
        sp.asD_sf = value( 3,4 ) ;
        sp.asA_sf = value( 3,5 ) ;
        sp.esD_sf = value( 3,6 ) ;
        sp.esA_sf = value( 3,7 ) ;
        sp.xpr_sf = value( 3,8 ) ;
        sp.asD_asA = value( 4,5 ) ;
        sp.esD_asD = value( 4,6 ) ;
        sp.esA_asD = value( 4,7 ) ;
        sp.xpr_asd = value( 4,8 ) ;
        sp.esD_asA = value( 5,6 ) ;
        sp.esA_asA = value( 5,7 ) ;
        sp.xpr_asa = value( 5,8 ) ;
        sp.esD_esA = value( 6,7 ) ;
        sp.xpr_esd = value( 6,8 ) ;
        sp.xpr_esa = value( 7,8 ) ;
        h_builder.Pscenpar = sp;
    end
end

methods(Static)
    [ scenarios , config_folder , file_names ] = supported_scenarios( parse_shortnames )
    varargout = call_private_fcn( functionName, varargin )
end

methods % Legacy parameter generation functions
    function gen_lsf_parameters( h_builder , force )   % Legacy "gen_lsf_parameters"
        if exist('force','var') && force(1)
            gen_parameters( h_builder, 0 );  % Reset all
        end
        gen_parameters( h_builder, 1 );
    end
    function gen_ssf_parameters( h_builder )  % Legacy "gen_ssf_parameters"
        gen_parameters( h_builder );
    end
end
end
