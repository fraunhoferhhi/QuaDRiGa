function generate_subpaths( h_builder )
%GENERATE_SUBPATHS Generates the subpaths
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

persistent did_warn

% Number of deterministic (LOS) paths
if logical( h_builder.scenpar.GR_enabled ) && ~h_builder.simpar.use_3GPP_baseline
    nLOS = 2;
else
    nLOS = 1;
end

n_clusters = h_builder.NumClusters;             % Number of clusters
n_mobiles = h_builder.no_rx_positions;          % Number of mobiles
n_freq = numel( h_builder.simpar.center_frequency );

% ACF
acf = h_builder.simpar.autocorrelation_function;

% Use tx-position for dual-mobility scenarios
if h_builder.dual_mobility
    tx_position = h_builder.tx_position;
else
    tx_position = [];
end

% Get the per-cluster DS in seconds (it is given in ns in the config files)
if h_builder.scenpar.PerClusterDS_gamma ~= 0
    f_GHz = h_builder.simpar.center_frequency / 1e9;
    PerClusterDS = h_builder.scenpar.PerClusterDS_gamma * log10( f_GHz ) +...
        h_builder.scenpar.PerClusterDS;
    PerClusterDS = max( PerClusterDS, h_builder.scenpar.PerClusterDS_min ) / 1e9;
else
    PerClusterDS = h_builder.scenpar.PerClusterDS / 1e9;
end

switch h_builder.scenpar.SubpathMethod
    
    case 'legacy'
        % This is the legacy subpath generation method for bandwidths below 100 MHz. 
        % See: 3GPP TR 38.901 V14.1.0 (2017-06), pp36
        
        % Determine if clusters will be split in sub-clusters with different delays
        % This depends on 3 conditions:
        %   1. The per-cluster delay spread is set to values > 0 in the scenario definition
        %   2. The number of sub-paths is set to 20
        %   3. There is at least 1 NLOS path (exluding ground reflection)
        use_cluster_DS = h_builder.scenpar.PerClusterDS ~= 0 &...
            h_builder.scenpar.NumSubPaths == 20 &...
            n_clusters > nLOS;
        
        % Plot a warning if NumSubPaths is npt 20
        if h_builder.scenpar.PerClusterDS ~= 0 && h_builder.scenpar.NumSubPaths ~= 20
            if isempty( did_warn )
                disp(' ');
                warning('QuaDRiGa:builder:generate_subpaths:wrongNumSubPaths',...
                    'Number of subpaths must be 20 for legacy cluster splitting. PerClusterDS was not applied.');
                did_warn = true;
            end
        end
            
        % According the 3GPP 38.901 v14.1.0, Table 7.5-5, p37, the 2 strongest clusters should
        % be split in sub-clusters with a delay-offset that is determined by the
        % "PerClusterDS". However, this would break the spatial consistency sice the strongest
        % clusters change with locations. To solve this, ALL clusters are split into sub-clusters.
        if use_cluster_DS
            
            % Get the number of clusters that need to be split
            % nSPL = min( n_clusters-nLOS , 2 );  % Original version with 2 clusters
            nSPL = n_clusters-nLOS;
            splt = sort( [ 1:nLOS, (1:nSPL)+nLOS, (1:nSPL)+nLOS, (1:nSPL)+nLOS, nLOS+nSPL+1:n_clusters ] );
            
            if numel( PerClusterDS ) > 1        % Frequency-dependent !!!
                taus = h_builder.taus( :,splt, ones(1,n_freq) );
            else
                taus = h_builder.taus( :,splt );
            end
            
            pow  = h_builder.pow( :,splt,: );
            nsplt = [ 10 , 6 , 4 ];
            NumSubPaths_new = ones(1,nLOS);
            for n = 0:nSPL-1
                NumSubPaths_new = [ NumSubPaths_new, nsplt ];
                for m = 1:3
                    ii = nLOS+n*3+m;
                    if numel( PerClusterDS ) > 1        % Frequency-dependent !!!
                        for iF = 1:n_freq
                            taus( :,ii,iF ) = taus( :,ii,iF ) + PerClusterDS(iF)*1.28*(m-1);
                        end
                    else
                        taus( :,ii ) = taus( :,ii ) + PerClusterDS*1.28*(m-1);
                    end
                    pow( :,ii,: )  = pow( :,ii,: ).*nsplt(m)/20;
                end
            end
            NumSubPaths_new = [ NumSubPaths_new, h_builder.NumSubPaths( nLOS+nSPL+1:n_clusters ) ];
            
            h_builder.taus = taus;
            h_builder.pow = pow;
            h_builder.AoD = h_builder.AoD( :,splt );
            h_builder.AoA = h_builder.AoA( :,splt );
            h_builder.EoD = h_builder.EoD( :,splt );
            h_builder.EoA = h_builder.EoA( :,splt );
            h_builder.NumClusters = numel( splt );
            h_builder.NumSubPaths = NumSubPaths_new;
        end
        
    case 'mmMAGIC'
        % This is the mmMAGIC subpath generation method.
        % See: H2020-ICT-671650-mmMAGIC/D2.2, pp86, Section 4.5.3 Mapping of paths to sub-paths
        
        n_paths = sum( h_builder.NumSubPaths );
        r_DS = h_builder.scenpar.r_DS;                  % Delay scaling coefficient

        % Reserve Memory
        pow  = zeros( n_mobiles, n_paths, n_freq );
        if size( PerClusterDS,2) == 1
            taus = zeros( n_mobiles, n_paths );
        else
            taus = zeros( n_mobiles, n_paths, n_freq );
        end
        AoD = zeros( n_mobiles, n_paths );
        AoA = zeros( n_mobiles, n_paths );
        EoD = zeros( n_mobiles, n_paths );
        EoA = zeros( n_mobiles, n_paths );
        
        % Copy LOS and GR data
        pow(:,1:nLOS,:) = h_builder.pow(:,1:nLOS,:);  
        for n = 1:size( PerClusterDS,2)        
            taus(:,1:nLOS,n) = h_builder.taus(:,1:nLOS);  
        end
        AoD(:,1:nLOS) = h_builder.AoD(:,1:nLOS);
        AoA(:,1:nLOS) = h_builder.AoA(:,1:nLOS);
        EoD(:,1:nLOS) = h_builder.EoD(:,1:nLOS);
        EoA(:,1:nLOS) = h_builder.EoA(:,1:nLOS);
        
        st = nLOS + 1;
        for n = nLOS+1 : n_clusters
            M = h_builder.NumSubPaths(n);
            oM = ones( 1,M );
            
            % Generate spatially correlated random numbers for the delay offsets
            if isempty( h_builder.clst_dl_sos )
                randC = rand( n_mobiles, M );
            else
                randC = val( h_builder.clst_dl_sos( n-nLOS,: ), h_builder.rx_positions, tx_position ).';
            end
            
            % Genreate delay offsets
            powers = zeros( n_mobiles, M, size( PerClusterDS,2) );
            for iF = 1 : size( PerClusterDS,2)
                delays = -r_DS * PerClusterDS(1,iF) * log( randC );                     % Calc delays
                delays = delays - min( delays,[],2 ) * oM;                              % Normalize delays
                if PerClusterDS(1,iF) > 1e-10
                    powers(:,:,iF) = exp( -delays * (r_DS-1)/(r_DS*PerClusterDS(1,iF)) );   % Calc powers
                else
                    powers(:,:,iF) = ones( n_mobiles,M );
                end
                powers(:,:,iF) = powers(:,:,iF) ./ ( sum( powers(:,:,iF) , 2 ) * oM );  % Normalize powers
                ds = sqrt( sum(powers(:,:,iF).*delays.^2,2) - sum( powers(:,:,iF).*delays,2).^2 );
                if PerClusterDS(1,iF) > 1e-10 
                    delays = (PerClusterDS(1,iF)./ds) * oM .*  delays;                  % Scale delays
                end
                taus(:,st:st+M-1,iF) = delays + h_builder.taus(:,n)*oM;                 % Add cluster delay
            end
            
            % Duplicate powers if cluster-DS is not freq.-dependent
            if size( PerClusterDS,2) < n_freq
                powers = powers(:,:,ones(1,n_freq));
            end
            
            % Apply cluster powers
            for iF = 1:n_freq
                pow(:,st:st+M-1,iF) = powers(:,:,iF) .* (h_builder.pow(:,n,iF)*oM);
            end
            
            % Generate angle offsets (laplacian angular offsets)
            ao = 1:floor( (M+0.1)/2 );
            if ao(end) > M/2-0.1
                ao = [ log(2*ao/(M+1)), -log(2*ao(end:-1:1)/(M+1)) ];   % Even number of subpaths
            else
                ao = [ log(2*ao/M), 0, -log(2*ao(end:-1:1)/M) ];        % Odd number of subpaths
            end
            ao = ones(n_mobiles,1) * ao;                        % Per MT angle offsets in [deg]
            P = sum( powers,3 );                                % Average power offsets for all frequencies
            P = P ./ ( sum( P , 2 ) * oM );                     % Normalize powers
            as = sqrt( sum(P.*ao.^2,2) - sum( P.*ao,2).^2 );    % angular spread in [deg]
            
            AoD(:,st:st+M-1) = ((h_builder.scenpar.PerClusterAS_D./as) * oM) .* ao * pi/180 + h_builder.AoD(:,n)*oM;
            AoA(:,st:st+M-1) = ((h_builder.scenpar.PerClusterAS_A./as) * oM) .* ao * pi/180 + h_builder.AoA(:,n)*oM;
            EoD(:,st:st+M-1) = ((h_builder.scenpar.PerClusterES_D./as) * oM) .* ao * pi/180 + h_builder.EoD(:,n)*oM;
            EoA(:,st:st+M-1) = ((h_builder.scenpar.PerClusterES_A./as) * oM) .* ao * pi/180 + h_builder.EoA(:,n)*oM;
            st = st + M;
        end
        
        % Restrict angle ranges
        AoD = angle(exp(1j*AoD));
        AoA = angle(exp(1j*AoA));
        EoD( EoD >  pi/2 ) =  pi - EoD( EoD >  pi/2 );
        EoD( EoD < -pi/2 ) = -pi - EoD( EoD < -pi/2 );
        EoA( EoA >  pi/2 ) =  pi - EoA( EoA >  pi/2 );
        EoA( EoA < -pi/2 ) = -pi - EoA( EoA < -pi/2 );
        
        h_builder.taus = taus;
        h_builder.pow = pow;
        h_builder.AoD = AoD;
        h_builder.AoA = AoA;
        h_builder.EoD = EoD;
        h_builder.EoA = EoA;
        h_builder.NumClusters = n_paths;
        h_builder.NumSubPaths = ones(1,n_paths);
        
    otherwise
        error('QuaDRiGa:builder:generate_subpaths:notSupported',...
            ['Subpath generation method "',h_builder.scenpar.SubpathMethod,'" is not supported.']);
        
end

end

