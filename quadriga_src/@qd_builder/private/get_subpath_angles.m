function [ aod,eod,aoa,eoa,delay ] = get_subpath_angles( h_builder,i_mobile )
%GET_SUBPATH_ANGLES Generate subpaths and perform random coupling according to the 3GPP baseline model
%
%   GET_SUBPATH_ANGLES generates the subpaths around the each path and
%   randomly couples the subpaths on the Tx- and Rx side. 
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

n_clusters  = h_builder.NumClusters;                     % no. clusters
n_subpath   = h_builder.NumSubPaths;
n_paths     = sum( h_builder.NumSubPaths );              % no. paths

% The origianl 20 offset angles
% offset = [0.0447 0.1413 0.2492 0.3715 0.5129 0.6797 0.8844 1.1481 1.5195 2.1551];
% offset = [ offset , -offset ]*0.017453292519943;  % in rad

% The optimized order of the 20 sub-paths for 
offset = [-0.8844 1.1481 0.6797 -1.1481 -0.6797 1.5195 -0.0447 0.3715 -1.5195,...
    -0.2492 -0.5129 0.5129 0.8844 -2.1551 -0.1413 0.1413 0.0447 2.1551 -0.3715 0.2492 ];
%offset = offset * 0.017453292519943;  % in rad

% The per cluster angular spread scaling coefficients
c_aod = h_builder.scenpar.PerClusterAS_D;    % in deg
c_aoa = h_builder.scenpar.PerClusterAS_A;    % in deg
c_eod = h_builder.scenpar.PerClusterES_D;    % in deg
c_eoa = h_builder.scenpar.PerClusterES_A;    % in deg

% Reserve some memory for the output
aod = zeros( 1,n_paths );    
aoa = aod;
eod = aod;
eoa = aod;

if size( h_builder.subpath_coupling, 3 ) == 1
    subpath_coupling = h_builder.subpath_coupling(:,:,1);
else
    subpath_coupling = h_builder.subpath_coupling(:,:,i_mobile);
end

ls = 1;
for l = 1 : n_clusters
    le = ls + n_subpath(l) - 1;
    
    % Get the offset angles
    if n_subpath(l) == 1
        of = 0;
    else
        of = offset( 1:n_subpath(l) );
        if n_subpath(l) < 20
            of = (of-mean(of));
            of = of./sqrt(mean(of.^2));
        end
        of = of .* 0.017453292519943;  % in rad
    end
    
    % Change order to account for random sub-path coupling
    [~,ind1] = sort( subpath_coupling(1,ls:le) );
    [~,ind2] = sort( subpath_coupling(2,ls:le) );
    [~,ind3] = sort( subpath_coupling(3,ls:le) );
    [~,ind4] = sort( subpath_coupling(4,ls:le) );
    aod(ls:le) = c_aod * of( ind1 ) + h_builder.AoD( i_mobile,l );
    aoa(ls:le) = c_aoa * of( ind2 ) + h_builder.AoA( i_mobile,l );
    eod(ls:le) = c_eod * of( ind3 ) + h_builder.EoD( i_mobile,l );
    eoa(ls:le) = c_eoa * of( ind4 ) + h_builder.EoA( i_mobile,l );
    
    ls = le + 1;
end

% Calculate delays
if nargout == 5
    dl = h_builder.taus(i_mobile,:);
    n_snapshots = h_builder.rx_track(1,i_mobile).no_snapshots;
    if h_builder.simpar.use_absolute_delays
        r_0 = h_builder.rx_positions(:,i_mobile) - h_builder.tx_position(:,i_mobile);
        D = sqrt(sum(r_0.^2)) / h_builder.simpar.speed_of_light;
        delay = repmat( dl+D , n_snapshots , 1  );
    else
        delay = repmat( dl , n_snapshots , 1  );
    end
end

end
