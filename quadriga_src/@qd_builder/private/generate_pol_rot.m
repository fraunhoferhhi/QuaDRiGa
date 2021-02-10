function [ gamma, kappa ] = generate_pol_rot( Ln, M, tx_pos, rx_pos, xpr_mu, xpr_sigma, xpr_sos )
%GENERATE_XPR Generates the NLOS polarization rotation for the subpaths 
%
% Input:
%   Ln              Number of NLOS clusters (excluding LOS and GR) [ 1 x 1 ]
%   M               Number of paths per clusters [ 1 x 1 ]
%   tx_pos          Fixed Tx (Single mobility): Tx-position (metric) [ 3 x 1 ]
%                   Dual mobility: Tx-position (metric) [ 3 x N ]
%   rx_pos          Rx-positions (metric) [ 3 x N ]
%   xpr_mu          Per-user XPR from the LSF model in (dB) [ F x N ]
%   xpr_sigma       XPR STD in (dB) [ F x N ]
%   xpr_sos         Pre-initialized SOS generators [ Ln x M+1 ]
%                   - Linear polarization is generated independently for each sub-path (M)
%                   - Elliptic polarization is generated per cluster (+1)
%
% Derived from input:
%   Ln = number of NLOS clusters (from onput)
%   M  = number of paths per clusters (from input)
%   N  = number of users (from rx positions)
%   F  = number of frequencies (xpr_mu)
%
% Output
%   gamma           Linear polarization rotation angle in (rad) [ N x Ln*M x F ]
%   kappa           Elliptic polarization rotation angle in (rad) [ N x Ln*M x F ]
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
if ~exist('Ln','var') || isempty( Ln ) || numel( Ln ) ~=1
    error('QuaDRiGa:qd_builder:generate_pol_rot','Numer of NLOS clusters must be given.')
end

if ~exist('M','var') || isempty( M ) || numel( M ) ~=1
    error('QuaDRiGa:qd_builder:generate_pol_rot','Numer of paths per cluster must be given.')
end
oM = ones(1,M);

if ~exist('rx_pos','var') || size( rx_pos,1) ~= 3
    error('QuaDRiGa:qd_builder:generate_pol_rot','Rx position is not given or has invalid format.')
end

% Number of users
N = size( rx_pos, 2 );

if ~exist('tx_pos','var') || size( tx_pos,1) ~= 3
    error('QuaDRiGa:qd_builder:generate_pol_rot','Tx position is not given or has invalid format.')
end

if size( tx_pos,2 ) == 1
    tx_pos_SOS = [];
elseif size( tx_pos,2 ) == N
    tx_pos_SOS = tx_pos;
else
    error('QuaDRiGa:qd_builder:generate_pol_rot','Invalid number of Tx-positions.')
end

if ~exist('xpr_mu','var') || size( xpr_mu,2) ~= N 
    error('QuaDRiGa:qd_builder:generate_pol_rot','XPR mu is not given or has invalid format.')
end

F = size( xpr_mu,1 );

if ~exist('xpr_sigma','var') || size( xpr_sigma,2) ~= N 
    error('QuaDRiGa:qd_builder:generate_pol_rot','XPR sigma is not given or has invalid format.')
end

if ~exist('xpr_sos','var') || isempty( xpr_sos )
    xpr_sos = [];
    
elseif isnumeric( xpr_sos ) && numel( xpr_sos ) == 1
    
    if xpr_sos == 0
        xpr_sos = [];
    else % Create sos_model
        SC_lambda = xpr_sos;
        acf = 'Comb300';
        xpr_sos = qd_sos([]);
        for i_freq = 1 : F
            for i_cluster = 1 : Ln
                for i_sub = 1 : M+1
                    xpr_sos(i_cluster,i_sub,i_freq) = qd_sos( acf, 'Normal', SC_lambda );
                    xpr_sos(i_cluster,i_sub,i_freq).sos_phase(:,2) = xpr_sos(i_cluster,i_sub).sos_phase(:,1);
                end
            end
        end
    end
    
elseif isa( xpr_sos, 'qd_sos' ) && size(xpr_sos,1) == Ln && size(xpr_sos,2) == M+1 && ( F == 1 || size(xpr_sos,3) == F )
    % Test if all generators are set to uniform distribution
    for n = 1 : size( xpr_sos, 1)
        for m = 1 : size( xpr_sos, 2)
            if ~strcmp( xpr_sos(n,m).distribution , 'Normal' )
                warning('QuaDRiGa:qd_builder:generate_pol_rot:Normal',...
                    'All random generators should be set to uniform distribution.');
                warning ('off','QuaDRiGa:qd_builder:generate_pol_rot:Normal');
            end
            tmp = xpr_sos(n,m).sos_phase;
            tmp = tmp(:,1,:) - tmp(:,2,:);
            if any(abs(tmp(:))>1e-6)
                warning('QuaDRiGa:qd_builder:generate_pol_rot:reciprocal','XPR SOS generators are not reciprocal.');
                warning ('off','QuaDRiGa:qd_builder:generate_pol_rot:reciprocal');
            end
        end
    end
    warning ('on','QuaDRiGa:qd_builder:generate_pol_rot:Normal');
    warning ('on','QuaDRiGa:qd_builder:generate_pol_rot:reciprocal');
else
    error('QuaDRiGa:qd_builder:generate_pol_rot','xpr_sos has invalid format.');
end

% Value generation
oLnM = ones(1,Ln*M);

gamma = zeros( N, Ln*M, F );
kappa = zeros( N, Ln*M, F );

for iF = 1 : F       % Polarization rotation is indepentently genertated for different frequencies
    
    % Generate spatially correlated random variables for the linear polarization offset
    if isempty( xpr_sos )
        randC = randn( N,Ln*M );
    else
        randC = val( xpr_sos(:,1:M,iF), rx_pos, tx_pos_SOS ).';   % Normal distribution
        randC = reshape( randC, N, Ln, M );
        randC = reshape( permute( randC, [1,3,2] ) , N,Ln*M );
    end
    
    xpr_linear = randC .* ( xpr_sigma(iF,:)' * oLnM ) + xpr_mu(iF,:).' * oLnM; % dB
    gamma(:,:,iF) = acot( sqrt( 10.^(0.1*xpr_linear) ) );
    
    % Generate spatially correlated random variables for the circular polarization offset
    if isempty( xpr_sos )
        randC = randn( N,Ln*M );
    else
        randC = val( xpr_sos(:,M+1,iF), rx_pos, tx_pos_SOS ).';   % Normal distribution
        randC = randC( :,:,oM );
        randC = reshape( permute( randC, [1,3,2] ) , N,Ln*M );
    end
    
    xpr_circular = randC .* ( xpr_sigma(iF,:)' * oLnM ) + xpr_mu(iF,:).' * oLnM; % dB
    kappa(:,:,iF) = acot( sqrt( 10.^(0.1*xpr_circular) ) );
end

end
