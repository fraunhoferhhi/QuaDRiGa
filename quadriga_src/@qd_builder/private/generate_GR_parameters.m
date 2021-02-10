function gr_epsilon_r = generate_GR_parameters( tx_pos, rx_pos, f_GHz, gr_sos )
%GENERATE_GR_PARAMETERS Generates the parameters for the Ground Reflection
%
% Input:
%   tx_pos          Fixed Tx (Single mobility): Tx-position (metric) [ 3 x 1 ]
%                   Dual mobility: Tx-position (metric) [ 3 x N ]
%   rx_pos          Rx-positions (metric) [ 3 x N ]
%   f_GHz           Vector of carrier frequencies [ 1 x F ]
%   gr_sos          Pre-initialized SOS generatorsfor the parameter generation [ 1 x 1 ]
%
% Derived from input:
%   N  = number of users (from rx positions)
%   F  = number of frequencies (from f_GHz)
%
% Output
%   gr_epsilon_r    Relative permeability for ground reflection [ F x N ]
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
if ~exist('rx_pos','var') || size( rx_pos,1) ~= 3
    error('QuaDRiGa:qd_builder:generate_GR_parameters','Rx position is not given or has invalid format.')
end

% Number of users
N = size( rx_pos, 2 );
oN = ones(1,N);

if ~exist('tx_pos','var') || size( tx_pos,1) ~= 3
    error('QuaDRiGa:qd_builder:generate_GR_parameters','Tx position is not given or has invalid format.')
end

if size( tx_pos,2 ) == 1
    tx_pos = tx_pos(:,oN);
    tx_pos_SOS = [];
elseif size( tx_pos,2 ) == N
    tx_pos_SOS = tx_pos;
else
    error('QuaDRiGa:qd_builder:generate_GR_parameters','Invalid number of Tx-positions.')
end

if any( rx_pos(3,:) <= 0 ) || any( tx_pos(3,:) <= 0 )
    warning('QuaDRiGa:qd_builder:generate_GR_parameters',...
        'Terminal heights should be larger than 0 when using ground reflection.');
end

if ~exist('rx_pos','var') || isempty( f_GHz )
    error('QuaDRiGa:qd_builder:generate_GR_parameters','Carrier frequencies are not given.')
end

F = numel( f_GHz );

if ~exist('gr_sos','var') || isempty( gr_sos )
    gr_sos = [];
    
elseif isnumeric( gr_sos ) && numel( gr_sos ) == 1
    
    if gr_sos == 0
        gr_sos = [];
    else % Create sos_model
        SC_lambda = gr_sos;
        acf = 'Comb300';
        gr_sos = qd_sos( acf, 'Uniform', SC_lambda );
        gr_sos.sos_phase(:,2) = gr_sos.sos_phase(:,1);
    end
    
elseif isa( gr_sos, 'qd_sos' ) && numel( gr_sos ) == 1 
    % Test if all generators are set to uniform distribution
    if ~strcmp( gr_sos.distribution , 'Uniform' )
        warning('QuaDRiGa:qd_builder:generate_GR_parameters',...
            'The ground reflection random generator should be set to uniform distribution.');
    end
   
    % Test if GR is reciprocal
    tmp = gr_sos.sos_phase;
    tmp = tmp(:,1,:) - tmp(:,2,:);
    if any(abs(tmp(:))>1e-6)
        warning('QuaDRiGa:qd_builder:generate_GR_parameters','Ground reflection is not reciprocal.');
    end
else
    error('QuaDRiGa:qd_builder:generate_paths','gr_sos has invalid format.');
end


% Parameter generation
% Generate spatially correlated random variables
if isempty( gr_sos )
    randC = rand( 1,N );
else
    randC = val( gr_sos, rx_pos,tx_pos_SOS );  % Uniform distribution
end

% There are 3 ground types defined: dry, medium dry, wet
% The complex-valued relative permittivity is frequency-dependent and given by:
g_dry = 3 + 1j * 0.003 * f_GHz.^1.34;
g_med = 30.4 * f_GHz.^-0.47 + 1j * 0.18 * f_GHz.^1.05;
g_wet = 31.3 * f_GHz.^-0.48 + 1j * 0.63 * f_GHz.^0.77;

% The reulting permittivity is obtained by a linear interpolation of the ground type
% depending on the random variable "randC". A value of 0 means "dry", 0.5 means "medium"
% and 1 means "wet".
gr_epsilon_r = zeros(F,N);

i1 = randC <= 0.5;        % Dry to Medium wet ground
i2 = randC >  0.5;        % Medium wet to wet ground

% The weights for the linear interpolation
w1 = 2 * randC(i1);
w2 = 2 * (randC(i2)-0.5);

% The relative permittivity of the ground at the different frequencies
if any(i1)
    gr_epsilon_r( :,i1 ) = g_dry.'*(1-w1) + g_med.'*w1;
end
if any(i2)
    gr_epsilon_r( :,i2 ) = g_med.'*(1-w2) + g_wet.'*w2;
end

end
