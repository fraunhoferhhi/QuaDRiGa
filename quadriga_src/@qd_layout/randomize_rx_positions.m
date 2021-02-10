function randomize_rx_positions( h_layout, max_dist, min_height, max_height, track_length, rx_ind )
%RANDOMIZE_RX_POSITIONS Generates random Rx positions and tracks. 
%
% Calling object:
%   Single object
%
% Description:
%   Places the users in the layout at random positions. Each user will be assigned a linear track
%   with random direction. The random height of the user terminal will be in between 'min_height'
%   and 'max_height'.
%
% Input:
%   max_dist
%   the maximum distance from the layout center in [m]. Default is 50 m.
%
%   min_height
%   the minimum user height in [m]. Default is 1.5 m.
%
%   max_height
%   the maximum user height in [m]. Default is 1.5 m.
%
%   track_length
%   the length of the linear track in [m]. Default is 1 m.
%
%   rx_ind
%   a vector containing the receiver indices for which the positions should be generated. Default:
%   All receivers
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

if numel( h_layout ) > 1 
   error('QuaDRiGa:qd_layout:randomize_rx_positions','randomize_rx_positions not definded for object arrays.');
else
    h_layout = h_layout(1,1); % workaround for octave
end

% Parse input variables
if ~exist( 'max_dist' , 'var' ) || isempty( max_dist )
    max_dist = 50;
elseif ~( all(size(max_dist) == [1 1]) &&...
        isnumeric(max_dist) &&...
        isreal(max_dist) &&...
        max_dist > 0 )
    error('??? "max_dist" must be a real scalar  > 0')
end

if ~exist( 'min_height' , 'var' ) || isempty( min_height )
    min_height = 1.5;
elseif ~( all(size(min_height) == [1 1]) &&...
        isnumeric(min_height) &&...
        isreal(min_height) )
    error('??? "min_height" must be a real scalar')
end

if ~exist( 'max_height' , 'var' ) || isempty( max_height )
    max_height = 1.5;
elseif ~( all(size(max_height) == [1 1]) &&...
        isnumeric(max_height) &&...
        isreal(max_height) )
    error('??? "max_height" must be a real scalar  > 0')
end

if ~exist( 'track_length' , 'var' ) || isempty( track_length )
    track_length = 1;
elseif ~( all(size(track_length) == [1 1]) &&...
        isnumeric(track_length) &&...
        isreal(track_length) )
    error('??? "track_length" must be a real scalar  > 0')
end

if ~exist( 'rx_ind' , 'var' ) || isempty( rx_ind )
    rx_ind = 1:1:h_layout.no_rx;
elseif islogical( rx_ind )
    rx_ind = find( rx_ind );
end

% Generate random positions and tracks
rx_track = h_layout.rx_track;
no_tx = h_layout.no_tx;
samples_per_meter = h_layout.simpar.samples_per_meter;
for i_rx = 1:numel( rx_ind )
    n = rx_ind( i_rx );

    a = (2*rand-1)*max_dist + 1j*(2*rand-1)*max_dist;
    while abs(a)>max_dist
        a = (2*rand-1)*max_dist + 1j*(2*rand-1)*max_dist;
    end
    b = rand * (max_height - min_height) + min_height;
    
    if track_length > 0
        trk = qd_track.generate( 'linear',track_length );
    else
        trk = qd_track([]);
        trk.orientation(3,1) = rand*2*pi-pi;
    end
    trk.name = ['Rx',sprintf('%04d',n)];
    trk.initial_position = [ real(a) ; imag(a) ; b ];
    if track_length>0
        trk.interpolate_positions( samples_per_meter );
        trk.calc_orientation;
    end
    scenarios = rx_track(1,n).scenario(1);
    trk.scenario = scenarios( ones(1,no_tx),1 );
    
    % Write new track to layout
    rx_track(1,n) = trk;
end
h_layout.rx_track = rx_track;

end