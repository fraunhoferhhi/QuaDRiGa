function freq_response = fr( h_channel, bandwidth, carriers, i_snapshot )
%FR Transforms the channel into frequency domain and returns the frequency response
%
% Calling object:
%   Single object
%
% Input:
%   bandwidth
%   The baseband bandwidth in [Hz]
%
%   carriers
%   The carrier positions. There are two options:
%      * Specify the total number of carriers. In this case, 'carriers' a scalar natural number >
%        0. The carriers are then equally spaced over the bandwidth.
%      * Specify the pilot positions. In this case, 'carriers' is a vector of carrier positions.
%        The carrier positions are given relative to the bandwidth where '0' is the begin of the
%        spectrum and '1' is the end. For example, if a 5 MHz channel should be sampled at 0, 2.5 and
%        5 MHz, then 'carriers' must be set to [0, 0.5, 1].
%
%   i_snapshot
%   The snapshot numbers for which the frequency response should be calculated. By default, i.e. if
%   'i_snapshot' is not given, all snapshots are processed.
%
% Output:
%   freq_response
%   The complex-valued channel coefficients for each carrier in frequency domain. The indices of
%   the 4-D tensor are: [ Rx-Antenna , Tx-Antenna , Carrier-Index , Snapshot ]
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


if numel( h_channel ) > 1
    error('QuaDRiGa:qd_channel:fr','??? "fr" is only defined for scalar objects.')
else
    h_channel = h_channel(1,1); % workaround for octave
end

check = true;
if nargin < 3
    error('??? You must specify the bandwidth and the number of carriers.')
elseif nargin < 4
    i_snapshot = 1:h_channel.no_snap;
    check = false;
end

if ~( size(bandwidth,1) == 1 && isnumeric(bandwidth) && all(size(bandwidth) ==...
        [1 1]) && min(bandwidth) > 0 )
    error('??? The bandwidth "bandwidth" must be scalar and > 0')
end

if isnumeric(carriers) && isreal(carriers)
    
    if all(size(carriers) == [1 1]) && mod(carriers,1)==0 && carriers>0
        pilot_grid = ( 0:carriers-1 )/carriers;
        
    elseif numel( size(carriers) ) == 2 && any(size(carriers)==1)
        if size(carriers,2) == 1
            pilot_grid = carriers.';
        else
            pilot_grid = carriers;
        end
        carriers = numel(pilot_grid);
        
    else
        error('??? Invalid input for "carriers".')
    end
else
    error('??? The no. of carriers must be numeric.')
end

if check
    if ~( any( size(i_snapshot)==1 ) && isnumeric(i_snapshot) &&...
            all( mod(i_snapshot,1)==0 ) ...
            && min(i_snapshot) > 0 && max(i_snapshot)<=h_channel.no_snap )
        error(['??? The snapshot range must be numeric,',...
            ' integer and can not exceed the numbers of snapshots']);
    end
end

% Get the dimension of the channel tensor
n_rx = h_channel.no_rxant;
n_tx = h_channel.no_txant;
n_i_snapshots = numel(i_snapshot);
n_taps = h_channel.no_path;

% Preallocate some memory and rearrange coefficients
freq_response = zeros(n_i_snapshots * n_rx * n_tx, carriers);
if h_channel.individual_delays
    m = reshape(permute(h_channel.delay(:, :, :, i_snapshot)*bandwidth,...
        [ 4 1 2 3 ]), n_i_snapshots*n_rx*n_tx, n_taps);
else
    m = repmat(h_channel.delay(:, i_snapshot)'*bandwidth, n_rx*n_tx, 1);
end
c = reshape(permute(h_channel.coeff(:, :, :, i_snapshot), [ 4 1 2 3 ]),...
    n_i_snapshots*n_rx*n_tx, n_taps);

% The arguments of the exponential function
v = -2 * pi * 1j * pilot_grid;

% The main calculation
o_carriers = ones(1, carriers);
for i_tap = 1:n_taps
    freq_response = freq_response + c(:, o_carriers*i_tap) .* exp(m(:, i_tap) * v);
end

% Reorder the output dimensions
freq_response = reshape(freq_response, n_i_snapshots, n_rx, n_tx, carriers);
freq_response = permute(freq_response, [ 2 3 4 1 ]);

end
