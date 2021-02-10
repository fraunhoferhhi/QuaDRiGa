function h_channel_quant = quantize_delays( h_channel, tap_spacing, verbose )
%QUANTIZE_DELAYS Fixes the path delays to a grid of delay bins 
%
% Calling object:
%   Single object
%
% Description:
%   For some applications, e.g. channel emulation, it is not possible to achieve an infinite delay
%   accuracy. However, when the delays are rounded to a fixed grid of delay-bins (also refereed to
%   as "taps"), the time-evolving channel is no longer smooth. When a delay "jumps" from one delay-
%   bin to the next, e.g. when a MT is moving away from the BS, the phases in the frequency domain
%   representation of the channel will suddenly change as well. Multi-carrier communications
%   systems with closed-loop channel adaption (e.g. OFDM, WiFi, LTE, etc.) will exhibit poor
%   performance in this case. This method corrects this problem by approximating the "real" delay
%   value by two delays at a fixed spacing.  For example: when we assume that the required tap
%   spacing is 5 ns (which corresponds to a 200 MHz sample-rate) and the distance between BS and MT
%   is 10 m, the delay of the LOS path would be 33.4 ns. However, the fixed tap spacing only allows
%   values of 30 or 35 ns. This method approximates the LOS delay by two taps (one at 30 and one at
%   35 ns) and linear interpolation of the path power. Hence, in the frequency domain, the
%   transition from one tap to the next is smooth. But note: this only works when the bandwidth of
%   the communication system is significantly less than the sample-rate.
%
% Input:
%   tap_spacing
%   The spacing of the delay-bin in [seconds]. Default: 5 ns
%
%   verbose
%   Enables (1, default) or disables (0) a progress bar.
%
% Output:
%   h_channel_quant
%   A qd_channel object containing the approximated delays and channel coefficients
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

% Test if we have a scalar channel object
if numel( h_channel ) > 1
    error('QuaDRiGa:qd_channel:quantize_delays','"propsim_export" is only defined for scalar objects.')
else
    h_channel = h_channel(1,1); % workaround for octave
end

if ~exist( 'tap_spacing' , 'var' ) || isempty( tap_spacing )
    tap_spacing = 5e-9;
end

if exist( 'verbose' , 'var' )
    if ~( isnumeric(verbose) || ~isnan(verbose) || all(size(verbose) == [1 1]) )
        error('QuaDRiGa:qd_channel:quantize_delays','"verbose" must be numeric.')
    else
        verbose = logical(verbose);
    end
else
    verbose = true;
end

if verbose
    fprintf('Delay Quant. [');
    vb_dots = 50;
    tStart = clock;
    m0=0;
end

% Read common variables
no_rxant = uint32( h_channel.no_rxant );
no_txant = uint32( h_channel.no_txant );
no_snap  = uint32( h_channel.no_snap );

% We need individual delays for each antenna
if ~h_channel.individual_delays
    h_channel = copy( h_channel );
    h_channel.individual_delays = true;
end

% Copy delays and cefficients
% Use single-precision processing to reduce memory and computation-time contraints
D = single( h_channel.delay );
C = single( h_channel.coeff );
tap_spacing = single( tap_spacing );

% Each delay is approximated by two taps, one below the target delay and one above.
Dl = floor( D./tap_spacing );                           % Delay for the lower tap
Do = (D - Dl.*tap_spacing)./tap_spacing;             	% The relative position of the reayl delay between the two

% The power weigthing exponenet must match the target bandwidth to maintain the desired channel
% power. A value of 0.5 corresponds to the full sampling bandwidth (e.. when 200 MHZ at 5 ns tap
% spacing).
C = cat( 3, ((1-Do).^0.5).*C, (Do.^0.5).*C );           % Coefficnets are weightd by delay offset. 

% Combine Delays
Dl = uint16( Dl );                                      % Convert to uint16
D = cat( 3, Dl, Dl+1 );
clear Dl Do                                             % Free Memory

% Reorder dimensions to reduce loops (faster processing)
no_coeff = no_rxant*no_txant*no_snap;
D = reshape( permute( D , [1,2,4,3] ), no_coeff, [] );
C = reshape( permute( C , [1,2,4,3] ), no_coeff, [] );

% Get the unique delays
dl = sort( D(:) );       
ii = diff( dl ) ~= 0;
dl = dl( ii ).';
clear ii                                                % Free memory

% Add coefficients that have the same delay
no_delay_quant = numel(dl);
coeff = complex( zeros( no_coeff, no_delay_quant,'single' ) );
for d = 1 : no_delay_quant
    % Update progress bar
    if verbose; m1=ceil(d/no_delay_quant*vb_dots); if m1>m0;
            for m2=1:m1-m0; fprintf('o'); end; m0=m1; end; end;
    
    di = D == dl(d);        % Matching delays
    dc = any(di,1);         % Column indices
    dr = any(di,2);         % Row indices
    
    coeff(dr,d) = sum( C(dr,dc) .* di(dr,dc) ,2 );
end
clear D C di                                            % Free memory

% Remove zeros
ii = coeff ~= 0;
si = sum( ii,2 );
delay = zeros( no_coeff, no_delay_quant,'uint16' );
for c = 1 : no_coeff
    ij = 1:si(c);
    coeff( c,ij ) = coeff(c,ii(c,:));
    delay( c,ij ) = dl(ii(c,:));
    coeff( c,si(c)+1:end ) = 0;
end
no_delay_quant = max(si);
ii = 1:no_delay_quant;
coeff = coeff( :,ii );
delay = delay( :,ii );

% Format output coefficients
coeff = reshape( coeff, no_rxant, no_txant, no_snap, no_delay_quant );
coeff = permute( coeff, [1,2,4,3] );
coeff = double( coeff );

% Format output delays
delay = reshape( delay, no_rxant, no_txant, no_snap, no_delay_quant );
delay = permute( delay, [1,2,4,3] );
delay = double( delay ) .* double( tap_spacing );

% Write data to output channel
h_channel_quant = qd_channel( coeff, delay );
h_channel_quant.name = h_channel.name;
h_channel_quant.center_frequency = h_channel.center_frequency;
h_channel_quant.par = h_channel.par;
h_channel_quant.tx_position = h_channel.tx_position;
h_channel_quant.rx_position = h_channel.rx_position;

if verbose
    fprintf('] %5.0f seconds\n',round( etime(clock, tStart) ));
end

end
