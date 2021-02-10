function h_qd_arrayant = gen_arrayant_parabolic( r, f, p_min, pol )
%GEN_ARRAYANT_PARABOLIC
%
%   An ideal parabolic reflector antenna with input parameters:
%
%      * r - Radius of the antenna aperture in [meters]
%      * f - Center frequency in [Hz]
%      * p_min - Min. sidelobe power relative to directivity in [dB] (default: -40 dB)
%      * pol - Polarization indicator
%           1. vertical (E-theta) polarization (default)
%           2. horizontal (E-phi) polarization
%           3. LHCP
%           4. RHCP
%           5. dual-polarized two-port antenna (LHCP,RHCP)
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

% Set inputs
if ~exist('r','var') || isempty( r )
    r = 1;
end

if ~exist('f','var') || isempty( f )
    f = 4e9;
end

if ~exist('p_min','var') || isempty( p_min )
    p_min = -40;
end

if ~exist('pol','var') || isempty( pol )
    pol = 1;
end

% Wave number
k = 2 * r * pi * f / qd_simulation_parameters.speed_of_light;

% Calc pattern vs az. angle at low resolution
azimuth_grid = 0:0.001:90;                  % AZ-Scan
z = k .* sin(azimuth_grid*pi/180);
e = 2 * besselj(1, z ) ./ z;                % Electric field
e(1) = 1;
p = 10*log10(abs(e).^2);                    % Power

dp = abs(diff(p));

ind_e = find( dp > 5, 1,'last' ) + 1;       % The last Null
ind_m = find( p > p_min , 1, 'last' );      % The last time wherer the pattern is > p_min
if ind_m < ind_e
    dp(1:ind_m) = 0;
    ind_m = find( dp > 5, 1,'first' ) + 1;  % The first Null after -p_min
else
    ind_m = ind_e;
end
az_max = azimuth_grid(ind_m);

% Use the SAGE algorithm to find the true az_max
a   = az_max;       % Initial value
x   = p(ind_m);
dm  = 0.001;        % Step size
delta = Inf; ddir = +1; lp = 2;
while lp<100 && delta > 1e-9
    if lp>1; an = a + ddir * dm; delta = abs(a-an); else an = a; end
    
    z = k .* sin(an*pi/180);
    e = 2 * besselj(1, z ) ./ z;    % Electric field
    xn = 10*log10(abs(e).^2);       % Power
    
    if xn < x; a = an; x = xn; else ddir = -ddir; dm = 0.2 * dm; end
    lp = lp + 1;
end
az_max = an;

% Calculate 2D pattern
no_steps = round( az_max*k/10 );            % Number of steps
step_size = az_max/no_steps;                % Step size

% Determine the reduced grid
gridA = (-az_max : step_size : az_max)*pi/180;

% Calculate pattern in 2D
cosG = cos(gridA);
omega = acos( cosG' * cosG );               % Angles between the x-axis and the point on the sphere

z = k .* sin(omega);
e = 2 * besselj(1, z ) ./ z;                % Electric field

ii = omega < 0.5*step_size*pi/180;          % Value at the x-axis
e(ii) = 1;

e( omega>az_max*pi/180 ) = 0;               % Cutoff

% Assemble Quadriga array antenna object
h_qd_arrayant = qd_arrayant([]);
h_qd_arrayant.name = 'parabolic';
h_qd_arrayant.center_frequency = f;
h_qd_arrayant.set_grid( gridA, gridA, 0 );
h_qd_arrayant.Fa = e;
h_qd_arrayant.normalize_gain;

switch pol
    case 2
        h_qd_arrayant.rotate_pattern(90,'x',1,2);
    case 3
        h_qd_arrayant.copy_element(1,2);
        h_qd_arrayant.rotate_pattern(90,'x',2,2);
        h_qd_arrayant.coupling = 1/sqrt(2) * [1;1j];
        h_qd_arrayant.combine_pattern;
    case 4
        h_qd_arrayant.copy_element(1,2);
        h_qd_arrayant.rotate_pattern(90,'x',2,2);
        h_qd_arrayant.coupling = 1/sqrt(2) * [1 ; -1j];
        h_qd_arrayant.combine_pattern;
    case 5
        h_qd_arrayant.copy_element(1,2);
        h_qd_arrayant.rotate_pattern(90,'x',2,2);
        h_qd_arrayant.coupling = 1/sqrt(2) * [1 1 ; 1j -1j];
        h_qd_arrayant.combine_pattern;
end

end
