function [ V, H, dist] = interpolate( h_qd_arrayant, azimuth, elevation, i_element,...
    azimuth_grid, elevation_grid, Fa, Fb, element_position )
%INTERPOLATE Interpolates the field pattern
%
% Calling object:
%   Single object
%
% Description:
%   Interpolation of the beam patterns is very computing intensive. It must be performed several
%   thousands of times during a simulation run. Therefore, optimized 2D linear interpolation is
%   used. There are additional input parameters specified in the .m-File that are not in the list
%   below. Those parameters correspond to the properties of the 'qd_arrayant' class. Passing those
%   variables during the function call takes less time than reading them from the object
%   properties. This is used internally in 'qd_builder.get_channels' but is irrelevant here.
%
% Input:
%   azimuth
%   A vector of azimuth angles in [rad]
%
%   elevation
%   A vector of elevation angles in [rad]
%
%   i_element
%   The element indices for which the interpolation is done. If no element index is given, the
%   interpolation is done for all elements in the array.
%
% Output:
%   V
%   The interpolated vertical field pattern (e-theta-component)
%
%   H
%   The interpolated horizontal field pattern (e-phi-component)
%
%   dist
%   The effective distances between the antenna elements when seen from the direction of the
%   incident path. The distance is calculated by an projection of the array positions on the normal
%   plane of the incident path. This is needed for the planar wave approximation.
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

if numel( h_qd_arrayant ) > 1
    error('QuaDRiGa:qd_arrayant:interpolate','interpolate not definded for object arrays.');
else
    h_qd_arrayant = h_qd_arrayant(1,1); % workaround for octave
end

% Determine the floating point precision
if isa( azimuth , 'single' )
    use_single_precision = true;
    precision = 'single';
else
    use_single_precision = false;
    precision = 'double';
end

% Parse input arguments
if nargin < 5
    if nargin < 3
        error('??? You need to specify azimuth and elevation for the pattern interpolation.');
    elseif nargin < 4
        i_element = 1:h_qd_arrayant.no_elements;
    end
    
    if use_single_precision
        [ Fa, Fb, azimuth_grid, elevation_grid, element_position ] = h_qd_arrayant.wrap_grid( [], 'single' );
    else
        [ Fa, Fb, azimuth_grid, elevation_grid, element_position ] = h_qd_arrayant.wrap_grid;
    end
    Fa = reshape( Fa, [], h_qd_arrayant.no_elements );
    Fb = reshape( Fb, [], h_qd_arrayant.no_elements );
end

% Get initial values
dims = size( azimuth );
no_values = numel( azimuth );
no_element = numel( i_element );
no_az = numel( azimuth_grid );
no_el = numel( elevation_grid );

% Reduce dimensions of arrays
azimuth = reshape( azimuth,1,no_values );
elevation = reshape( elevation,1,no_values );

% elevation may be outside of (-pi/2,pi/2).
% Mapping of angles to (-pi/2,pi/2) and changing the azimuthal orientation of rays
% with original elevation angles within +-(pi/2,3pi/2)+n*2*pi. 29.9.2008 EsaK

ind = elevation > pi | elevation <= -pi;
elevation(ind) = mod( elevation(ind)+pi , 2*pi )-pi;

ind = elevation > pi/2;
elevation(ind) = pi - elevation(ind);
azimuth(ind) = azimuth(ind) + pi;

ind = elevation < -pi/2;
elevation(ind) = -pi - elevation(ind);
azimuth(ind) = azimuth(ind) + pi;

ind = azimuth > pi | azimuth <= -pi;
azimuth(ind) = mod( azimuth(ind)+pi , 2*pi )-pi;

% Interpolate field patterns
tmp_1_no_vals = uint32( (1:no_values) );

% Determine the nearest location of xi in x and the difference to
% the next point
[tmp,b] = sort( azimuth );
[~,a]   = sort( [azimuth_grid,tmp] );
ui      = uint32( 1:(no_az + no_values) );
ui(a)   = ui;
ui      = ui(no_az+1:end) - tmp_1_no_vals;
ui(b)   = ui;
ui( ui==no_az ) = no_az-1;
ui( ui==0 ) = 1;
if numel( azimuth_grid ) > 1
    uin = ui+1;
    u = (azimuth-azimuth_grid(ui))./( azimuth_grid(uin)-azimuth_grid(ui) );
    u = u';
else
    uin = 1;
    u = 0;
end

% Determine the nearest location of yi in y and the difference to
% the next point
[tmp,b] = sort( elevation );
[~,a]   = sort( [elevation_grid,tmp] );
vi      = uint32( 1:(no_el + no_values) );
vi(a)   = vi;
vi      = vi(no_el+1:end) - tmp_1_no_vals;
vi(b)   = vi;
vi( vi==no_el ) = no_el-1;
vi( vi==0 ) = 1;
if numel( elevation_grid ) > 1
    vin = vi+1;
    v = (elevation-elevation_grid(vi))./( elevation_grid(vin)-elevation_grid(vi) );
    v = v';
else
    vin = 1;
    v = 0;
end

% Illustration of the Interpolation procedure
% 
%      c -----f------------------- d
%      |      |                    |
%      |      |(1-v)               |
%      |      |                    |
%      |      |                    |
%      |  u   |         (1-u)      |
%      g----- x ------------------ h
%      |      |                    |
%      |      |v                   |
%      |      |                    |
%      a------e--------------------b

% Determine the indices of the elements
pa = vi  + ( ui  -1 )*no_el;
pb = vi  + ( uin -1 )*no_el;
pc = vin + ( ui  -1 )*no_el;
pd = vin + ( uin -1 )*no_el;

V = phase_interpolation_1d( Fa(pa,i_element), Fa(pb,i_element), u(:,ones(1,no_element)) ); % Point e
if numel( elevation_grid ) > 1
    f = phase_interpolation_1d( Fa(pc,i_element), Fa(pd,i_element), u(:,ones(1,no_element)) );
    V = phase_interpolation_1d( V, f, v(:,ones(1,no_element)) );
end

H = phase_interpolation_1d( Fb(pa,i_element), Fb(pb,i_element), u(:,ones(1,no_element)) ); % Point e
if numel( elevation_grid ) > 1
    f = phase_interpolation_1d( Fb(pc,i_element), Fb(pd,i_element), u(:,ones(1,no_element)) );
    H = phase_interpolation_1d( H, f, v(:,ones(1,no_element)) );
end

% Remap output to match input dimensions
V = reshape( V, [dims,no_element] );
H = reshape( H, [dims,no_element] );

% The effective distances are only needed when we calculate the channel
% coefficients. Therefore, it is also provided by the interpolation
% function.

if nargout > 2
    % Compute the effective distances for each path
    
    % Transform incoming angles from spherical coordinates to Cartesian
    % coordinates.
    if any(any( element_position~=0 ))
        A = zeros( 3, no_values, precision );
        [A(1, :), A(2, :), A(3, :)] = sph2cart( azimuth, elevation, 1 );
    end
    
    % The distances are calculated by a parallel projection.
    % See: http://de.wikipedia.org/wiki/Parallelprojektion
    % The normal vector of the projection plane is given in X.
    % The original point P is the element position defined by the array antenna
    % geometry. The origin of the projection plane is the origin of the
    % array antenna coordinate system. The projection result is B.
    
    dist = zeros( no_values, no_element, precision );
    
    for n = 1 : no_element
        E = element_position(:,n);
        
        if any( E ~= 0 )            % The distance is 0 if P is in the origin
            tmp = E'*A;             % cos(E,A)
            D = tmp([1 1 1],:) .* A;
            dist(:,n) = -sign( A(1,:).*D(1,:) + A(2,:).*D(2,:) + A(3,:).*D(3,:) ) .* ...
                sqrt(sum( D.^2 ));
            
            % This is the same as
            % for m = 1:no_values
            %    D = -( -E.' * A(:,m) ) * A(:,m) ;
            %    dist(m,n) = sign(D'*A(:,m)) * sqrt(sum(  D.^2 ));
            % end
        end
    end
    
    % Return the output
    dist = reshape( dist, [dims,no_element] );
end

end
