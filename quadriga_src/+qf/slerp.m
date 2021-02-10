function [ phiI, thetaI, pI ] = slerp( x, phi, theta, xi )
%SLERP Spherical linear interpolation optimized for performace
%
% Description:
%   Slerp is shorthand for spherical linear interpolation and refers to constant-speed motion along
%   a unit-radius great circle arc. This function is needed when interpolating angles in spheric or
%   circular coordinates, e.g. when interpolating phase information (qd_channel.interpolate) or
%   orientations (qd_track.interpolate).  
%
%   Circular linear interpolation can be done by using
%   phiI = slerp( x, phi, 0, xi )
%
% Input:
%   x
%   Vector of x sample points; size [ 1, nx ] or [ nx, 1 ]
%
%   phi
%   Vector of phi angles in [rad]; values from -pi to pi; size [ nx, ne ]; the 2nd dimension allows
%   for interpolations of multiple data-sets
%
%   theta
%   Vector of theta angles in [rad]; values from -pi/2 (down) to pi/2 (up); size [ nx, ne ]; the
%   2nd dimension allows for interpolations of multiple data-sets
%
%   xc
%   Vector of x sample points after interpolation; size [ 1, nxi ] or [ nxi , 1 ]
%
% Output:
%   phiI
%   Vector of interpolated phi angles in [rad]; values from -pi to pi; size [ nxi, ne ]
%
%   thetaI
%   Vector of interpolated theta angles in [rad]; values from -pi/2 (down) to pi/2 (up); 
%   size [ nxi, ne ]
%
%   pI
%   Cartesian coordinates of the interpolated points on the unit-sphere; size [ 3, nxi, ne ]
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

if isa(phi,'double')
    phi = single( phi );
    phi_is_double = true;
else
    phi_is_double = false;
end

x  = single( x(:).' );
xi = single( xi(:).' );

nx  = numel( x );
ni  = numel( xi );

if size( phi,1 ) == 1 && size( phi,2 ) == nx
    % OK
elseif size( phi,1 ) ~= nx || numel( size(phi) ) > 2
    error('QuaDRiGa:qf:slerp','Size of "phi" does not match size of "x".');
else
    phi = permute( phi, [3,1,2] );
end

nph = size( phi,2 );
ne = size( phi,3 );
oe = ones(1,ne,'uint8');

nth = numel( theta );
if nth == 1 % Option for circular interpolation
    theta = zeros(1,nx,ne,'single');
    sinth = theta;
    costh = ones(1,nx,ne,'single');
elseif size( theta,1 ) == 1 && size( theta,2 ) == nx
    theta = single( theta );
    costh = cos(theta);
    sinth = sin(theta);
elseif size( theta, 1 ) ~= nx || size( theta, 2 ) ~= ne || numel( size(theta) ) > 2
    error('QuaDRiGa:qf:slerp','Size of "theta" does not match size of "phi".');
else
    theta = single( permute( theta, [3,1,2] ) );
    costh = cos(theta);
    sinth = sin(theta);
end

% Transform orientations to Cartesian coordinates
p = zeros( 3,nx,ne,'single' );
p(1,:,:) = costh .* cos(phi);
p(2,:,:) = costh .* sin(phi);
p(3,:,:) = sinth;

% Calculate the angles "omega" between two consecutive points in 3D space
tmp         = sum( p(:,1:end-1,:) .* p(:,2:end,:) );
tmp(tmp>1)  = 1;
tmp(tmp<-1) = -1;
omega       = acos( tmp );
sinomega    = sin( omega );

% Determine the nearest location of xi in x and the difference to the next point
ii      = uint32( 1:ni );
[tmp,b] = sort( xi );
[~,a]   = sort( [x,tmp] );
ui      = uint32( 1:(nx + ni) );
ui(a)   = ui;
ui      = ui(nx+1:end) - ii;
ui(b)   = ui;
ui( ui==nx ) = nx-1;
ui( ui==0 ) = 1;
uin     = ui+1;
u       = (xi-x(ui))./( x(uin)-x(ui) );

% Calculate the slerp scaling coefficients
% Form small offset angles, the slerp interpolations becomes numerically insteble, because
% sin(omega) = 0 (division by zero). Hence, we approximate this by a linear interpolation, assuming
% that sin(u) ~ u  for small u.
iL = omega < 1e-3;
omega( iL ) = 1e-3;
sinomega( iL ) = 1e-3;

uN = 1-u;
c1 = sin( uN(1,:,oe) .* omega( 1,ui,: ) ) ./ sinomega( 1,ui,: );
c2 = sin(  u(1,:,oe) .* omega( 1,ui,: ) ) ./ sinomega( 1,ui,: );

% Interpolate points
o3 = ones(1,3,'uint8');
pI = c1(o3,:,:) .* p(:,ui,:) + c2(o3,:,:) .* p(:,uin,:);

% Transform orientations to spherical coordinates
phiI = permute( atan2(pI(2,:,:),pI(1,:,:)) , [2,3,1] );
if nargout > 1
    if nth > 1
        hypotxy = hypot(pI(1,:,:),pI(2,:,:));
        thetaI  = permute( atan2(pI(3,:,:),hypotxy) , [2,3,1] );
        if phi_is_double
            thetaI = double( thetaI );
        end
    elseif phi_is_double
        thetaI  = zeros( ni,ne );
    else
        thetaI  = zeros( ni,ne,'single' );
    end
end

% Typecast outputs to double if needed
if phi_is_double
    phiI = double( phiI );
    if nargout > 2
        pI = double( pI );
    end
end

end
