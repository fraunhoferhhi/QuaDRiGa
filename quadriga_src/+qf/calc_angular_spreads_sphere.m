function [ as, es, orientation, phi, theta ] = calc_angular_spreads_sphere( az, el, pow, calc_roll_angle )
%CALC_ANGULAR_SPREADS Calculates the angular spread in [rad]
%
% Description:
%   This function calculates the RMS angular spread from a given set of angles. It is used by the
%   "qd_builder" to map the path powers to angles.
%
% Input:
%   ang
%   A vector of angles in [rad]. Dimensions: [ n_ang x n_path ]
%
%   pow
%   A vector of path powers in [W]. Dimensions: [ n_ang x n_path ] or [ 1 x n_path ]
%
%   wrap_angles
%   A logical variable. If set to 1, angles will be wrapped around +/- pi. If set to 0, no wrapping
%   is applied. Default: 1 (with wrapping)
%
% Output:
%   as
%   The RMS angular spread in [rad] for each angle vector. Dimensions: [ n_ang x 1 ]
%
%   mean_angle
%   The mean angles in [rad]. Dimensions: [ n_ang x 1 ]
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

if ~exist('calc_roll_angle','var') || isempty( calc_roll_angle );
    calc_roll_angle = true;
end

N = size( az,1 );
L = size( az,2 );
oL = ones( 1,size(az,2 ));

if size( el,1 ) < N
    el = zeros( size( az ));
end
if size( pow,1 ) < N
    pow = pow( ones(1,N), : );
end

% Normalize powers
pt = sum( pow,2 );
pow = pow./pt( :,oL  );

% Pre-process data
powX  = permute( pow, [3,2,1] );
phi   = permute( az,  [3,2,1] );
theta = permute( el,  [3,2,1] );

% Placeholder
C     = zeros( 3,L,N );
R     = zeros( 3,3,N );
yaw   = angle( sum( powX.*exp( 1j*phi ) , 2 ) );              	% Calculate yaw angle
pitch = zeros( 1,1,N );
roll  = zeros( 1,1,N );
orientation = zeros( 3,1,N );

upd = true(1,N);
lp  = 1;
while any( upd ) && lp < 10
    
    phi(:,:,upd) = phi(:,:,upd) - yaw( :,oL,upd );              % Apply yaw angle (z-rotation)
    orientation(3,1,upd) = orientation(3,1,upd) + yaw(1,1,upd); % Save orientation change
    
    pitch(:,:,upd) = angle( sum( powX(:,:,upd).*exp( 1j*theta(:,:,upd) ) , 2 ) );  % Calculate pitch angle
    
    R(:) = 0;                                                   % Calculate rotation matrix
    co = cos( pitch(:,:,upd) );
    si = sin( pitch(:,:,upd) );
    R(1,1,upd) = co;
    R(1,3,upd) = si;
    R(3,1,upd) = -si;
    R(3,3,upd) = co;
    R(2,2,upd) = 1;
    
    C(1,:,upd) = cos( phi(:,:,upd) ) .* cos( theta(:,:,upd) ); % Transform to Cartesian coordinates
    C(2,:,upd) = sin( phi(:,:,upd) ) .* cos( theta(:,:,upd) );
    C(3,:,upd) = sin( theta(:,:,upd) );
    
    for n = 1 : N                                               % Apply pitch angle
        if upd( n )
            C(:,:,n) = R(:,:,n) * C(:,:,n);
        end
    end
    
    theta(:,:,upd) = atan2(C(3,:,upd),hypot( C(1,:,upd),C(2,:,upd) ));  % Transform to Spheric coordinates
    phi(:,:,upd)   = atan2(C(2,:,upd),C(1,:,upd));
    orientation(2,1,upd) = orientation(2,1,upd) + pitch(1,1,upd);       % Save orientation change
    
    yaw(:,:,upd) = angle( sum( powX(:,:,upd).*exp( 1j*phi(:,:,upd) ) , 2 ) );      % Calculate yaw angle
    pitch(:,:,upd) = angle( sum( powX(:,:,upd).*exp( 1j*theta(:,:,upd) ) , 2 ) );  % Calculate pitch angle
    
    upd = abs(yaw(1,:)) > 0.001 | abs(pitch(1,:)) > 0.001;
    lp = lp + 1;
end

% Calculate angular spreads
as = qf.calc_angular_spreads( permute(phi,[3,2,1]), pow );
es = qf.calc_angular_spreads( permute(theta,[3,2,1]), pow );

% Try to find the roll angle that maximized the AS and minimized the ES
if calc_roll_angle
    
    cst   = es./as;                                 % Initial cost function at 0 deg roll angle
    cstN  = cst;                                    % Placeholder for cost function updated
    rollN = roll;                                   % Placeholder for roll angle updates
    thetaN = theta;
    phiN  = phi;
    dir   = ones(size(roll));                       % Initial search direction
    step  = ones(size(roll)) * 10 * pi/180;         % Initial step-size
    
    % Transform to Cartesian coordinates
    C(1,:,:) = cos( phi(1,:,:) ) .* cos( theta(1,:,:) );
    C(2,:,:) = sin( phi(1,:,:) ) .* cos( theta(1,:,:) );
    C(3,:,:) = sin( theta(1,:,:) );
    D = C;                                          % Placeholder for rotated angles
    
    lp = 1;                                         % Loop counter
    upd = true(N,1);                                % Indicator which values nnet to be updated
    while any( upd ) && lp < 50
        
        % Set new roll angle
        rollN(:,:,upd) = roll(:,:,upd) + dir(:,:,upd) .* step(:,:,upd);
        
        R(:) = 0;                                  	% Calculate rotation matrix
        co = cos( rollN(:,:,upd) );
        si = sin( rollN(:,:,upd) );
        R(1,1,upd) = 1;
        R(2,2,upd) = co;
        R(2,3,upd) = -si;
        R(3,2,upd) = si;
        R(3,3,upd) = co;
        
        for n = 1 : N                               % Apply roll angle
            if upd( n )
                D(:,:,n) = R(:,:,n) * C(:,:,n);
            end
        end
        
        % Transform to Spheric coordinates and calculate updated costs
        thetaN(:,:,upd) = atan2(D(3,:,upd),hypot( D(1,:,upd),D(2,:,upd) ));
        phiN(:,:,upd) = atan2(D(2,:,upd),D(1,:,upd));
        asN = qf.calc_angular_spreads( permute(phiN(:,:,upd),[3,2,1]), pow(upd,:) );
        esN = qf.calc_angular_spreads( permute(thetaN(:,:,upd),[3,2,1]), pow(upd,:) );
        
        cstN(upd) = esN ./ asN;
        
        ind = cstN > cst & upd;                     % Results is worse
        dir(:,:,ind)  = -dir(:,:,ind);              % Direction change
        step(:,:,ind) = 0.37 * step(:,:,ind);       % Reduce step size
        
        ind = cstN < cst & upd;                     % Result is better
        roll(:,:,ind) = rollN(:,:,ind);             % Set new roll
        theta(:,:,ind) = thetaN(:,:,ind);           % Store new theta values
        phi(:,:,ind) = phiN(:,:,ind);               % Store new phi values
        cst(ind) = cstN(ind);                       % Update cost
        
        upd = step(:) > 1e-6;                       % Resolution limit
        lp = lp + 1;
    end
    
    orientation(1,1,:) = roll;                      % Save roll angle to orientation
    
    if 0 % Visualize alignment
        % Test angles
        test = [ 0 , -90 : 1 : 90 ];
        nT = numel( test );
        oT = ones( 1,nT );
        
        % Extend variables
        phi   = phi( oT,:,: );
        theta = theta( oT,:,: );
        as = as(:,oT);
        es = es(:,oT);
        
        for t = 2 : nT
            R(:) = 0;                                                   % Calculate rotation matrix
            co = cos( test(t)*pi/180 );
            si = sin( test(t)*pi/180 );
            R(1,1,:) = 1;
            R(2,2,:) = co;
            R(2,3,:) = -si;
            R(3,2,:) = si;
            R(3,3,:) = co;
            
            C(1,:,:) = cos( phi(t,:,:) ) .* cos( theta(t,:,:) ); % Transform to Cartesian coordinates
            C(2,:,:) = sin( phi(t,:,:) ) .* cos( theta(t,:,:) );
            C(3,:,:) = sin( theta(t,:,:) );
            
            for n = 1 : N                                               % Apply pitch angle
                C(:,:,n) = R(:,:,n) * C(:,:,n);
            end
            
            theta(t,:,:) = atan2(C(3,:,:),hypot( C(1,:,:),C(2,:,:) ));  % Transform to Spheric coordinates
            phi(t,:,:)   = atan2(C(2,:,:),C(1,:,:));
            
            as(:,t) = qf.calc_angular_spreads( permute(phi(t,:,:),[3,2,1]), pow );
            es(:,t) = qf.calc_angular_spreads( permute(theta(t,:,:),[3,2,1]), pow );
        end
        imagesc(as./es);
    end
end

orientation = permute( orientation, [1,3,2] );
orientation = angle(exp( 1j*orientation ));

end
