function [ as, mean_angle ]  = calc_angular_spreads( ang, pow, wrap_angles )
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

if ~exist('wrap_angles','var')
    wrap_angles = true;
end

N = size( ang,1 );
if size( pow,1 ) < N
    pow = pow( ones(1,N), : );
end

% Normalize powers
pt = sum( pow,2 );
pow = pow./pt( :,ones(1,size(pow,2))  );

if wrap_angles
    mean_angle = angle( sum( pow.*exp( 1j*ang ) , 2 ) ); % [rad]
else
    mean_angle = sum( pow.*ang,2 ); 
end

phi = ang - mean_angle(:,ones( 1,size(ang,2) ) );

if wrap_angles
    phi = angle( exp( 1j*phi ) );
end

as = sqrt( sum(pow.*(phi.^2),2) - sum( pow.*phi,2).^2 );

end
