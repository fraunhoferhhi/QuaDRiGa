function [ ds, mean_delay ] = calc_delay_spread( taus, pow )
%CALC_DELAY_SPREAD Calculates the delay spread in [s]
%
% Description:
%   This function calculates the RMS delay spread from a given set of delays and powers. It is used
%   by the "qd_builder" to map the path powers to delays.
%
% Input:
%   ang
%   A vector of deays [s]. Dimensions: [ n_taus x n_path ]
%
%   pow
%   A vector of path powers in [W]. Dimensions: [ n_taus x n_path ]
%
% Output:
%   as
%   The RMS delay spread for each delay vector. Dimensions: [ n_taus x 1 ]
%
%   mean_angle
%   The mean delay in [s]. Dimensions: [ n_taus x 1 ]
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

N = size( taus,1 );
if size( pow,1 ) < N
    pow = pow( ones(1,N), : );
end

% Normalize powers
pt = sum( pow,2 );
pow = pow./pt( :,ones(1,size(pow,2))  );

mean_delay = sum( pow.*taus,2 ); 

tmp = taus - mean_delay(:,ones( 1,size(taus,2) ) );

ds = sqrt( sum(pow.*(tmp.^2),2) - sum( pow.*tmp,2).^2 );

end
