function calc_orientation( h_track )
%CALC_ORIENTATION Calculates the orientation of the mobile device from the positions
%
% Calling object:
%   Object array
%
% Description:
%   This function calculates the orientations of the mobile device based on the snapshot positions.
%   If we assume that the receive array antenna is fixed on a car and the car moves along the
%   track, then the antenna turns with the car when the car is changing direction. This needs to be
%   accounted for when generating the channel coefficients. This function calculates the
%   orientation based on the positions and stores the output in the orientation property of the
%   track object.
%
%   The 3-element orientation vector describes the orientation of the radio device
%   for each position on the track. The reference system for aircraft principal axes is used. The
%   first value describes the "roll angle", i.e. the rotation around an axis drawn through the body
%   of the vehicle from tail to nose in the normal direction of movement. Positive rotation is
%   clockwise (seen from the pilot/drivers perspective). This values is set to 0 when using
%   "calc_orientation". The second value describes the "pitch angle", i.e. the vertical (tilt)
%   angle relative to the horizontal plane; positive rotation is up. The third value describes the
%   bearing or "yaw angle", in mathematic sense. East corresponds to 0, and the angles increase
%   counter-clockwise, so north is 90 degree degrees, south is -90 degree, and west is 180 degree.
%   All values are given in [rad]. Note that by default, QuaDRiGa antennas face east (roll = 0,
%   pitch = 0, yaw = 0).
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

if numel(h_track) > 1
    
    sic = size( h_track );
    prc = false( sic );
    for n = 1 : prod( sic )
        if ~prc( n )
            [ i1,i2,i3,i4 ] = qf.qind2sub( sic, n );
            calc_orientation( h_track(i1,i2,i3,i4) );
            prc( qf.eqo( h_track(i1,i2,i3,i4), h_track ) ) = true;
        end
    end
    
else
    h_track = h_track(1,1);     % Fix for Octave

    if h_track.no_snapshots < 2
        h_track.orientation = [0;0;0];
        
    else
        P = h_track.Ppositions(:,2:end) - h_track.Ppositions(:,1:end-1);
        [ az,el ] = cart2sph(P(1, :), P(2, :), P(3, :));
        
        n = h_track.no_snapshots;
        if h_track.closed
            az(n) = az(1);
            el(n) = el(1);
        else
            az(n) = az(n-1);
            el(n) = el(n-1);
        end
        
        h_track.orientation = [ zeros(1,n) ; el ; az ];
    end
end

end
