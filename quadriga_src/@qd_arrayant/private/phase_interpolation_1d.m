function a = phase_interpolation_1d( a, b, u )
%PHASE_INTERPOLATION_1D Phase continuous interpolation
%   
%   This functions interpolates the complex values given by [ a, b ] to obtain the values x
%   which lies in between the two valus. The amplitudes are interpolated by linear interpolation
%   and the phases are interpolated using circular linear interpolation. 
%
% Input:
%   a   Value 1 [ N x T ]
%   b   Value 2 [ N x T ]
%   u   Relative distance between a and b ( 0 = a, 1 = b ) [ N x 1 ] or [ N x T ]
%
% Output:
%   x   Interpolated value [ N x T ]
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

% Neagive u
uN = 1-u;

% Calculate absolute values
aa = abs( a );
ba = abs( b );

% Check if none of the values are 0
ia = aa > 1e-6;
aa( ~ia ) = 1e-6;

ib = ba > 1e-6;
ba( ~ib ) = 1e-6;

% Separate real and imaginary parts of the phase
ar = real(a)./aa; 
ai = imag(a)./aa;
br = real(b)./ba; 
bi = imag(b)./ba;

% Do we need circular interpolation for point x ???
xc = ia & ib;                           % Values cannot be 0
xc(xc) = u(xc) > 1e-6;                  % Grid points
xc(xc) = u(xc) < 1-1e-6;

% Only use circular interpolation if there are large phase differences
% The value 0.25 corrsponds to roughly 0.25*180/pi ~ 15 degree
Dab     = (ar(xc)-br(xc)).^2 + (ai(xc)-bi(xc)).^2;
xc(xc)  = Dab > 0.0625;                   

if any(xc(:))
    % Calculate the angle omega between points a and b
    omega = acos( ar(xc).*br(xc) + ai(xc).*bi(xc) );
    
    % Check if the angle is close to pi (phase jump)
    io = omega < pi-1e-6;
    xc(xc) = io;
    
    if any(xc(:))
        % Calculate the amplitude pf the output (linear interpolation)
        xa = uN(xc).*aa(xc) + u(xc).*ba(xc);
        
        % Calculate the weights for the linear combinations
        sinomega  = sin( omega(io) );
        uN(xc) = sin( uN(xc) .* omega(io) ) ./ sinomega ;
        u(xc)  = sin(  u(xc) .* omega(io) ) ./ sinomega ;
        
        % Calculate the real and imaginary parts of point x
        xr = uN(xc).*ar(xc) + u(xc).*br(xc);  % Real part
        xi = uN(xc).*ai(xc) + u(xc).*bi(xc);  % Imaginary part
        
        % Calculate x
        a(xc) = xa.*(xr + 1j*xi);               % Circular interpolation
    end
    xci = ~xc;
    a(xci) = uN(xci).*a(xci) + u(xci).*b(xci);  % Complex linear interpolation
else
    a = uN.*a + u.*b;                           % Complex linear interpolation (all values)
end

end
