function [ az, el, J, pow, accuracy, resolved ] = calc_angles( cf, h_qd_arrayant, no_subpath, range, verbose )
%CALC_ANGLES Estimates the arrival angles from the channel coefficients
%
% Description:
%   This function estimates the arrival angles from a given coefficient matrix cf. This is done by
%   comparing the coefficients with the antenna response of the given antenna object h_qd_arrayant
%   and determining the most likely arrival direction in azimuth and elevation relative to the
%   local antenna coordinate system. Note that the direction finding results depend on the spatial
%   resolution capabilities of the antenna. If the array antenna has no spatial resolution (e.g. a
%   single dipole), then the returned values will be incorrect.
%
% Input:
%   cf
%   The channel coefficients; size [ n_rx, n_tx, n_coeff ]
%
%   h_qd_arrayant
%   The receive antenna object. The number of elements must match  [ n_rx ]
%
%   no_subpath
%   The maximum number of sub-paths per path, default: 1
%
%   range
%   The allowed angle range in degree [ -el +el -az +az ], Default: [-90 90 -180 180]
%
%   verbose
%   Show the estimation progress. 0 = Disable, 1 = Show progress bar, 2 = Show spatial signature
%
% Output:
%   az
%   The azimuth angles in [rad], dimension [ n_coeff, no_subpath ]
%
%   el
%   The elevation angles in [rad], dimension [ n_coeff, no_subpath ]
%
%   J
%   The Jones vectors of the path, dimension [ 2, n_tx, n_coeff, no_subpath ]
%
%   pow
%   The power of the path without antenna pattern [ n_tx, n_coeff, no_subpath ]
%
%   accuracy
%   The estimation MSE [ n_coeff, 1 ]
%
%   resolved
%   The resolved spatial power [ n_coeff, 1 ]
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

if ~exist('no_subpath','var') || isempty( no_subpath )
    no_subpath = 1;
end

if ~exist( 'range','var' ) || isempty( range )
    range = [-90 90 -180 180];
end

if ~exist('verbose','var') || isempty( verbose )
    verbose = 1;
end

% Progress bar
tStart = clock;
if verbose==1; fprintf('Calc Angles  ['); end; m0=0;

% Number of paths
L = size(cf,3);

% Number of transmit antennas
no_tx = size( cf,2 );

% Check h_qd_arrayant 
if exist('h_qd_arrayant','var') && ~isempty( h_qd_arrayant ) ...
        && isa(h_qd_arrayant,'qd_arrayant') && numel(h_qd_arrayant) == 1
    if any( abs(h_qd_arrayant.element_position(:) ) > 1e-6)
        h_qd_arrayant = h_qd_arrayant.copy;
        h_qd_arrayant.combine_pattern;
    end
else
    error('QuaDRiGa:QF:calc_angles:h_qd_arrayant','??? "h_qd_arrayant" must be a qd_arrayant object.')
end

% Range mask
el_mask = h_qd_arrayant.elevation_grid*180/pi + 1e-4 > range(1) & ...
    h_qd_arrayant.elevation_grid*180/pi - 1e-4 < range(2);
if ~any( el_mask )
    [~,tmp] = min( abs(h_qd_arrayant.elevation_grid*180/pi - mean(range(1:2))) );
    el_mask(tmp) = true;
end
az_mask = h_qd_arrayant.azimuth_grid*180/pi + 1e-4 > range(3) & ...
    h_qd_arrayant.azimuth_grid*180/pi - 1e-4 < range(4);
if ~any( az_mask )
    [~,tmp] = min( abs(h_qd_arrayant.azimuth_grid*180/pi - mean(range(3:4))) );
    az_mask(tmp) = true;
end

elevation_grid = single( h_qd_arrayant.elevation_grid( el_mask ) );
azimuth_grid   = single( h_qd_arrayant.azimuth_grid(az_mask) );

no_az = sum( az_mask );
no_el = sum( el_mask );
no_rx = h_qd_arrayant.no_elements;

element_position = single( h_qd_arrayant.element_position );

F_theta_vec = reshape( single( h_qd_arrayant.Fa(el_mask,az_mask,:) ), [], no_rx );
F_phi_vec   = reshape( single( h_qd_arrayant.Fb(el_mask,az_mask,:) ), [], no_rx );

% Allocate output variables
az = zeros( L,no_subpath ) ;
el = zeros( L,no_subpath ) ;
J  = zeros( 2,no_tx,L,no_subpath ) ;
pow = zeros( no_tx,L,no_subpath ) ;
accuracy = Inf( L,1 ) ;
resolved = zeros( L,1 ) ;

% Placeholder for the antenna responses
F = complex( zeros( no_rx, 2 ,'single') );

for l = 1 : L
    if verbose==1; m1=ceil(l/L*50); if m1>m0; for m2=1:m1-m0; fprintf('o'); end; m0=m1; end; end;
    
    % The [nr x nt] coefficient matrix in time domain
    G  = single( cf(:,:,l) );
    
    % The power for each tx antenna
    pG = sum(abs(G).^2,1);
    
    % Placeholder for the reconstructed spatial sub-paths
    Ghat = complex( zeros( no_rx, no_tx, no_subpath ,'single') );
    
    for m = 1 : no_subpath
        
        if m == 1
            Gtilde = G;
            pGtilde = pG;
        else
            % Remove already detected sub-paths
            Gtilde = G - sum(Ghat,3);
            pGtilde = sum(abs(Gtilde).^2,1); % Numeric precision issue
        end
        
        % Initialize Jones Vector
        Js = ones(2,no_tx,'single')./sqrt(2);
        
        % First, use the initial Jones vector to calculate an initial estimate of the angles, then
        % estimate the correct Jones vector from the pattern and repeat the angle estimation.
        for iJ = 1:2
            % Calculate combined pattern using channel coefficients as
            % beamforming weights
            conjG      = conj( Gtilde ); % Saves time
            Fhat_theta = F_theta_vec * conjG;
            Fhat_phi   = F_phi_vec * conjG;
            
            % Estimate angles
            Fhat_sq = abs( Fhat_theta ).^2 * abs(Js(1,:).').^2 + abs( Fhat_phi ).^2 * abs(Js(2,:).').^2;
            Fhat_sq = Fhat_sq - min(Fhat_sq);
            Fhat_sq = Fhat_sq./max(Fhat_sq);
            Fhat_sq = reshape( Fhat_sq , no_el , no_az , [] );
            
            [ xEl, xAz ] = maxi( Fhat_sq, elevation_grid, azimuth_grid, 0.01 );
            
            % Calculate pattern response for the extimated angles
            [ F(:,1),F(:,2) ] = h_qd_arrayant.interpolate(xAz,xEl,...
                1:no_rx,azimuth_grid, elevation_grid,F_theta_vec,F_phi_vec,element_position);
            
            % Estimate Jones vector
            [U,S,V] = svd(F,0);
            Si = S;
            if S(2,2) > no_rx * eps(S(1))
                S([1,4]) = 1./S([1,4]);
            else
                S(1) = 1./S(1);
                S(4) = 0;
            end
            pinvF = V*S*U';
            
            Js = pinvF * Gtilde;
            Pmpc   = sum( abs(Js).^2 , 1 );
            Js     = Js./([1;1]*sqrt(Pmpc));
            Js(1,pGtilde == 0) = 1;
            Js(2,pGtilde == 0) = 0;
            
            if 10*log10(Si(1) / Si(4)) > 10
                break
            end
        end
        
        % Reconstructed pattern response from the estimated angles and Jones vector
        Ghat(:,:,m) = ones(no_rx,1)*sqrt(Pmpc) .* (F * Js);
        
        tmp = abs( G-sum(Ghat,3) ).^2;
        MSE = sum(tmp(:)) / sum(abs(G(:)).^2);
        resolved_tmp = sum(abs(Ghat(:)).^2) / sum(abs(G(:)).^2);
        
        % Multiple detections
        if m==1 || ( MSE < accuracy( l ) && MSE > 0)
            % Save results
            el( l, m )     = double( xEl );
            az( l, m )     = double( xAz );
            J( :,:,l, m )  = double( Js );
            pow( :, l, m ) = double( Pmpc );
            accuracy( l )  = double( MSE );
        else
            Ghat(:,m) = 0;  
            m = m - 1; %#ok<FXSET>   
            break
        end
        
        if resolved_tmp >= 0.95
            break
        end
    end
    resolved( l ) = resolved_tmp;
    
    if verbose == 2
        
        [~,ii] = max(sum(pow( :, l, m ),3));

        if ~exist('han1','var')
            han1 = figure;
        end
              
        % Spatial fingerprint
       
        conjG      = conj( G ); % Saves time
        Fhat_theta = F_theta_vec * conjG;
        Fhat_phi   = F_phi_vec * conjG;
            
        Fhat_sq = sum(abs( Fhat_theta ).^2 + abs( Fhat_phi ).^2,2);
        Fhat_sq = Fhat_sq - min(Fhat_sq);
        Fhat_sq = Fhat_sq./max(Fhat_sq);
        Fhat_sq = reshape( Fhat_sq , no_el , no_az , [] );
        
        figure(han1)
        imagesc( azimuth_grid*180/pi, elevation_grid*180/pi, Fhat_sq )
        hold on
        for mm = 1 : m
            plot(az(l,mm)*180/pi,el(l,mm)*180/pi,'sk','Markerfacecolor','m','Markersize',10)
        end
        hold off
        set(gca,'Ydir','Normal');
        title('Spatial signature')
        xlabel('Azimuth angle [deg]')
        xlabel('Elevation angle [deg]')
            
        if ~exist('han2','var')
            han2 = figure;
        end
        
        figure(han2)
        xxx = [  G(:,ii) , sum(Ghat(:,ii,:),3) ];
        plot(xxx.','-bo')
        hold on
        plot(G(:,ii),'o','Markerfacecolor','r')
        hold off
        axis([-1 1 -1 1]*max(abs(xxx(:)))*1.1)
        title(['Coefficients ; MSE ',num2str(10*log10(accuracy(l)),'%1.1f'),' dB'])
        
        disp(' ')
        disp(['Path ',num2str(l),' : MSE ',num2str(10*log10(accuracy(l)),'%1.1f'),' dB; Resolved  ',...
            num2str( resolved( l ) * 100,'%1.1f'),' %'])
        for iPp = 1:m
            POW = 10*log10(mean(pow(:,l,iPp),1));
            str = ['Subpath ',num2str(iPp),' : Pow. ',num2str( POW ,'%1.1f'),...
                ' dB, El. ',num2str(el( l,iPp )*180/pi,'%1.1f'),...
                ' deg, Az. ',num2str(az( l,iPp )*180/pi,'%1.1f'),' deg'];
            disp(str)
        end
        
        if l < L
            pause
        end
    end
end

if verbose == 1
    fprintf('] %5.0f seconds\n',round( etime(clock, tStart) ));
end

end

% =========== Subfunction: MAXI ==========================
function [ el, az ] = maxi( P, elevation_grid, azimuth_grid, res )
% MAXI Spline interpolation between angle estimates to find the most lilely angle

% Interpolation targets
res = 1/ceil(1/res);
u   = 0:res:1;
nu  = numel(u);

nEl = numel(elevation_grid);
nAz = numel(azimuth_grid);

if abs(azimuth_grid(1)+pi) < 1e-6 && abs(azimuth_grid(end)-pi) < 1e-6
    wrap_around = true;
else
    wrap_around = false;
end

% Get initial Values
[~,ind] = max( P(:) );
[iEl,iAz] = ind2sub(size(P),ind);

% Get the neighboring elevation positions
pEl1 = int16( [ iEl-2 ; iEl-1 ; iEl ; iEl+1 ] );
pEl2 = int16( [ iEl-1 ; iEl ; iEl+1 ; iEl+2 ] );
if pEl1(1) < 1
    pEl2 = pEl2-pEl1(1)+1;
    pEl1 = pEl1-pEl1(1)+1;
elseif pEl2(end) > nEl
    pEl1 = pEl1 - (pEl2(end)-nEl);
    pEl2 = pEl2 - (pEl2(end)-nEl);
end
pEl1( pEl1 > nEl ) = nEl;
pEl2( pEl2 > nEl ) = nEl;
pEl1( pEl1 < 1 )   = 1;
pEl2( pEl2 < 1 )   = 1;
vEl = elevation_grid(pEl1(2:4));

% Get the neighboring azimuth positions
pAz1 = int16( [ iAz-2 ; iAz-1 ; iAz ; iAz+1 ] );
pAz2 = int16( [ iAz-1 ; iAz ; iAz+1 ; iAz+2 ] );
if pAz1(1) < 1
    if wrap_around
        pAz1(pAz1==0)  = nAz-1;
        pAz1(pAz1==-1) = nAz-2;
        pAz2(pAz2==0)  = nAz-1;
    else
        pAz2 = pAz2-pAz1(1)+1;
        pAz1 = pAz1-pAz1(1)+1;
    end
elseif pAz2(end) > nAz
    if wrap_around
        pAz1(pAz1==nAz+1) = 2;
        pAz2(pAz2==nAz+1) = 2;
        pAz2(pAz2==nAz+2) = 3;
    else
        pAz1 = pAz1 - (pAz2(end)-nAz);
        pAz2 = pAz2 - (pAz2(end)-nAz);
    end
end
pAz1( pAz1 > nAz ) = nAz;
pAz2( pAz2 > nAz ) = nAz;
pAz1( pAz1 < 1 )   = 1;
pAz2( pAz2 < 1 )   = 1;
vAz = unwrap(azimuth_grid( pAz1(2:4) ));

% Kernel for the spline interpolation function
K = 0.5 .* ...
    [0,  2,  0,  0; ...
    -1,  0,  1,  0; ...
     2, -5,  4, -1; ...
    -1,  3, -3,  1].';

% Get the offset between the actual user position and the next point
% on the map.
uK = K * [ ones(1,numel(u)) ; u ; u.^2 ; u.^3];

% Interpolate
Pi = zeros(2*nu-1);
Pi( 1:nu , 1:nu  ) = uK.' * P(pEl1,pAz1) * uK;
Pi( 1:nu ,nu:end ) = uK.' * P(pEl1,pAz2) * uK;
Pi(nu:end, 1:nu  ) = uK.' * P(pEl2,pAz1) * uK;
Pi(nu:end,nu:end ) = uK.' * P(pEl2,pAz2) * uK;

% Find the maximum
[~,ind] = max( Pi(:) );
[iEli,iAzi] = ind2sub(size(Pi),ind);

% Determine final angle
if vEl(1) ~= vEl(3)
    v  = vEl(1) : mean(diff(vEl))/(nu-1) : vEl(3);
    el = v(iEli);
else
    el = vEl(2);
end

if vAz(1) ~= vAz(3)
    v = vAz(1) : mean(diff(vAz))/(nu-1) : vAz(3);
    az = angle(exp(1j*v(iAzi)));
else
    az = vAz(2);
end

end
