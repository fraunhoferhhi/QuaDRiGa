function [ mu, epsilon, gamma, Ri, sigma, kappa, delta ] =...
    log2dfit( vi, ai, bi, round_digits, show_plot, a_linear, b_linear, a_name, b_name, v_name )
%LOG2DFIT Logarithmic curve fitting with two variables
%
% Calling object:
%   None (static method)
%
% Description:
%   This function fits the input data "vi" sampled at points "ai" and "bi" to the following model
%   in a least square sense:
%
%       V = R + X * S  
%       R = mu + epsilon * log10( a ) + gamma * log10( b ) 
%       S = sigma + kappa * log10( a ) + delta * log10( b ) 
%
%   V is assumed to be a random, normal-distributed variable with a mean value R and standard
%   deviation S. The reference value R depends on the variables a and b. X is a normal-distributed
%   random variable with zero-mean and unit-variance. The scaling of the STD of V is done by S,
%   which also depends on the variables a and b. 
%
% Input:
%   vi
%   Vector of input data values
%
%   ai
%   Sample points for variable a (linear). If ai is empty, 1D fitting is done for bi only.
%
%   bi
%   Sample points for variable b (linear). If bi is empty, 1D fitting is done for ai only. If ai
%   and bi are empty, only mu and sigma are returned.
%
%   round_digits
%   Rounds the output to this number of decimal digits after the coma.
%
%   show_plot
%   Shows a plot with the results and outputs the fitted values to the console. The plot contains
%   the average and the 1-sigma interval above and below the average as well as the data points.
%   The options are: (0) disables the plot; (1, default) shows the plot and the text for avg and
%   std; (2) shows the plot and the text for avg only; (3) shows the plot and the text for std
%   only; (4) shows the plot only
%
%   a_linear
%   If set to true, fitting is done for linear values of ai instead of log10( ai ).
%
%   b_linear
%   If set to true, fitting is done for linear values of bi instead of log10( bi ).
%
%   a_name
%   Alternative variable name for a in the figure title and console text
%
%   b_name
%   Alternative variable name for a in the figure title and console text
%
%   v_name
%   Alternative name for the z-axis in the figure
%
% Output:
%   mu
%   Reference value R at a = 1 and b = 1
%
%   epsilon
%   Scaling of the reference vale R with a [R/log10(a)] or [R/a]
%
%   gamma
%   Scaling of the reference value R with b [R/log10(b)] or [R/b]
%
%   Ri
%   The reference value R for the input sample points ai and bi
%
%   sigma
%   The STD of V at a = 1 and b = 1
%
%   kappa
%   Scaling of the STD S with a [S/log10(a)] or [S/a]; If kappa is not requested as an output
%   variable, no STD scaling is assumed.
%
%   delta
%   Scaling of the STD S with b [S/log10(b)] or [S/b]; If delta is not requested as an output
%   variable, STD scaling is only done for a.
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

% Check if ai is given
if exist( 'ai','var' ) && ~isempty( ai ) && any(abs(ai-mean(ai)) > 1e-3)
    use_ai = true;
else
    use_ai = false;
    ai = ones(size(vi));
end

% Check if bi is given
if exist( 'bi','var' ) && ~isempty( bi ) && any(abs(bi-mean(bi)) > 1e-3)
    use_bi = true;
else
    use_bi = false;
    bi = ones(size(vi));
end

if exist( 'a_linear','var' ) && ~isempty( a_linear )
    a_linear = logical( a_linear(1) );
else
    a_linear = false;
end

if exist( 'b_linear','var' ) && ~isempty( b_linear )
    b_linear = logical( b_linear(1) );
else
    b_linear = false;
end

if ~exist( 'a_name','var' ) || isempty( a_name )
    a_name = 'a';
end
if ~exist( 'b_name','var' ) || isempty( b_name )
    b_name = 'b';
end
if ~exist( 'v_name','var' ) || isempty( v_name )
    v_name = 'v';
end

if ~exist( 'round_digits','var' ) || isempty( round_digits )
    num_format = '%g';
else
    num_format = ['%1.',num2str(round_digits),'f'];
end

if ~exist( 'show_plot','var' ) || isempty( show_plot )
    show_plot = false;
end

% Check if "kappa" is requested
if nargout > 5 && use_ai
    calc_kappa = true;
else
    calc_kappa = false;
end

% Check if "delta" is requested
if nargout > 6 && use_bi
    calc_delta = true;
else
    calc_delta = false;
end

no_vals = numel( vi );

vi = reshape( vi , no_vals, 1 );
ai = reshape( ai , no_vals, 1 );
bi = reshape( bi , no_vals, 1 );

% Sort out errors
ii = ~isinf(vi) & ~isnan(vi) & ~isinf(ai) & ~isnan(ai) & ~isinf(bi) & ~isnan(bi);
vi = vi(ii);
ai = ai(ii);
bi = bi(ii);

% Logarithmic variables
if a_linear
    lai = ai;
else
    lai = log10( ai );
end
if b_linear
    lbi = bi;
else
    lbi = log10( bi );
end

if use_bi && use_ai
    %cst = [];
    
    % Low-Res Global search
    gamma = -50:5:50;
    xn = zeros(numel(gamma),1);
    for iC = 1 : numel(gamma)
        tmp = polyfit( lai, vi-gamma(iC)*log10( bi ) ,1);
        epsilon = tmp(1);
        mu = tmp(2);
        v = vi - epsilon*lai - mu - gamma(iC)*lbi;
        xn(iC) = sum( abs( v ).^2 );
        %cst(end+1) = xn(iC);
    end
    [x,iC] = min(xn);
    gamma = gamma(iC);
    
    % Use the SAGE algorithm fo fit the exponential function to the data
    a   = gamma;    % Initial value
    dm  = 1;    % Step size
    delta = Inf; ddir = +1; lp = 2;
    while lp<100 && delta > 1e-9
        if lp>1; an = a + ddir * dm; delta = abs(a-an); else an = a; end
        
        tmp = polyfit( lai, vi-an*log10( bi ) ,1);
        epsilon = tmp(1);
        mu = tmp(2);
        v = vi - epsilon*lai - mu - an*lbi;
        
        xn = sum( abs( v ).^2 );
        %cst(end+1) = xn;
        
        if xn < x; a = an; x = xn; else ddir = -ddir; dm = 0.2 * dm; end
        lp = lp + 1;
    end
    gamma = an;
    
    %plot(10*log10(cst - cst(end) + 1e-9))
    
elseif use_ai
    tmp = polyfit( lai, vi ,1);
    epsilon = tmp(1);
    mu = tmp(2);
    gamma = 0;
    
elseif use_bi
    tmp = polyfit( lbi, vi ,1);
    epsilon = 0;
    mu = tmp(2);
    gamma = tmp(1);
    
else
    epsilon = 0;
    mu = mean(vi);
    gamma = 0;
    
end

if a_linear && b_linear
    Ri = mu + epsilon * ai  + gamma *  bi;
elseif a_linear
    Ri = mu + epsilon * ai + gamma * log10( bi );
elseif b_linear
    Ri = mu + epsilon * log10( ai ) + gamma * bi;
else
    Ri = mu + epsilon * log10( ai ) + gamma * log10( bi );
end

di = abs(vi-Ri).^2;  % Variance
if calc_kappa && calc_delta
    
    % Fit the variance
    [ sigma, kappa, delta ] = qf.log2dfit( di, ai, bi, [], 0, a_linear, b_linear );
    
    x = rand(1,10000)*(max(ai)-min(ai)) + min(ai);
    y = rand(1,10000)*(max(bi)-min(bi)) + min(bi);
    if a_linear && b_linear
        z = sigma + kappa.*x + delta.*y;
    elseif a_linear
        z = sigma + kappa.*x + delta.*log10(y);
    elseif b_linear
        z = sigma + kappa.*log10(x) + delta.*y;
    else
        z = sigma + kappa.*log10(x) + delta.*log10(y);
    end
    z(z<0) = 0;
    z = sqrt(z);
    
    % Fit the STD
    [ sigma, kappa, delta ] = qf.log2dfit( z, x, y, [], 0, a_linear, b_linear );  
    
elseif calc_kappa
    
    % Fit the variance
    [ sigma, kappa ] = qf.log2dfit( di, ai, [], [], 0, a_linear );
    
    % Transform variance to STD
    x = (0:0.0001:1)*(max(ai)-min(ai)) + min(ai);
    if a_linear
        z = sigma + kappa.*x;
    else
        z = sigma + kappa.*log10(x);
    end
    z(z<0) = 0;
    z = sqrt(z); 
    
    % Fit the STD
    [ sigma, kappa ] = qf.log2dfit( z, x, [], [], 0, a_linear );
    delta = 0;
    
elseif calc_delta 
    
    % Fit the variance
    [ sigma, delta ] = qf.log2dfit( di, bi, [], [], 0, b_linear );  
    
    % Transform variance to STD
    x = (0:0.0001:1)*(max(bi)-min(bi)) + min(bi);
    if b_linear
        z = sigma + delta.*x;
    else
        z = sigma + delta.*log10(x);
    end
    z(z<0) = 0;
    z = sqrt(z); 
    
    % Fit the STD
    [ sigma, delta ] = qf.log2dfit( z, x, [], [], 0, b_linear );
    kappa = 0;
    
else
    sigma = std( vi-Ri );
    delta = 0;
    kappa = 0;
end

if exist( 'round_digits','var' ) && ~isempty( round_digits )
    mu = round( mu * 10^round_digits )/10^round_digits;
    epsilon = round( epsilon * 10^round_digits )/10^round_digits;
    gamma = round( gamma * 10^round_digits )/10^round_digits;
    sigma = round( sigma * 10^round_digits )/10^round_digits;
    kappa = round( kappa * 10^round_digits )/10^round_digits;
    delta = round( delta * 10^round_digits )/10^round_digits;
end

str ={};
str{1} = ['avg = ',num2str(mu,num_format),' '];
if epsilon > 0.001
    str{1} = [str{1},'+ ',num2str(epsilon,num_format)];
elseif epsilon < -0.001
    str{1} = [str{1},'- ',num2str(-epsilon,num_format)];
end
if a_linear && abs( epsilon ) > 0.001
    str{1} = [str{1},' * ',a_name,' '];
elseif ~a_linear && abs( epsilon ) > 0.001
    str{1} = [str{1},' * log10( ',a_name,' ) '];
end
if gamma > 0.001
    str{1} = [str{1},'+ ',num2str(gamma,num_format)];
elseif gamma < -0.001
    str{1} = [str{1},'- ',num2str(-gamma,num_format)];
end
if b_linear && abs( gamma ) > 0.001
    str{1} = [str{1},' * ',b_name,' '];
elseif ~b_linear && abs( gamma ) > 0.001
    str{1} = [str{1},' * log10( ',b_name,' ) '];
end
str{1} = str{1}(1:end-1);

str{2} = ['std = ',num2str(sigma,num_format),' '];
if kappa > 0.001
    str{2} = [str{2},'+ ',num2str(kappa,num_format)];
elseif kappa < -0.001
    str{2} = [str{2},'- ',num2str(-kappa,num_format)];
end
if a_linear && abs( kappa ) > 0.001
    str{2} = [str{2},' * ',a_name,' '];
elseif ~a_linear && abs( kappa ) > 0.001
    str{2} = [str{2},' * log10( ',a_name,' ) '];
end
if delta > 0.001
    str{2} = [str{2},'+ ',num2str(delta,num_format)];
elseif delta < -0.001
    str{2} = [str{2},'- ',num2str(-delta,num_format)];
end
if b_linear && abs( delta ) > 0.001
    str{2} = [str{2},' * ',b_name,' '];
elseif ~b_linear && abs( delta ) > 0.001
    str{2} = [str{2},' * log10( ',b_name,' ) '];
end
str{2} = str{2}(1:end-1);


if show_plot
    
    % Display fitted function on the command line
    if show_plot == 1 || show_plot == 2
        disp(str{1});
    end
    if show_plot == 1 || show_plot == 3
        disp(str{2});
    end
    
    if use_bi && use_ai
        tmp = 0.1*(max(lai)-min(lai));  as = min(lai)-tmp : tmp : max(lai)+1.5*tmp;
        tmp = 0.1*(max(lbi)-min(lbi));  bs = min(lbi)-tmp : tmp : max(lbi)+1.5*tmp;
        [as,bs] = meshgrid( as , bs  );
        
        vs  = mu + epsilon * as + gamma * bs;
        vsu = vs + sigma + kappa * as + delta * bs;
        vsl = vs - sigma - kappa * as - delta * bs;
        
        if ~a_linear
            as = 10.^as;
        end
        if ~b_linear
            bs = 10.^bs;
        end
        
        surf( as, bs, vs, 'facealpha',1)
        hold on
        surf( as, bs, vsu,'facealpha',0.3,'LineStyle',':')
        surf( as, bs, vsl,'facealpha',0.3,'LineStyle',':')
        plot3( ai,bi,vi,'.m' )
        hold off
        grid on
        
        title(str)
        xlabel(a_name);
        ylabel(b_name);
        zlabel(v_name);
        
    elseif use_ai
        
        tmp = 0.1*(max(lai)-min(lai));  as = min(lai)-tmp : tmp : max(lai)+1.5*tmp;
        vs  = mu + epsilon * as;
        vsu = vs + sigma + kappa * as;
        vsl = vs - sigma - kappa * as;
        if ~a_linear
            as = 10.^as;
        end
        plot( ai,vi,'.m' )
        hold on
        plot( as,vs,'-+k','Linewidth',3 )
        plot( as,vsu,'--^k','Linewidth',2 )
        plot( as,vsl,'--vk','Linewidth',2 )
        hold off
        grid on
        title(str)
        xlabel(a_name);
        ylabel(v_name);
        legend('Data','AVG','STD')
        
    elseif use_bi
        
        tmp = 0.1*(max(lbi)-min(lbi));  bs = min(lbi)-tmp : tmp : max(lbi)+1.5*tmp;
        vs  = mu + gamma * bs;
        vsu = vs + sigma + delta * bs;
        vsl = vs - sigma - delta * bs;
        if ~b_linear
            bs = 10.^bs;
        end
        plot( bi,vi,'.m' )
        hold on
        plot( bs,vs,'-+k','Linewidth',3 )
        plot( bs,vsu,'--^k','Linewidth',2 )
        plot( bs,vsl,'--vk','Linewidth',2 )
        hold off
        grid on
        title(str)
        xlabel(b_name);
        ylabel(v_name);
        legend('Data','AVG','STD')
        
    end
end

end
