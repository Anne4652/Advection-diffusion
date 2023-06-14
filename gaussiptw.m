function [ipt, W] = gaussiptw(nipt, nsd, intrule)
%
% [ipt, W] = gaussiptw(nipt, nsd, rule)
%
% Integration points and corresponding weights for gaussian numerical
% integration.
%
% Input:
%   nipt    number of integration points (in one spatial direction)
%   nsd     number of spatial dimensions (e.g. 1D, 2D, 3D)
%   rule    (optional) integration rule:
%           (1) Gauss-Legendre integration points & weights
%               (default, nipt <= 5)
%           (2) Gauss-Lobatto integration points & weights
%               (forces ipts at the element end points, nipt <= 7)
%           (3) Gaussian integration points & weights for triangular
%               elements (nsd == 2  only, nipt <= 7)
% Output:
%   ipt     location of integration points in parent element
%           (in local coordinates, -1:1)
%   W       weight factors corresponding to to the integration points in
%           ipt
%
% See also elemshp

% set defaults
if nargin < 3, intrule = 1; end
if nargin < 2, intrule = 3; end % single argument must be a triangle call


if intrule == 1
    % default Gaussian-Legendre integration points and weights
    
    if nipt == 1
        
        ksi1d    = 0;
        weight1d = 2;
        
    elseif nipt == 2
        
        ksi1d    = [ -1/sqrt(3) 1/sqrt(3) ];
        weight1d = [      1         1     ];
        
    elseif nipt == 3
        
        ksi1d    = [ -sqrt(3/5)  0  sqrt(3/5) ];
        weight1d = [     5/9    8/9     5/9   ];
        
    elseif nipt == 4
        
        a  = sqrt(3/7 + 2/7*(sqrt(6/5)));
        b  = sqrt(3/7 - 2/7*(sqrt(6/5)));
        wa = (18 - sqrt(30))/36;
        wb = (18 + sqrt(30))/36;
        ksi1d    = [ -a -b  b a  ];
        weight1d = [ wa wb wb wa ];
        
    elseif nipt == 5
        
        a  = 1/3*sqrt(5 + 2*sqrt(10/7));
        b  = 1/3*sqrt(5 - 2*sqrt(10/7));
        w0 = 128/255;
        wa = (322 - 13*sqrt(70))/900;
        wb = (322 + 13*sqrt(70))/900;
        ksi1d    = [ -a  -b  0  b  a ];
        weight1d = [ wa  wb w0 wb wa ];
        
    else
        
        error(['%ind order Gauss-Legendre integration (intrule = 1)', ...
            ' not implemented'], nipt)
    end
elseif intrule == 2
    % Gauss-Lobatto integration points and weights (ipts at the nodes)
    
    if nipt == 1
        
        ksi1d    = 0;
        weight1d = 2;
        
    elseif nipt == 2
        
        ksi1d    = [ -1 1 ];
        weight1d = [  1 1 ];
        
    elseif nipt == 3
        
        ksi1d    = [-1   0    1];
        weight1d = [1/3 4/3 1/3];
        
    elseif nipt == 4
        
        ksi1d    = [-1 -sqrt(1/5) sqrt(1/5) 1];
        weight1d = [1   5         5         1]/6;
        
    elseif nipt == 5
        
        ksi1d    = [-1 -sqrt(3/7) 0 sqrt(3/7) 1];
        weight1d = [1/10 49/90 32/45 49/90 1/10];
        
    elseif nipt == 6
        
        a  = sqrt(1/3 - (2*sqrt(7)/21) );
        aw = (14 + sqrt(7))/30;
        b  = sqrt(1/3 + (2*sqrt(7)/21) );
        bw = (14 - sqrt(7))/30;
        
        ksi1d    = [-1 -b -a a b 1];
        weight1d = [1/15 bw aw aw bw 1/15];
        
    elseif nipt == 7
        
        a  = sqrt(5/11 - 2*sqrt(5/3)/11 );
        aw = (124 + 7*sqrt(15))/350;
        b  = sqrt(5/11 + 2*sqrt(5/3)/11 );
        bw = (124 - 7*sqrt(15))/350;
        
        ksi1d    = [-1 -b -a 0 a b 1];
        weight1d = [1/21 bw aw 256/525 aw bw 1/21];
        
    else
        
        error(['%ind order Gauss-Lobatto integration (intrule = 2)', ...
            ' not implemented'], nipt)
    end
    
elseif intrule == 3
    % Gaussian integration points and weights for triangular elements
    W = zeros(nipt, 1);
    
    if nipt == 1
        
        ipt = [1/3 1/3 1/3]';
        W   = 1;
        
    elseif nipt == 3
        
        ipt = ( ones(3, 3) + eye(3) )/6;
        W   = [1/3 1/3 1/3]';
        
    elseif nipt == 4
        
        ipt = 1/15*[%
            5     5     5
            9     3     3
            3     9     3
            3     3     9];
        
        W   = 1/48*[-27 25 25 25]';
        
    elseif nipt == 6
        
        a1 = 0.816847572980459;
        a2 = 0.091576213509771;
        b1 = 0.108103018168070;
        b2 = 0.445948490915965;
        
        ipt = [%
            a1 a2 a2
            a2 a1 a2
            a2 a2 a1
            b1 b2 b2
            b2 b1 b2
            b2 b2 b1];
        
        W(1:3) = 0.109951743655322;
        W(4:6) = 0.223381589678011;
        
    elseif nipt == 7
        
        a  = 1/3;
        a1 = 0.0597158717;
        a2 = 0.4701420641;
        b1 = 0.7974269853;
        b2 = 0.1012865073;
        
        ipt = [%
            a  a  a
            a1 a2 a2
            a2 a1 a2
            a2 a2 a1
            b1 b2 b2
            b2 b1 b2
            b2 b2 b1];
        
        W(1)   = 0.225;
        W(2:4) = 0.1323941527;
        W(5:7) = 0.1259391805;
        
    else
        
        error(['%ind order Gauss integration rule not implemented', ...
            ' for triangular elements (intrule = 3)'], nipt)
        
    end
    
end


if intrule < 3
    % extend to multiple dimensions & collect output
    
    if nsd == 1
        
        ipt = ksi1d';
        W   = weight1d;
        
    elseif nsd == 2
        
        weight2d = weight1d'*weight1d;
        W        = weight2d(:)';
        
        eta = ones(nipt,1)*ksi1d;
        ksi = eta';
        ipt = [ksi(:), eta(:)];
        
    elseif nsd == 3
        
        weight2d = weight1d'*weight1d;
        weight2d = weight2d(:);
        weight3d = weight2d*weight1d;
        W        = weight3d(:)';
        
        dzeta  = ones(nipt,1)*ksi1d;
        eta    = dzeta';
        dzeta  = ones(nipt,1)*(dzeta(:)');
        dzeta  = dzeta(:);
        
        eta    = ones(nipt, 1)*(eta(:)');
        ksi    = eta';
        eta    = eta(:);
        ksi    = ksi(:);
        
        ipt    = [ksi eta dzeta];
        
    end
else
    W = W/2; % triangle weights should sum to 1/2?
end


% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
end