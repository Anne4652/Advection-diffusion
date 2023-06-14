function [n, dn] = shape(intcoord, inttype)
%
% [n, dn] = shape(intcoord, inttype)
%
% definition of shape functions for quadrilateral elements
%
% input:
%   intcoord  : ksi, eta, zeta coordinates of integration points
%   inttype   : 1   -  constant interpolation
%
%             : 2   -  2-noded (linear) line element
%             : 3   -  3-noded (quadratic) line element
%
%             : 14  -  4-noded (bi-linear) quad (default)
%             : 15  -  5-noded quad (14 with bubble)
%             : 18  -  8-noded (bi-quadratic) quad
%             : 19  -  9-noded (quadratic) quad
%
%             : 23  -  3-noded (linear) triangle
%             : 24  -  4-noded triangle (23 with bubble)
%             : 26  -  6-noded (quadratic) triangle
%             : 27  -  7-noded triangle (26 with bubble)
%
%             : 38  -  8-noded (linear) brick
%             : 320 - 20-noded (quadratic) brick
%
% output:
%   n(1:nint, 1:nlnodes)        : shape function values
%
%   dn(1:nlnodes, 1:nint*nint)  : shape function derivative with respect to
%                                 xi, eta and zeta
%   with:
%       nint    = number of integration points: length(intcoord(:, 1))
%       nlnodes = number of nodes

if nargin==1, inttype = 14; end


% nsd = length(intcoord(1, :));

if inttype == 1
    
    nint = length(intcoord(:, 1));
    n    = ones(nint, 1);
    dn   = zeros(1, nint*2);
    
elseif inttype == 2
    
    [n, dn] = shline2(intcoord);
    
elseif inttype == 3
    
    [n, dn] = shline3(intcoord);
    
elseif inttype == 14
    
    [n, dn] = shquad4(intcoord);
    
elseif inttype == 15
    
    [n, dn] = shquad5(intcoord);
    
elseif inttype == 16
    
    [n, dn] = shquad4p1(intcoord);
    
elseif inttype == 18
    
    [n, dn] = shquad8(intcoord);
    
elseif inttype == 19
    
    [n, dn] = shquad9(intcoord);
    
elseif inttype == 23
    
    [n, dn] = shtri3(intcoord);
    
elseif inttype == 24
    
    [n, dn] = shtri4(intcoord);
    
elseif inttype == 26
    
    [n, dn] = shtri6(intcoord);
    
elseif inttype == 27
    
    [n, dn] = shtri7(intcoord);
    
elseif inttype == 38
    
    [n, dn] = shbrick8(intcoord);
    
elseif inttype == 320
    
    [n, dn] = shbrck20(intcoord);
    
    
else
    error([' This shape function type = ', ...
        num2str(inttype) ' is not supported']);
    
end


% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
end