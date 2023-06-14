function [Me, Ce, Ke, rhse] = elm1dcd(ielem, coord, top, mat)
%
% linear convection diffusion element
%
% du/dt + a du/dx = d/dx(c du/dx) + f
%
%
% Input
%  ielem    : element number
%  coord    : global node coordinates array
%  top      : global topology array
%  mat      : material properties
%
% Output
%  Me      : element mass matrix
%  Ce      : element damping matrix
%  Ke       : element stifness matrix
%  rhse     : element right hand side
%
% See also fem1dcd, elm1dcd_i, elm1dcd_s, elm1dcd_d


% get the number of nodes on this element
nlnodes = length(nonzeros(top(ielem, 1:end-2)));
% get the coordinates of the nodes of this element
nodcoord = coord(top(ielem, 1:nlnodes), :);

% pick up some material parameters
imat = top(ielem, nlnodes + 1);
c = mat(imat, 1);     % diffusion coefficient
a = mat(imat, 2);     % convection velocity
f = mat(imat, 3);     % volume force
fsupg = mat(imat, 4); % if fsupg ~= 0, use supg, otherwise galerkin

% element topology
etop = top(ielem, :);

% position of the local degrees of freedom in the global solution array
[vpos, vshp]     = el1Du_i( 1, [], etop, mat, 5);
[n, dndxi, intw] = elemshp([], coord, etop, mat, vshp(1, :));

nint = length(intw);
ndofu = length(nonzeros(vpos));

% initialize Ke, Me and Ce and rhse
Ke = zeros(ndofu, ndofu);
Ce = Ke;
Me = Ke;

rhse = zeros(ndofu, 1);
if fsupg > 0
    h = max(nodcoord)-min(nodcoord);
    beta  = abs(a)*h/(2*c);
    alpha = coth(beta)-1/beta;
end

% loop over integration points
for int = 1:nint
    
    % get shape function derivatives with respect to global coordinates
    [dndxy, detj] = getderiv(dndxi, nodcoord, int);
    
    % compute the strain displacement matrix
    b = dndxy(:, 1)';
    
    % diffusion part
    Ke = Ke + b'*c*b*intw(int)*detj;
    
    % convection part
    if fsupg > 0
        Ce = Ce + (n(int, :)' + ...
            h*alpha*a/(2*abs(a))*dndxy(:, 1))*a*dndxy(:, 1)'*intw(int)*detj;
    else
        Ce = Ce + n(int, :)'*a*dndxy(:, 1)'*intw(int)*detj;
    end
    
    % mass matrix
    if fsupg > 0
        Me = Me + (n(int, :)' + ...
            h*alpha*a/(2*abs(a))*dndxy(:, 1))*n(int, :)*intw(int)*detj;
    else
        Me = Me + n(int, :)'*n(int, :)*intw(int)*detj;
    end
    
    % apply volume forces
    rhse = rhse + n(int, :)'*f*intw(int)*detj;
    
end

% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
end