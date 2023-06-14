function cdudx = el1Du_d(ielem, coord, top, mat, pos, sol)
% cdudx = el1Du_d(ielem, coord, top, mat, pos, sol)
%
% Flux c du/dx for 1D convection diffusion elements
%
% input:
%
%  ielem    : element number
%  coord    : global node coordinates array
%  top      : global topology array
%  mat      : material properties (diffusion constant, c = mat(imat, 1);
%  pos      : global pos array
%  sol      : array with calculated concentrations (solutions)
%
% output:
%   cdudx   : row with flux c du/dx for each node of the element
%
% See also elm1d, elm1dcd, fem1d, fem1dcd


% get the coordinates of the nodes of this element
nodcoord = coord(top(ielem, 1:end-2), :);

% pick up some material parameters
imat = top(ielem, end-1);
c = mat(imat, 1);

% element topology
etop = top(ielem,:);

% position of the local degrees of freedom in the global solution array
[~, vshp]        = el1Du_i( 1,[], etop, mat, 5);
[n, dndxi, intw] = elemshp( [], 1, etop, mat, vshp(1, :));

% nodal solutions
unod = sol(pos(ielem, :));

% allocate for flux at the integration points
cdudxint = zeros(2, 1);
for ipt = 1:length(intw)
    
    % get shape function derivatives with respect to global coordinates
    dndxy = getderiv(dndxi, nodcoord, ipt);
    
    % compute the strain displacement matrix
    b = dndxy(:, 1)';
    
    % stress (c du/dx) at this integration point
    cdudxint(ipt) = c*b*unod;
    
end

% cdudxint can be extrapolated to the nodes
cdudx = (n\cdudxint)';

% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
end
