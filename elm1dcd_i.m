function [par1, par2] = elm1dcd_i(ielem, ~, top, mat, ichois)
% [par1, par2] = elm1dcd_i(ielem, ~, top, mat, ichois)
%
% Wrapper for el1Du_i:
%   [par1, par2] = el1Du_i(ielem, 0, top, mat, ichois)
%
% input:
%   ielem   : element number
%   ~       : not used (coord, nodal coordinates)
%   top     : mesh topology array
%   mat     : material properties (mat.mat), mat(iimat, 5) = ietype
%
% output:
%   ichois = 1  par1: number of degrees of freedom for each node
%               par2 = [];
%   ichois = 3  par1: sequence of nodes for plotting purposes
%               par2 = [];
%   ichois = 4  par1: shape function that defines the geometry
%               par2 = [];
%   ichois = 5  par1: definition of location of dofs (vpos)
%               par2: definition of their shape functions (vshp)
%
% See also el1Du_i, elm1dcd, elm1d

% make the call
[par1, par2] = el1Du_i(ielem, 0, top, mat, ichois);

% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
end
