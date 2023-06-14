function cdudx = elm1dcd_d(ielem, coord, top, mat, pos, sol)
% cdudx = elm1dcd_d(ielem, coord, top, mat, pos, sol)
%
% Wrapper for el1Du_d:
%   cdudx = el1Du_d(ielem, coord, top, mat, pos, sol)
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
% See also el1Du_d, elm1d, elm1dcd

% make the call
cdudx = el1Du_d(ielem, coord, top, mat, pos, sol);

% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
end