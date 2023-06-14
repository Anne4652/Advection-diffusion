function [N, dN, intw, intcoord, detj] = elemshp(ielem, coord, top, ...
    ~, shp, intcoord, idet)
% [N, dN, intw, intcoord, detj] = elemshp(ielem, coord, top, ~, shp, ...
%    intcoord, idet)
%
% input:
%   ielem    : element number
%   coord    : nodal coordinates
%   top      : mesh topology array
%   ~        : mat structure (not used)
%   shp      : array with shape function options
%             shp(1) inttype: integration type parameter for shape
%             shp(2) integration rule (generates ipts for missing intcoord)
%             shp(3) nint: number of integration points
%             shp(4) (optional) overwrites nsd = size(coord, 2)
%   intcoord : (optional) local coordinates of integration points
%   idet     : (optional) Jacobian
%
% output:
%   N        : rows with shape function values
%   dN       : columns with shape function derivative values
%   intw     : weighing factors for the integration points
%   intcoord : (local) coordinates for the integration points
%   detj     : value of Jacobian per integration point
%
% See also shape, gaussiptw, getderiv

% parse arguments
if nargin < 6, intcoord = []; end
if nargin < 7, idet = []; end
if isempty(ielem), ielem = 1; end

% parse shp
shfunction = shp(1);
intrule    = shp(2);
nint       = shp(3);
if numel(shp) > 3
    % overwrite nsd
    nsd = shp(4);
else
    nsd = size(coord, 2);
end

% initialise detj & intw
detj = [];
intw = [];


if isempty(intcoord)
    % get integration points
    [intcoord, intw] = gaussiptw(nint, nsd, intrule);
end
% get shape functions and derivatives
[N, dN] = shape(intcoord, shfunction);

% get Jacobian
if ~isempty(idet)
    % nodal coordinates
    nodcoord = coord(nonzeros(top(ielem, 1:end-2)), :);
    % number of integration points
    nipt = size(intcoord, 1);
    % allocate
    detj = zeros(1, nipt);
    for ipt = 1:nipt
        [~, detj(ipt)] = getderiv(dN, nodcoord, ipt);
    end
end

% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
end
