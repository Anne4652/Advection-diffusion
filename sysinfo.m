function [nnodes, nsd, nelem, maxnlnodes, iimat, iitype] = ...
    sysinfo(sizecoord, sizetop)
% [nnodes, nsd, nelem, maxnlnodes, iimat, iitype] = ...
%     sysinfo(sizecoord, sizetop)
%
% SYSINFO returns general system information
%
% [nnodes,nsd,nelem,maxnlnodes,iimat,iitype] = sysinfo(size(coord),size(top));
%
% Input :
% SIZE(COORD) = dimensions of the COORD array
% SIZE(TOP)   = dimensions of the TOP   array
%
% Output :
% NNODES      = total number of nodes in the mesh
% NSD         = dimension of the coordinate space
% NELEM       = total number of elements in the mesh
% MAXNLNODES  = maximum number of nodes in an element
% IIMAT       = location of material type in TOP array
% IITYPE      = location of element  type in TOP array
%
% system function

nnodes     = sizecoord(1);
nsd        = sizecoord(2);
nelem      = sizetop(1);
maxnlnodes = sizetop(2) - 2;
iimat      = maxnlnodes + 1;
iitype     = maxnlnodes + 2;



% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
end