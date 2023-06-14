
function [dndx,detjac,jac,jaci]=getderiv(dndksi,nodcoord,ip);
%
% [dndx,detjac,jac,jaci]=getderiv(dndksi,nodcoord,ip);
%
% input:
%   dndksi    : derivatives of the shape function n with respect to local coord.
%                 dndksi(1:nlnodes,1:nint*nsd),
%                   nlnodes   : number of nodes of the element
%                   nint      : number if integration points
%                   nsd       : space dimension = length(coord(1,:));
%   nodcoord  : coordinates of the nodal points
%                 nodcoord(1:nlnodes,1:nsd)
%   ip        : integration point number
%
% output:
%   dndx     : derivatives of the shape function n with respect to global coor.
%               in integration point ip
%                 dndx(1:nlnodes,1:nsd)
%   detjac    : determinant of the jacobian matrix
%   jac       : jacobian of the transformation matrix from local to gloval coor.
%   jaci      : inv(jac)

nsd = length(nodcoord(1,:));

k=nsd*(ip-1);

dn(:,1:nsd) = dndksi(:,k+1:k+nsd);

jac         = dn'*nodcoord; 
detjac      = det(jac); 
jaci        = inv(jac);

dndx        = dn*jaci';




% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
