
function [n,dn]=shline2(intcoord);
%
% [n,dn] = shline2(intcoord);
%
% shape function routine for 2 noded line element
%
% input:
%   intcoord :  local coordinates (ksi,eta) of integration points
%
% output:
%   n        :  shape functions at integration points
%   dn       :  derivative of shape functions with resp. to ksi and eta
%

nint    = length(intcoord(:,1));
nsd     = length(intcoord(1,:));
nlnodes = 2;
n       = zeros(nint,nlnodes);
dn      = zeros(nlnodes,nint*nsd);


if nsd~=1,
   error('This shape function only valid for 1D problems');
end

for ip=1:nint
%
% get intcoord (k1) and eta (k2) value
%
  k1 = intcoord(ip,1); 
%
% n(i,j): at integration point i value of shape function of node j
%
  n(ip,1) =-0.5*(k1-1);
  n(ip,2) = 0.5*(k1+1);
%
  k = 2*ip-1;
%
% dn(i,j): derivative of shape fuction of node i in the direction
%               intcoord and eta at all int. points
%
  dn(1,ip) = -0.5;
  dn(2,ip) = +0.5;

end

% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
