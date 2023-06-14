function [iu,ip] = indexbnd(bndcon,dest)

% INDEXBND returns indices of prescribed and unknown dof's
%
% [iu,ip] = indexbnd(bndcon,dest)
%
% creates index arrays IU and IP, containing the equation numbers
% of the unknown (IU) and prescrided (IP) DOFs
%
% Input
% BNDCON = definition of essential boundary conditions
% DEST   = global equation numbers per nodal point
%
% Output
% IU     = equation numbers of unkown (free) DOFs
% IP     = equation numbers of unkown (free) DOFs
%
% system function

[nbndcon,dummy] = size(bndcon);
ndof            = max(max(dest));
fixed           = zeros(ndof,1);

if nbndcon>0
  for ibc =1:nbndcon
    node     = bndcon(ibc,1);
    dof      = bndcon(ibc,2);
    P        = dest(node,dof);
    if P~=0
       fixed(P) = ibc;
    end
  end; %loop over boundary conditions %
end

%** find prescribed DOFs 
ip = find(fixed);

%** find free DOFs 
iu = find(~fixed);

% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
