function [sol] = solvestat(stif,rhs,bndcon,nodfrc,dest,solpar)

% sol = solvestat(stif,rhs,bndcon,nodfrc,dest,solpar);
%
% returns the solution of a static problem using direct solver
%
% solpar    : solpar(1)=1, (default) boundary condtions are put in sol(ip)
%             solpar(1)~=1, sol(ip)=zeros(size(ip));
%
% Modifications:
% 03-11-95:   check for emptiness of ip

if (nargin==5), solpar=1; end

iter=solpar(1);

%** get global paramaters
ndof  = max(max(dest));
nsol  = length(rhs(1,:));
sol   = zeros(ndof,nsol);

%** create index arrays
[iu,ip] = indexbnd(bndcon,dest);

%** add nodal forces **
if ~isempty(nodfrc)
   rhs = rhs + mkvec(nodfrc,dest,0)*ones(1,nsol);
end

%** fill essential boundary conditions
if ~isempty(bndcon)
   bndvec    = mkvec(bndcon,dest,0)*ones(1,nsol);
   if (iter==1), 
       sol(ip,:) = bndvec(ip,:);
   else
       sol(ip,:) = zeros(size(ip));
   end
end

%* solve **
if ~isempty(ip),
   [n,m]=size(rhs);
   disp(['Condition number of matrix: ' num2str(1/condest(stif(iu,iu))) ]);
   for j=1:m
      rhs(iu,j)=rhs(iu,j)-stif(iu,ip)*sol(ip,j);
   end
   
   sol(iu,:) = stif(iu,iu) \ rhs(iu,:);
else
   sol = stif\rhs;
end

return
% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
