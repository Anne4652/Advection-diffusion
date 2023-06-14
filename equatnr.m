function [pos,dest,tsed] = equatnr(coord,top,material,ichois);
%
% EQUATNR returns FEM-bookkeeping arrays
%
% [pos,dest,tsed] = equatnr(coord,top,material,ichois);
%
% On input, COORD and TOP must be defined as
%
% COORD(INODE,ISD)    = X     : coordinate of node INODE in dimension ISD
% TOP(IELEM,ILNODE)   = INODE : global node number INODE for each element IELEM
% TOP(IELEM,NLNODE+1) = IMAT  : material type of element IELEM
% TOP(IELEM,NLNODE+2) = ITYPE : element type of element IELEM
%
% On output POS and DEST are filled according to
%
% POS[ielem,ielq]   = IEQ : global equation number for local equation IELQ of
%                           element IELEM
% DEST(inode,ildof) = IEQ : global equation number for local dof ILDOF in
%                           node INODE
%
% Where : INODE  = global node number
%         IEQ    = global equation number
%         IELEM  = element number
%         ILNODE = local node number in an element
%         NLNODE = maximum number of local nodes in an element
%         ILDOF  = local degree of freedom number in a node
%         IELQ   = local equation number in an element
%
% system function
%
% External calls:
%  sysinfo
%  eleminfo


[nnodes,nsd,nelem,maxnlnodes,iimat,iitype] = sysinfo(size(coord),size(top));
if nargin<4
   ichois=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POS & DEST array 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%** initialize maximum number of local DOF's to 1 (will be increased 
%   later on if required)
maxnldof = 1; 
dest    = zeros(nnodes,maxnldof);

%** initialize global equation number
eqnr = 0;

for ielem=1:nelem

    %** get local nodes and local DOFs
    connect = top(ielem,1:maxnlnodes);
    lnodes  = find(connect);
    nlnodes = length(lnodes);
    nodes   = top(ielem,lnodes);
    ldof    = eleminfo(ielem,coord,top,material,ichois(1));
    if length(ichois)>1
       ldof = ichois(2)*ones(size(ldof));
    end
    
    %** update maximum number of local DOF's 
    maxnldofnew = max(ldof);
    if maxnldofnew > maxnldof
       dest     = [dest zeros(nnodes,maxnldofnew-maxnldof) ];
       maxnldof = maxnldofnew;
    end
     
    %** create global equation numbers
    for ilnode=1:nlnodes
        for ildof=1:ldof(ilnode)

            ignode = nodes(ilnode);
            if dest(ignode,ildof) == 0 % otherwise use previous equation number
                eqnr               = eqnr + 1;
                dest(ignode,ildof) = eqnr;
                tsed(eqnr,1:2)     = [ ignode ildof ];
            end
            jdof            = sum(ldof(1:ilnode)) - ldof(ilnode) + ildof;
            pos(ielem,jdof) = dest(ignode,ildof);            

        end;
    end;

end;

ndof = eqnr;

% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
