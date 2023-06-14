function [par1, par2] = el1Du_i(ielem, ~, top, mat, ichois)
% [par1, par2] = el1Du_i(ielem, ~, top, mat, ichois)
%
% Element information for 1D elements with 1 DoF (u)
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
% See also elm1d, elm1dcd

% get some element info
imat   = top(ielem, end-1);
ietype = mat(imat, 5);
nnodes = numel(nonzeros(top(ielem, 1:end-2)));

% initialise
par1 = []; par2 = [];

if ichois == 1
    % definition of number of degrees of freedom for each node
    
    if ietype == 1
        % linear element
        par1 = [1 1];
    elseif ietype == 2
        % quadratic element
        par1 = [1 1 1];
    end
    
    % sanity check
    npar = length(par1);
    if npar ~= nnodes
        string_err = ['Number of nodes = ' num2str(nnodes) ...
            ' for element type ' num2str(ietype) ...
            ' in ielem ' num2str(ielem) ', should be '];
        error([string_err num2str(npar)]);
    end
    
    
elseif ichois == 2
    
    if ietype == 1
        % linear element
        par1 = 3*ones(1, nnodes);
    elseif ietype == 2
        % quadratic element
        par1 = 3*ones(1, nnodes);
    end
    
    
elseif ichois == 3
    % definition of sequence of nodes for plotting purposes
    
    if ietype == 1
        % linear element
        par1 = [1 2];
    elseif ietype == 2
        % quadratic element
        par1 = [1 2; 2 3];
    end
    
    
elseif ichois == 4
    % shape function that defines the geometry
    
    if ietype == 1
        % linear element
        gshp = [2 1 2];
    elseif ietype == 2
        % quadratic element
        gshp(1, :)   = [3 1 3];
    end
    par1 = gshp;
    
    
elseif ichois == 5
    % definition of location of dofs (vpos) and their shape functions
    
    if ietype == 1
        % linear element
        vpos       = zeros(1, 2);
        vpos(1, :) = 1:2;
        
        vshp       = zeros(1, 3);
        vshp(1, :) = [2 1 2];
        
    elseif ietype == 2
        % quadratic element
        vpos       = zeros(1, nnodes);
        vpos(1, :) = 1:3;
        
        vshp       = zeros(1, 3);
        vshp(1, :) = [3 1 3];
    end
    
    par1 = vpos; par2 = vshp;
    
else
    
    error([' Ichois = ' num2str(ichois) ' not supported ']);
    
end


% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
end
