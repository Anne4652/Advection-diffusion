% Script to solve the 1D instationary convection diffusion equation
%
%   du/dt + v*du/dx = d/dx (c du/dx) + f
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('sol', 'var')
    % initialise sol
    ndof = max(dest(:));
    if istat == 1 % steady state solution
        sol = zeros(ndof, 1);
    elseif istat == 2 % unsteady solution
        sol = zeros(ndof, ntime + 1);
    end
end
if ~exist('nodfrc', 'var'), nodfrc = []; end


% pos and dest
disp('* create equation numbers ...');
[pos, dest] = equatnr(coord, top, mat);


% global M, C, K and f
disp('* start of matrix-building process ...')

% initialize
ndof = max(dest(:)); % the total number of degrees of freedom
ldof = length(pos(1, :));
% global stiffness matrix (spalloc allocated space for a sparse matrix)
K = spalloc(ndof, ndof, sum(ldof)^2); % global stiffness matrix
M = K; % global mass matrix
C = K; % global damping matrix
rhs = zeros(ndof, 1); % global right-hand-side column


% assemble q and rhs
for ielem = 1:nelem
    % compute the element stiffness matrix and right-hand-side column
    [Me, Ce, Ke, rhse] = elm1dcd(ielem, coord, top, mat.mat);
    
    % assemble the global stiffness matrix and the right-hand-side colomn
    ii = pos(ielem, :);
    K(ii, ii) = K(ii, ii) + Ke;
    M(ii, ii) = M(ii, ii) + Me;
    C(ii, ii) = C(ii, ii) + Ce;
    
    rhs(ii) = rhs(ii) + rhse;
end


% sum & partition and solve
disp('* solving the equations ...')
if istat == 1
    % solve the system of equations, stationairy problem
    q = K + C;
    sol = solvestat(q, rhs, bndcon, [], dest);
    
    % calculate flux (c du/dx)
    cdudx = zeros(size(top, 1), size(top, 2)-2);
    for ielem = 1:size(top, 1)
        cdudx(ielem, :) = elm1dcd_d(ielem, coord, top, mat.mat, ...
            pos, sol);
    end
    
elseif istat == 2
    
    % unsteady solution
    q = M/dt + theta*C + theta*K;
    soln = sol(:, 1);
    
    % calculate flux (c du/dx)
    cdudx = zeros(size(top, 1), size(top, 2)-2, ntime + 1);
    for ielem = 1:size(top, 1)
        cdudx(ielem, :, 1) = elm1dcd_d(ielem, coord, top, ...
            mat.mat, pos, soln);
    end
    
    for itime = 1:ntime
        fprintf('** Start of time step %i\n', itime);
        rhsn = rhs + (M/dt - (1 - theta)*(C + K))*soln;
        sol(:, itime + 1) = solvestat(q, rhsn, bndcon, [], dest);
        soln = sol(:, itime + 1);
        
        % calculate flux (c du/dx)
        for ielem = 1:size(top, 1)
            cdudx(ielem, :, itime + 1) = elm1dcd_d(ielem, coord, top, ...
                mat.mat, pos, soln);
        end
    end
    
    
end

% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
