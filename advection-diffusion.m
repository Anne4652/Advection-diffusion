%%%%%%%%%%%%%%%%%%%%%%%% demo of fem1dcd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% This demo solves the 1D convection diffusion problem:                 %
%                   du/dt + v du/dx = d/dx(c du/dx) + f                 %
%                                                                       %
% The demo uses fem1dcd as engine with the element elm1dcd.             %
%                                                                       %
% The problem that is solved is Example 15.1 from the book:             %
%            Biomechanics: Concepts and Computation                     %
% See also demo_2d_mimic_conv_diffusion for a `2D' implementation       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all  %#ok<CLALL>

istat = 1;      % 1: steady state problem, 2: unsteady problem
dt    = 0.01; 	% magnitude of time step
ntime = 10;  	% number of time steps
theta = 0.5;   	% theta parameter of theta-time integration scheme

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Pre-processing                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%    DEFINE THE MESH    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
xmin = 0; xmax = 1; % domain = [xmin xmax]
nelem = 5;	        % number of elements
norder = 1;         % element order (1 - linear, 2 - quadratic)

% build coord and top
if norder == 1
    % linear elements, 2 nodes / element
    dx = (xmax-xmin)/nelem;
    coord = (xmin:dx:xmax)';
    top = [(1:nelem)' (2:nelem + 1)' ones(nelem, 2)];
elseif norder == 2
    % quadratic elements, 3 nodes / element
    dx = (xmax-xmin)/(2*nelem);
    coord = (xmin:dx:xmax)';
    top = [(1:2:nelem*2)' (2:2:(nelem*2 + 1))', ...
        (3:2:(nelem*2 + 2))' ones(nelem, 2)];
end


%%%%%%%%%%%%%%%%%% DEFINE THE MATERIAL PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%
mat.mat(1) = 0.02;		   % c: diffusion coefficient
mat.mat(2) = 0.5;  	   % v: convective velocity
mat.mat(3) = 0;		   % f: source term
mat.mat(5) = norder;   % element order (1 - linear, 2 - quadratic)
% define element type(s)
mat.types = 'elm1dcd'; % element type

%%%%%%%%%%%%%%%%%%%%%%   Initial Condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sol = zeros(size(coord, 1), 1);


%%%%%%%%%%%%%%%%% DEFINE THE BOUNDARY CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%
bndcon = [1        1 0;
    size(coord, 1) 1 1];
nodfrc = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     Create and solve equations                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% Main program FEM1DCD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fem1dcd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Post-processing                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% post process the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = mat.mat(1);
v = mat.mat(2);
h = xmax/nelem;
Pe = (v*(xmax-xmin))/c; % Peclet number for title

analitinis=(1-exp(v*coord/c))/(1-exp(v*xmax/c));


figure 1
xlabel('x [-]')
ylabel('u [-]')
plot(coord,analitinis,'r-p',coord, sol,'-*');legend('exact solution','numerical solution');
grid on
title(['Pe = ', num2str(Pe)])

error = max(analitinis - sol)
% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
