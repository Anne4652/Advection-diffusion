function [vec]=mkvec(def,dest,t)

% MKVEC generates a vector containing boundary conditions or nodal
% forces.
% 
% [vec]=mkvec(def,dest,t);
%
% VEC     = vector to be created
% DEF     = boundary condition or nodal force definition array
% DEST    = destination array (see DESTIN)
% T   (0) = time level
%
% Array DEF must be filled on input according to :              
%   b        = boundary condition/nodal force number        
%   def(b,1) = INODE    : global node number                       
%   def(b,2) = idof     : nodal degree of freedom number       
%   def(b,3) = amplitud : amplitude of load
%   def(b,4) = ichois   : type of load time dependency
%   def(b,5) = tau      : characteristic time of load time-dependency
%
% On output VEC = AMPLITUD * FUNC(T/TAU);
%
%  ichois = 1 :  func(t) = 1                        constant/step
%  ichois = 2 :  func(t) = t                        ramp
%  ichois = 3 :  func(t) = 1      0 <= t <=1        pulse
%                func(t) = 0      t > 1 
%  ichois = 4 :  func(t) = 2*t    0   <= t <= 0.5   triangle
%                func(t) = 2-2*t  0.5 <  t <= 1
%                func(t) = 0      t > 1
%
% sytem function

if nargin<3, t=0; end

ndof = max(max(dest));
ndef = length(def(:,1));

vec=zeros(ndof,1);

if length(def(1,:))==3  % old type of (constant) boundary definition

  for idef=1:ndef
    A    = def(idef,1);
    idof = def(idef,2);
    ipos = dest(A,idof);
    if ipos~=0
       vec(ipos) = def(idef,3);
    end
  end

else   % new type (time dependent) boundary condition

  for idef=1:ndef
    A    = def(idef,1);
    idof = def(idef,2);
    ipos = dest(A,idof);
    amplit = def(idef,3);
    ichois = def(idef,4);
    if ichois==1
       tau = 1;
    else
       tau    = def(idef,5);
    end;
    vec(ipos) = amplit * mkfunc(ichois,t/tau);
  end

end


% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
