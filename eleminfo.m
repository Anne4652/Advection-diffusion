function [elinfo1, elinfo2] = eleminfo(ielem, coord, top, mat, ichois)
% [elinfo1, elinfo2] = eleminfo(ielem, coord, top, mat, ichois)
%
% Wrapper for <elem>_i:
%   [elinfo1, elinfo2] = <elem>_i(ielem, coord, top, mat, ichois)
%
% With <elem> = mat.types(eltype, :), and
%                         eltype = top(ielem, end)
%
% See also elemtype, elemshp, elembld, elemdrv, el1Du_i, el2Du_i, el2Duv_i,
% el2Duvp_i


% get current element type / name
elmname = elemtype(ielem, top, mat);
% append _i
elmname = strcat(elmname, '_i');
% and call element information function
[elinfo1, elinfo2] = feval(elmname, ielem, coord, top, mat.mat, ichois);

% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
end