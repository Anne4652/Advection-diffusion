function elmname = elemtype(ielem, top, mat)
% elmname = elemtype(ielem, top, mat)
%
% Get element name from mat.types(top(ielem, end), :) if it exists, or
% spawn a warning and return mat.types(1, :)
%
% See also eleminfo, elemshp, elembld, elemdrv


% element type index from topology
eltype = top(ielem, end);

% sanity check for corresponding entry in mat.types(eltype, :)
if eltype > size(mat.types, 1)
    warning(['Element type %i not defined in mat.types ', ...
        '(size(mat.types, 1) = %i), falling back to mat.types(1, :)'], ...
        eltype, size(mat.types, 1));
    eltype = 1;
end

% return element name (type)
elmname = deblank(mat.types(eltype, :));

% part of mlfem_nac: https://gitlab.tue.nl/STEM/mlfem_nac
end