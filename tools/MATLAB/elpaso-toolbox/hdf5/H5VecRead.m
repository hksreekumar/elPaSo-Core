function [vec] = H5VecRead(filename, accesspath)
% H5VecRead - Reads a complex vector from HDF5
%
% Syntax: [vec] = H5VecRead(filename, accesspath)
%
% Inputs:
%   - filename: elpaso hdf5 filename WITH extension
%   - accesspath: the hdf5 path to the vector dataset
%
% Outputs:
%   - vec: read vector

%display(['Reading ' fileprefix '...'])
valstruct = h5read(filename,accesspath);

if(~length(valstruct.imag))
    %display('real data');
else
    %display('complex data');
    vec= valstruct.real.' + 1i*valstruct.imag.';
end

end

