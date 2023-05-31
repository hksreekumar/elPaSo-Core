function WriteVariableLengthStringToGroup(fileName,argpath,attr, value)
% WriteVariableLengthStringToGroup - Writes a variable length string as
% attribute to a group
%
% Syntax: WriteVariableLengthStringToGroup(filename, accesspath)
%
% Inputs:
%   - fileName: elpaso hdf5 filename WITH extension
%   - argpath: the hdf5 path to the group
%   - attr: attribute name to create
%   - value: string value
%
% Outputs:
%   - none

fid = H5F.open (fileName, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
gid  = H5G.open (fid, argpath);

filetype = H5T.copy ('H5T_FORTRAN_S1');
H5T.set_size (filetype,'H5T_VARIABLE');
memtype = H5T.copy ('H5T_C_S1');
H5T.set_size (memtype, 'H5T_VARIABLE');
H5T.set_cset (memtype, 'H5T_CSET_UTF8');

space = H5S.create ('H5S_SCALAR');

H5A.delete(gid,attr);
aid = H5A.create (gid, attr, memtype, space, 'H5P_DEFAULT');
H5A.write (aid, memtype, {value});

H5S.close(space);
H5A.close(aid);
H5G.close(gid);
H5F.close(fid);
end

