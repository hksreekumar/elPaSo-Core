function [noa] = H5GetNumberOfAttributesFromGroup(filename, grouppath)
%H5GetNumberOfAttributesFromGroup Returns the number of attributes in a group
%
% Syntax: H5GetNumberOfAttributesFromGroup(filename, accesspath)
%
% Inputs:
%   - fileName: hdf5 filename WITH extension
%   - grouppath: the hdf5 path to the group
%
% Outputs:
%   - noa: the number of attributes


fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
gid  = H5G.open (fid, grouppath);

info = H5G.get_info(gid);
noa = info.nlinks;

H5G.close(gid);
H5F.close(fid);

end

