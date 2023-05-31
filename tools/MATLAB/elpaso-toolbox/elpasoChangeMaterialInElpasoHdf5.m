function elpasoChangeMaterialInElpasoHdf5(file,materials)
% ChangeMaterialInElpasoHdf5 - Changes the material parameters in the hdf5
% file
%
% Syntax: ChangeMaterialInElpasoHdf5(file,materials)
%
% Inputs:
%   - file: elpaso hdf5 filename without extension
%   - materials: array of material struct holding the various material
%   types and their parameters
%
% Outputs:
%   - none

file = erase(file,'.hdf5');
filenametemp = [file , '_temp.hdf5' ];
filename = [file , '.hdf5' ];
copyfile(filename,filenametemp);

fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
% Change the material by first deleting the attributes
for i = 1:length(materials)
    dsid = H5D.open(fid,['/Materials/material',num2str(materials(i).Id)]);
    if strcmp(materials(i).type,'STR_LIN_VIS_ISO_DIR')
        H5A.delete(dsid,'E');
        H5A.delete(dsid,'rho');
        H5A.delete(dsid,'t');
        H5A.delete(dsid,'nu');
        H5A.delete(dsid,'eta');
    elseif strcmp(materials(i).type,'STR_LIN_ELA_ISO_DIR')
        H5A.delete(dsid,'E');
        H5A.delete(dsid,'rho');
        H5A.delete(dsid,'t');
        H5A.delete(dsid,'nu');
    elseif strcmp(materials(i).type,'AF_LIN_UAF_ISO_DIR')
        H5A.delete(dsid,'c');
        H5A.delete(dsid,'rho');
    elseif strcmp(materials(i).type,'isotropic') || strcmp(materials(i).type,'STR_LIN_ELA_ISO_DIR')
        H5A.delete(dsid,'E');
        H5A.delete(dsid,'rho');
        H5A.delete(dsid,'t');
        H5A.delete(dsid,'nu');
    else
        error('Material unknown')
    end
    
    H5D.close(dsid);    
end
H5F.close(fid);

% change values
for i = 1:length(materials)
    if strcmp(materials(i).type,'STR_LIN_VIS_ISO_DIR')
        WriteVariableLengthStringToDataSet(filename,['/Materials/material',num2str(materials(i).Id)],'E', num2str(materials(i).E));
        WriteVariableLengthStringToDataSet(filename,['/Materials/material',num2str(materials(i).Id)],'rho', num2str(materials(i).rho));
        WriteVariableLengthStringToDataSet(filename,['/Materials/material',num2str(materials(i).Id)],'t', num2str(materials(i).t));
        WriteVariableLengthStringToDataSet(filename,['/Materials/material',num2str(materials(i).Id)],'nu', num2str(materials(i).nu));
    elseif strcmp(materials(i).type,'STR_LIN_ELA_ISO_DIR')
        WriteVariableLengthStringToDataSet(filename,['/Materials/material',num2str(materials(i).Id)],'E', num2str(materials(i).E));
        WriteVariableLengthStringToDataSet(filename,['/Materials/material',num2str(materials(i).Id)],'rho', num2str(materials(i).rho));
        WriteVariableLengthStringToDataSet(filename,['/Materials/material',num2str(materials(i).Id)],'t', num2str(materials(i).t));
        WriteVariableLengthStringToDataSet(filename,['/Materials/material',num2str(materials(i).Id)],'nu', num2str(materials(i).nu));
    elseif strcmp(materials(i).type,'AF_LIN_UAF_ISO_DIR')
        WriteVariableLengthStringToDataSet(filename,['/Materials/material',num2str(materials(i).Id)],'c', num2str(materials(i).cf));
        WriteVariableLengthStringToDataSet(filename,['/Materials/material',num2str(materials(i).Id)],'rho', num2str(materials(i).rhof));
    else
        error('Material unknown')
    end   
end

end

