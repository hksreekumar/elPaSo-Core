function [KMAT, MMAT] = elpasoComputeAndExportSystemMatrices(elpasoControlStruct, elpasoProblemHdf5, materials)
% elpasoComputeAndExportSystemMatrices - call elpaso and get the system
% matrices
%
% Inputs:
%   - elpasoControlStruct: elpaso controls to run the elpasoC executable
%   - elpasoProblemHdf5: elpaso hdf5 filename WITH extension
%   - materials: material struct
%
% Outputs:
%   - KMAT: Stiffness matrix
%   - MMAT: Mass matrix

% create export file
[pathstr, name, ext] = fileparts(elpasoProblemHdf5);

elpasoProblemHdf5_forExport = [name, '_systemExport.hdf5'];
eGenSystemFileHdf5 = ['eGenSystem_', name, '_systemExport', ext];

MD = pwd;
cd(pathstr)
copyfile(elpasoProblemHdf5, elpasoProblemHdf5_forExport);

% change analysis type to eigen
WriteVariableLengthStringToGroup(elpasoProblemHdf5_forExport, '/Analysis','type','eigen');

%%% change material
elpasoChangeMaterialInElpasoHdf5(elpasoProblemHdf5_forExport,materials);

% compute
tic
[stat, cmdOut] = system([elpasoControlStruct.environment ' && ' elpasoControlStruct.executable ' -c -inp ' elpasoProblemHdf5_forExport ' -exportSysH5 1' ' ' elpasoControlStruct.flags],'-echo');
exportTimer = toc;

% extract info
KMAT = H5CSRToSparseMatFOM(eGenSystemFileHdf5, '/Stiffness');
MMAT = H5CSRToSparseMatFOM(eGenSystemFileHdf5, '/Mass');

cd(MD);

end
