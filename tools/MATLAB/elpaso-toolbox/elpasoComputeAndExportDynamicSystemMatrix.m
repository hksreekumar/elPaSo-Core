function [KDYN] = elpasoComputeAndExportDynamicSystemMatrix(elpasoControlStruct, elpasoProblemHdf5, materials, compute_for_freq)
% elpasoComputeAndExportDynamicSystemMatrix - call elpaso and get the dynamic system
% matrix
%
% Inputs:
%   - elpasoControlStruct: elpaso controls to run the elpasoC executable
%   - elpasoProblemHdf5: elpaso hdf5 filename WITH extension
%   - materials: material struct
%   - compute_for_freq: frequency step to compute
%
% Outputs:
%   - KDYN: Dynamic stiffness matrix

% create export file
[pathstr, name, ext] = fileparts(elpasoProblemHdf5);

elpasoProblemHdf5_forExport = [name, '_systemExport.hdf5'];
eGenSystemFileHdf5 = ['eGenSystem_', name, '_systemExport', ext];

MD = pwd;
cd(pathstr)
copyfile(elpasoProblemHdf5, elpasoProblemHdf5_forExport);

% change analysis type to frequency
% use type frequency-basic for eCore and frequency for eResearch
WriteVariableLengthStringToGroup(elpasoProblemHdf5_forExport, '/Analysis','type','frequency'); 

% set frequency
h5writeatt(elpasoProblemHdf5_forExport,'/Analysis','start', compute_for_freq);

%%% change material
elpasoChangeMaterialInElpasoHdf5(elpasoProblemHdf5_forExport,materials);

% compute
tic
[stat, cmdOut] = system([elpasoControlStruct.environment ' && ' elpasoControlStruct.executable ' -c -inp ' elpasoProblemHdf5_forExport ' -exportSysH5 1' ' ' elpasoControlStruct.flags],'-echo');
exportTimer = toc;

% extract info
KDYN = H5CSRToSparseMatFOM(eGenSystemFileHdf5, '/DynamicStiffness');

cd(MD);

end
