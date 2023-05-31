function [solution, node2dofMap, freq] = elpasoComputeAndExportSystemSolution(elpasoControlStruct, elpasoProblemHdf5, materials, compute_for_freq)
% elpasoComputeAndExportSystemMatrices - call elpaso and get the system
% response
%
% Inputs:
%   - elpasoControlStruct: elpaso controls to run the elpasoC executable
%   - elpasoProblemHdf5: elpaso hdf5 filename WITH extension
%   - materials: material struct
%   - compute_for_freq: frequency range to compute
%
% Outputs:
%   - solution: Solution after elpaso solve
%   - node2dofmap: Matrix that relates row number to node and respective dof
%   - freq: The frequency steps for which the solution is computed

% create export file
[pathstr, name, ext] = fileparts(elpasoProblemHdf5);

elpasoProblemHdf5_forSolving = [name, '_systemSolve.hdf5'];
eGenOutputFileHdf5 = ['eGenOutput_', name, '_systemSolve', ext];

MD = pwd;
cd(pathstr)
copyfile(elpasoProblemHdf5, elpasoProblemHdf5_forSolving);

% change analysis type to frequency
% use type frequency-basic for eCore and frequency for eResearch
WriteVariableLengthStringToGroup(elpasoProblemHdf5_forSolving, '/Analysis','type','frequency');

%%% change material
elpasoChangeMaterialInElpasoHdf5(elpasoProblemHdf5_forSolving,materials);

% set frequency
h5writeatt(elpasoProblemHdf5_forSolving,'/Analysis','start', min(compute_for_freq));
h5writeatt(elpasoProblemHdf5_forSolving,'/Analysis','steps', length(compute_for_freq));
if length(compute_for_freq) == 1
    h5writeatt(elpasoProblemHdf5_forSolving,'/Analysis','delta', uint64(1));
else
    h5writeatt(elpasoProblemHdf5_forSolving,'/Analysis','delta', uint64((max(compute_for_freq)-min(compute_for_freq))/(length(compute_for_freq)-1)));
end
% compute
tic
[stat, cmdOut] = system([elpasoControlStruct.environment ' && ' elpasoControlStruct.executable ' -c -inp ' elpasoProblemHdf5_forSolving ' ' elpasoControlStruct.flags],'-echo');
exportTimer = toc;

% extract info
[solution,freq,node2dofMap] = elpasoReadSolutionFile(eGenOutputFileHdf5);

cd(MD);

end
