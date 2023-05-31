function [solution,freq,node2dofMap] = elpasoReadSolutionFile(solutionFileHdf5)
%elpasoReadSolutionFile Reads the solution exported in eGenOutput_*.hdf5
%
% Syntax: [solution,freq] = elpasoReadSolutionFile(solutionFileHdf5)
%
% Inputs:
%   - solutionFileHdf5: elpaso output hdf5 filename WITH extension
%
% Outputs:
%   - solution: Solution after elpaso solve
%   - node2dofmap: Matrix that relates row number to node and respective dof
%   - freq: The frequency steps for which the solution is computed


% Get frequency info
nFreq = H5GetNumberOfAttributesFromGroup(solutionFileHdf5, '/Solution/State');

% Read solution
solution = [];
freq = [];
for iFreq = 1:nFreq
    valstruct = h5read(solutionFileHdf5,['/Solution/State/vecFemStep' num2str(iFreq)]);
    val= valstruct.real + 1i*valstruct.imag;
    solution = [solution, val];
    
    currfreq = str2num(cell2mat(h5readatt(solutionFileHdf5,'/Solution/State',['FreqStep' num2str(iFreq)])));
    freq = [freq, currfreq];
end

node2dofMap = h5read(solutionFileHdf5,'/Solution/Maps/mtxDofMap')';

end

