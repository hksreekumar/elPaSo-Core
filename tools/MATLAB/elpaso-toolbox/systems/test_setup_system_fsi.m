function [s] = test_setup_system_fsi(p)

% define elpaso controls for the system
elpasoControls.executable = '/home/sreekumar/software/repos/elPaSo-Research/bin/elpasoC';
elpasoControls.flags = '';
elpasoControls.environment = 'export LD_LIBRARY_PATH=none && source /software/intel2020/compilers_and_libraries/linux/bin/compilervars.sh -arch intel64 -platform linux'; % IMPORTANT!

% define problem 
elpasoProblemHdf5 = '/mnt/scratch/sreekumar/elpaso_temp/plateCavityModel/plateCavityModel_freq.hdf5'; % add the absolute path
s.domains = [1890,8505]; %[<dof in plate>, <dof in fluid>]

% define materials for p in [0,1]^d
n_d = length(p);

%%% material 1
materials(1).name = 'steel';
materials(1).type = 'STR_LIN_ELA_ISO_DIR';
materials(1).Id = 1;     
materials(1).E = 210e9*(1 + 0.2*(p(1)-0.5));     % sE(1);
materials(1).nu = 0.3*(1 + 0.2*(p(2)-0.5));      % sNu(1); 
materials(1).t = 0.005*(1 + 0.05*(p(3)-0.5));     % sH(1);
materials(1).rho = 7850*(1 + 0.05*(p(4)-0.5));    % sRho(1);

%%% material 2
materials(2).name = 'fluid';
materials(2).type = 'AF_LIN_UAF_ISO_DIR';
materials(2).Id = 2; 
materials(2).cf = 150*(1 + 0.1*(p(5)-0.5));      % sCf(1);
materials(2).rhof = 1.21*(1 + 0.1*(p(6)-0.5));   % sRhof(1);;


%%% load stiffness and mass matrix
[s.KMAT, s.MMAT] = elpasoComputeAndExportSystemMatrices(elpasoControls, elpasoProblemHdf5, materials);
% KMAT and MMAT require postprocessing for fsi problems
[s.KMAT, s.MMAT] = elpasoPostProcessSystemMatricesForFSI(s.KMAT, s.MMAT, s.domains);
[s.KMAT, s.MMAT, s.rdomains] = elpasoRemoveStructInplaneDofsForFSI(s.KMAT, s.MMAT, s.domains);

%%% some function handles
s.getSystemResponseElpaso = @(f) elpasoComputeAndExportSystemSolution(elpasoControls, elpasoProblemHdf5, materials, f);
s.getDynamicSystemMatElpaso = @(f) elpasoComputeAndExportDynamicSystemMatrix(elpasoControls, elpasoProblemHdf5, materials, f);
s.getDynamicSystemMat = @(f) -(2*pi*f)*(2*pi*f)*s.MMAT + s.KMAT;

end