function [KMAT, MMAT] = elpasoPostProcessSystemMatricesForFSI(KMAT, MMAT, domains)
%ELPASOPOSTPROCESSSYSTEMMATRICESFORFSI Does the postprocessing for KMAT,
%and MMAT
%   
%
% Inputs:
%   - KMAT: Stiffness matrix
%   - MMAT: Mass matrix
%   - domains: [num dofs in structure, num dofs in fluid]
%
% Outputs:
%   - KMAT: Corrected stiffness matrix
%   - MMAT: Corrected mass matrix

% adjustments
ns = domains(1);
nf = domains(2);

% Extract coupling matrix
KCMAT = KMAT(1:ns,ns+1:ns+nf);
KCMAT_T = KMAT(ns+1:ns+nf,1:ns);

% Delete coupling terns from stiffness
KMAT(1:ns,ns+1:nf) = 0;
KMAT(ns+1:nf,1:ns) = 0;

% Add coupling to KMAT and MMAT
KMAT(1:ns,ns+1:ns+nf) = KCMAT;
MMAT(ns+1:ns+nf,1:ns) = -KCMAT';

end

