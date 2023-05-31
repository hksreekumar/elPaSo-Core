function [KMAT, MMAT, rdomains] = elpasoRemoveStructInplaneDofsForFSI(KMAT, MMAT, domains)
%elpasoRemoveStructInplaneDofsForFSI Removes the in-place dofs from the
%structural domain
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
%   - rdomains: Reduced domain [num dofs in structure/2, num dofs in fluid]

% Set new domain sizes
rdomains = [domains(1)/2, domains(2)];
plateIndDel = [1:6:domains(1) 2:6:domains(1) 6:6:domains(1)];   % Plate's indices that are to be deleted

KMAT(plateIndDel,:) = [];      % Deletes the respective rows
KMAT(:,plateIndDel) = [];      % Deletes the respective cols

if ~isscalar(MMAT)
    MMAT(plateIndDel,:) = [];      % Deletes the respective rows
    MMAT(:,plateIndDel) = [];      % Deletes the respective cols
end

end

