%% An example matlab script to extract system matrix from elPaSo

clear; close all; clc;
addpath(genpath('..'))

%% Set the parametric system
my_system = @(p) test_setup_system_fsi(p);

%% Call and export
sys = my_system(0.5*ones(6,1));

%% dof test
freq = 100;
KDYN_100Hz = sys.getDynamicSystemMat(freq);
sol_100Hz = sys.getSystemResponseElpaso(freq);

test_load = zeros(length(KDYN_100Hz),1);
test_load(157,1) = 1;

a = KDYN_100Hz\test_load;
a(1)
sol_100Hz(3,1)
display('Should be -1.2961e-06')
a(sys.rdomains(1)+1)
sol_100Hz(sys.domains(1)+1,1)
display('Should be 0.1380')

%% frf test
freq = 100:10:400;

% compute in matlab
sol_frf_matlab = [];
for ifreq = 1:length(freq)
    sol_curr_freq = sys.getDynamicSystemMat(freq(ifreq))\test_load;
    sol_frf_matlab = [sol_frf_matlab, sol_curr_freq];
end
% compute with elpaso
sol_frf_elpaso = sys.getSystemResponseElpaso(freq);
% reference
load('reference_frf.mat');

figure;
semilogy(freq, abs(sol_frf(sys.domains(1)+1,:)),'-','LineWidth',2);
hold on;
semilogy(freq, abs(sol_frf_elpaso(sys.domains(1)+1,:)),'--','LineWidth',2);
hold on;
semilogy(freq, abs(sol_frf_matlab(sys.rdomains(1)+1,:)),'-.','LineWidth',2);
hold on;