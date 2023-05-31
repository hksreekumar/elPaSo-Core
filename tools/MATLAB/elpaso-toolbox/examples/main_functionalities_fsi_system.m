%% An example matlab script to extract system matrix from elPaSo

clear; close all; clc;
addpath(genpath('..'))

%% Set the parametric system
my_system = @(p) setup_system_fsi(p);

%% Call and export
sys = my_system(zeros(6,1));

%% plot sparsity
figure;
subplot(1,2,1);
spy(sys.KMAT,'k');
title('Stiffness matrix')
hold on;
subplot(1,2,2);
spy(sys.MMAT,'k');
title('Mass matrix')
hold on;

%% plot frf
compute_freq = 100:100:400;
solution = sys.getSystemResponseElpaso(compute_freq);
figure;
plot(compute_freq, solution(end,:), 'LineWidth',2)
hold on;
xlabel('Frequency (Hz)')
ylabel('FRF')

%% plot sparsity
figure;
KDYN = sys.getDynamicSystemMatElpaso(100);
spy(KDYN,'k');
title('Dynamic stiffness matrix')
hold on;