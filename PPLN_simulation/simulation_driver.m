% simulation_driver.m
% Top-level runner: loads params, computes R(dli), shows plots and optionally saves results.

clear; close all; clc;
if ~exist('outputs','dir'); mkdir('outputs'); end

P = params(); % load parameter struct

% compute idler central wavelength (energy conservation)
P.lambda_i = 1/(1/P.lambda_p - 1/P.lambda_s);
fprintf('Computed central idler lambda = %.2f nm\n', P.lambda_i*1e9);

% integrate R(dli)
[dli_vec, Rvec] = integrate_R0(P);

% save results
save(fullfile('outputs','R_vs_dli.mat'),'dli_vec','Rvec','P');

% plot
figure('Name','R vs dli','NumberTitle','off');
plot(dli_vec*1e3, Rvec, '-','LineWidth',1.6);
xlabel('\Delta l_i (mm)'); ylabel('R (normalized)');
title('Simulated interferogram R(\Delta l_i)');
grid on;

% animate beams & overlay R if requested
if P.doAnimation
    animate_beams(P, dli_vec, Rvec);
end

fprintf('simulation_driver: done. Results saved to outputs/. \n');
