% ppln_model_JSI.m
% Simulate PPLN phase-matching and JSI for SPDC (1D along energy-conservation)
% Paste this into MATLAB and run. Replace Sellmeier & poling placeholders below.

clear; close all; clc;

%% ------------------- USER PARAMETERS (replace with PDF values) -------------------
% wavelengths (meters)
lambda_p = 532e-9;       % pump central wavelength (m)
lambda_s_center = 810e-9;% target signal central (m) (lambda_i computed)
% Crystal parameters
L = 20e-3;               % crystal length (m) (replace)
Lambda_poling = 7.5e-6;  % poling period (m) (replace)
d_eff = 2.0e-11;         % effective nonlinear (m/V) placeholder
% Pump
pump_FWHM_nm = 1.0;      % pump FWHM in nm (spectral envelope)
% Temperature
T0 = 25; T = 25;         % reference T0 and simulation T (°C)
dn_dT = 3e-5;            % thermo-optic coef (1/°C) - placeholder
% Pump waist (for angular accept approx)
w0 = 50e-6;              % meters

% Sellmeier placeholder (units: coefficients arranged for microns)
% REPLACE these with the Sellmeier coeffs from your PDF (important!)
Sell.A = 5.35583; Sell.B = 0.100473; Sell.C = 0.20692^2; % placeholder example
Sell.D = 0; Sell.E = 0; Sell.F = 0; Sell.G = 0;

% Numerical settings
nLam = 800;              % resolution for lambda_s sampling
dLam = 8e-9;             % half-window around lambda_s_center to sample (m)
lambda_s_grid = linspace(lambda_s_center - dLam, lambda_s_center + dLam, nLam);

%% ------------------- derived & helper functions -------------------
c = 299792458;
omega = @(lam) 2*pi*c./lam;

% n(lambda,T) using placeholder Sellmeier (assumes Sell coefficients in micron units)
sellmeier_n = @(lam_m) sellmeier_placeholder(lam_m, Sell) + dn_dT*(T - T0);

k_from_n = @(lam, n) 2*pi .* n ./ lam;

% energy conservation: given lambda_s, compute lambda_i
lambda_i_from_s = @(lam_s) 1 ./ (1/lambda_p - 1./lam_s);

% pump spectral envelope (Gaussian in wavelength)
sigma_lambda = (pump_FWHM_nm/2.355)*1e-9;
pump_env = @(lam) exp( - (lam - lambda_p).^2 ./ (2*sigma_lambda^2) );

% QPM grating vector
Kqpm = 2*pi / Lambda_poling;

%% ------------------- compute phase-mismatch + PM amplitude -------------------
Delta_k = zeros(size(lambda_s_grid));
Phi = zeros(size(lambda_s_grid));
for ii = 1:length(lambda_s_grid)
    ls = lambda_s_grid(ii);
    li = lambda_i_from_s(ls);
    % ensure physically meaningful positive lambda_i
    if li <= 0 || imag(li)~=0
        Delta_k(ii) = NaN;
        Phi(ii) = 0;
        continue;
    end
    np = sellmeier_n(lambda_p);
    ns = sellmeier_n(ls);
    ni = sellmeier_n(li);
    kp = k_from_n(lambda_p, np);
    ks = k_from_n(ls, ns);
    ki = k_from_n(li, ni);
    dk = kp - ks - ki - Kqpm;
    Delta_k(ii) = dk;
    % PM amplitude (sinc with argument dk*L/2)
    arg = dk * L / 2;
    Phi(ii) = sinc_safe(arg); % sinc(x) = sin(x)/x
end

% JSA (1D along energy-conserving curve)
JSA = pump_env( (lambda_s_grid + lambda_i_from_s(lambda_s_grid))/2 ) .* Phi;
JSI = abs(JSA).^2;

%% ------------------- angular acceptance (approx) -------------------
% Approx: pump waist -> transverse k spread sigma_k ~ 1/w0; approximate angular spread sigma_theta ~ lambda_p/(pi*w0)
sigma_theta = lambda_p / (pi * w0);
% Angular weighting factor for each pair (assuming paired opposite small angles)
% For coarse visualization assume angular acceptance multiplies JSI by Gaussian in theta mapped from phase mismatch transverse part
theta_grid = linspace(0, 5*sigma_theta, 200); % rad
% We compute a simple angular acceptance gamma(theta) ~ exp(-(k_i*theta*w0)^2/2)
gamma_theta = @(theta, lam) exp( - ( (2*pi*sellmeier_n(lam).*theta./lam).*w0 ).^2 / 2 );

%% ------------------- visualizations -------------------
figure('Units','normalized','Position',[0.05 0.05 0.9 0.8]);

% 1) Delta_k vs lambda_s
subplot(2,3,1);
plot(lambda_s_grid*1e9, Delta_k, '-k','LineWidth',1.2);
xlabel(' \lambda_s (nm)'); ylabel('\Delta k (m^{-1})'); grid on;
title('Phase mismatch \Delta k(\lambda_s)');

% 2) PM amplitude sinc^2 and pump envelope
subplot(2,3,2);
yyaxis left;
plot(lambda_s_grid*1e9, (abs(Phi)).^2, '-r','LineWidth',1.4); ylabel('sinc^2(|\Phi|)');
yyaxis right;
plot(lambda_s_grid*1e9, pump_env((lambda_s_grid + lambda_i_from_s(lambda_s_grid))/2)./max(pump_env((lambda_s_grid + lambda_i_from_s(lambda_s_grid))/2)), '--b');
ylabel('pump env (norm)'); xlabel('\lambda_s (nm)'); grid on; title('Pump × PM');

% 3) JSI (1D)
subplot(2,3,3);
plot(lambda_s_grid*1e9, JSI./max(JSI),'k','LineWidth',1.4);
xlabel('\lambda_s (nm)'); ylabel('JSI (norm)'); grid on; title('1D JSI along energy-cons.');

% 4) JSI heatmap (approx 2D by smearing)
subplot(2,3,4);
LS = repmat(lambda_s_grid*1e9, length(theta_grid),1);
LI = repmat(lambda_i_from_s(lambda_s_grid)*1e9, length(theta_grid),1);
% apply simple angular weighting using idler lam
W = zeros(size(LS));
for j = 1:length(theta_grid)
    W(j,:) = JSI .* gamma_theta(theta_grid(j), lambda_i_from_s(lambda_s_grid));
end
imagesc(lambda_s_grid*1e9, theta_grid, W);
axis xy; xlabel('\lambda_s (nm)'); ylabel('theta (rad)'); title('JSI smeared by angular acceptance');
colorbar;

% 5) Sweep poling period Lambda vs phase-match center (compute lambda_s_peak)
subplot(2,3,5);
Lambda_list = linspace(Lambda_poling*0.8, Lambda_poling*1.2, 60);
peak_ls = zeros(size(Lambda_list));
for j = 1:length(Lambda_list)
    Ktmp = 2*pi / Lambda_list(j);
    dk_temp = zeros(size(lambda_s_grid));
    for ii = 1:length(lambda_s_grid)
        ls = lambda_s_grid(ii); li = lambda_i_from_s(ls);
        np = sellmeier_n(lambda_p); ns = sellmeier_n(ls); ni = sellmeier_n(li);
        dk_temp(ii) = 2*pi*np/lambda_p - 2*pi*ns/ls - 2*pi*ni/li - Ktmp;
    end
    % find lambda where dk crosses zero (closest)
    [~, idxmin] = min(abs(dk_temp));
    peak_ls(j) = lambda_s_grid(idxmin);
end
plot(Lambda_list*1e6, peak_ls*1e9, '-o'); xlabel('Poling period \Lambda (\mum)'); ylabel('phase-match \lambda_s peak (nm)');
grid on; title('Phase-match vs poling period');

% 6) Sweep temperature T effect (if dn/dT nonzero)
subplot(2,3,6);
Tlist = T0 + linspace(-40,40,41);
peak_ls_T = zeros(size(Tlist));
for j = 1:length(Tlist)
    % simple thermo-optic shift by dn/dT
    Sell2 = Sell; % same Sell sheet but we'll adjust via dn_dT in sellmeier_n wrapper
    for ii = 1:length(lambda_s_grid)
        ls = lambda_s_grid(ii); li = lambda_i_from_s(ls);
        np = sellmeier_n(lambda_p) + dn_dT*(Tlist(j)-T0);
        ns = sellmeier_n(ls) + dn_dT*(Tlist(j)-T0);
        ni = sellmeier_n(li) + dn_dT*(Tlist(j)-T0);
        dk_temp(ii) = 2*pi*np/lambda_p - 2*pi*ns/ls - 2*pi*ni/li - Kqpm;
    end
    [~, idxmin] = min(abs(dk_temp));
    peak_ls_T(j) = lambda_s_grid(idxmin);
end
plot(Tlist, peak_ls_T*1e9,'-s'); xlabel('Temperature (°C)'); ylabel('phase-match \lambda_s peak (nm)');
grid on; title('Temperature tuning of phase-match');

%% ------------------- some printed outputs -------------------
% central idler for exact energy cons at desired signal
lambda_i_center = lambda_i_from_s(lambda_s_center);
fprintf('Central idler (for lambda_s = %.1f nm, pump %.1f nm) = %.1f nm\n', lambda_s_center*1e9, lambda_p*1e9, lambda_i_center*1e9);

% FWHM estimate of JSI (in lambda_s)
maxJSI = max(JSI);
half = maxJSI/2;
% find indices around peak
[~, idxpeak] = max(JSI);
left = find(JSI(1:idxpeak) <= half, 1, 'last');
right = idxpeak - 1 + find(JSI(idxpeak:end) <= half, 1, 'first');
if isempty(left); left = 1; end
if isempty(right); right = length(JSI); end
FWHM_nm = (lambda_s_grid(right) - lambda_s_grid(left))*1e9;
fprintf('Approx JSI FWHM (signal) = %.3f nm\n', FWHM_nm);

%% ------------------- helper functions -------------------
function y = sinc_safe(x)
    y = ones(size(x));
    nz = x~=0;
    y(nz) = sin(x(nz))./x(nz);
end

function n = sellmeier_placeholder(lambda_m, SellStruct)
    % Very simple placeholder Sellmeier: assumes coefficients given for microns
    % lam in meters -> convert to microns
    lam_um2 = (lambda_m*1e6).^2;
    A = SellStruct.A;
    if isfield(SellStruct,'B') && isfield(SellStruct,'C')
        B = SellStruct.B; C = SellStruct.C;
        n2 = A + B./(lam_um2 - C);
    else
        n2 = A;
    end
    n = sqrt(abs(n2));
end
