% params.m
% Central parameter file for PPLN SPDC simulation (units SI unless noted)
% Edit values here or create params_custom.mat to override.

function P = params()
    % physical constants
    P.c = 299792458;
    P.hbar = 1.054571817e-34;
    P.kB = 1.380649e-23;

    % wavelengths (meters)
    P.lambda_p = 532e-9;       % pump (default)
    P.lambda_s = 810e-9;       % target signal (you can change)
    % idler computed from energy conservation inside code

    % PPLN / nonlinear
    P.L = 20e-3;               % crystal length (m)
    P.Lambda_poling = 7.5e-6;  % poling period (m) - placeholder
    P.d_eff = 2e-11;           % effective nonlinear (m/V) placeholder
    P.w0 = 50e-6;              % pump waist (m)

    % Sellmeier placeholder (replace with real coefficients from PDF)
    P.Sellmeier.A = 5.35583;
    P.Sellmeier.B = 0.100473;
    P.Sellmeier.C = 0.20692^2;
    P.Sellmeier.type = 'placeholder';

    % thermo-optic dn/dT (per degC)
    P.dn_dT = 3e-5;
    P.T0 = 25;                 % reference temp (degC)
    P.T = 25;                  % current temp (degC)

    % numeric grids and ranges
    P.theta_i_max = deg2rad(5);   % max idler angle (rad)
    P.N_theta = 120;
    P.N_omega = 180;
    P.omega_rel_span = 0.03; % fractional bandwidth around central idler frequency

    % idler path length sweep (for interferogram)
    P.dli_min = -5e-3;
    P.dli_max = 5e-3;
    P.N_dli = 200;

    % animation / output controls
    P.doAnimation = true;
    P.animationFrames = 200;
    P.videoFile = fullfile('outputs','PPLN_simulation_demo.mp4');

    % numerical tolerances
    P.verbose = true;

end
