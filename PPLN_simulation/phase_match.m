% phase_match.m
% Compute longitudinal phase mismatch Delta_kz and PM amplitude sinc(Delta_kz*L/2)
% Inputs: lambda_p, lambda_s, lambda_i (m), theta_s, theta_i (rad), P struct
% Outputs: Delta_kz (1/m), PM_amp (complex) where PM_amp = sinc(Delta_kz * L/2)

function [Delta_kz, PM_amp] = phase_match(lambda_p, lambda_s, lambda_i, theta_s, theta_i, P)
    % compute refractive indices (use extraordinary / ordinary as appropriate)
    Sell = P.Sellmeier;
    np = refractive_index(lambda_p, Sell);
    ns = refractive_index(lambda_s, Sell);
    ni = refractive_index(lambda_i, Sell);

    % thermo-optic linear approximation
    dn = P.dn_dT * (P.T - P.T0);
    np = np + dn;
    ns = ns + dn;
    ni = ni + dn;

    % wavevectors
    kp = 2*pi*np ./ lambda_p;
    ks = 2*pi*ns ./ lambda_s;
    ki = 2*pi*ni ./ lambda_i;

    % longitudinal components approx (paraxial): kz = k * cos(theta) ~ k*(1 - theta^2/2)
    kpz = kp .* (1 - 0.5*theta_s.*0); % pump assumed collinear (theta~0)
    ksz = ks .* (1 - 0.5*theta_s.^2);
    kiz = ki .* (1 - 0.5*theta_i.^2);

    % QPM grating vector
    Kqpm = 2*pi / P.Lambda_poling;

    Delta_kz = kpz - ksz - kiz - Kqpm;
    PM_amp = sinc(Delta_kz .* P.L / 2); % MATLAB's sinc(x) = sin(pi*x)/(pi*x); but here we use sin(x)/x, so define below

    % convert to sin(x)/x style (safe)
    PM_amp = sinc_custom(Delta_kz .* P.L ./ 2);
end

function y = sinc_custom(x)
    y = ones(size(x));
    nz = (x~=0);
    y(nz) = sin(x(nz))./x(nz);
end
