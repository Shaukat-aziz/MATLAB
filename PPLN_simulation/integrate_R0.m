% integrate_R0.m
% Numerically integrates the SPDC counts R(dli) over omega_i and theta_i
% Implements a simplified version of the integral in the paper.
%
% Usage:
% [dli_vec, Rvec] = integrate_R0(P)
% returns Rvec normalized (arb. units) across dli sweep

function [dli_vec, Rvec] = integrate_R0(P)
    % compute idler central wavelength from energy conservation
    lambda_p = P.lambda_p; lambda_s = P.lambda_s;
    lambda_i0 = 1/(1/lambda_p - 1/lambda_s);

    % build omega (or lambda) and theta grids
    theta_vec = linspace(0, P.theta_i_max, P.N_theta);
    lam_i_span = lambda_i0 * (1 + P.omega_rel_span*[-1 1]);
    lam_i_vec = linspace(lam_i_span(1), lam_i_span(2), P.N_omega);

    % dli sweep
    dli_vec = linspace(P.dli_min, P.dli_max, P.N_dli);
    Rvec = zeros(size(dli_vec));

    % Precompute pump and PM along grid quickly -> build corresponding lambda_s per lam_i
    % energy conservation: 1/lam_p = 1/lam_s + 1/lam_i => lam_s = 1/(1/lam_p - 1/lam_i)
    lam_p = lambda_p;
    [LI, TH] = meshgrid(lam_i_vec, theta_vec); % TH rows = theta, columns lam_i
    LS = 1./(1/lam_p - 1./LI);

    % compute PM amplitude and pump envelope and trans factor
    W = two_photon_amplitude(P, LS(:)', LI(:)', TH(1,:) ); % this call is generic; to keep dims consistent we'll recompute below

    % Instead compute manually to avoid mismatch: loop reasonably sized grids
    Wmat = zeros(length(theta_vec), length(lam_i_vec));
    for jj = 1:length(lam_i_vec)
        li = lam_i_vec(jj);
        ls = 1/(1/lam_p - 1/li);
        for ii = 1:length(theta_vec)
            ti = theta_vec(ii);
            % pump envelope at mean lambda (ls+li)/2
            sigma_lam = P.omega_rel_span * lam_p;
            pump_env = exp(-(( (ls+li)/2 - lam_p).^2)/(2*sigma_lam^2));
            % PM amplitude
            [Delta_kz, PM_amp] = phase_match(lam_p, ls, li, 0, ti, P);
            PMsq = abs(PM_amp).^2;
            % transverse factor
            ni = refractive_index(li, P.Sellmeier);
            ki = 2*pi*ni/ li;
            trans = exp( - (P.w0.^2/2) * (ki * ti)^2 );
            Wmat(ii,jj) = pump_env .* PMsq .* trans;
        end
    end

    % numeric quadrature weights (simple trapezoid)
    dtheta = theta_vec(2)-theta_vec(1);
    dlam = lam_i_vec(2)-lam_i_vec(1);
    for k = 1:length(dli_vec)
        dli = dli_vec(k);
        % interference term cos(phi0 + omega_i/c * dli). We'll set phi0=0 for simplicity.
        omega_i = 2*pi * P.c ./ LI; % but LI size mismatch; rebuild omega_i grid
        omega_grid = 2*pi * P.c ./ repmat(lam_i_vec, length(theta_vec), 1);
        costerm = cos(omega_grid * dli / P.c); % dimensionless
        integrand = Wmat .* (1 + costerm); % 1 + cos(omega/c * dli)
        Rvec(k) = trapz(lam_i_vec, trapz(theta_vec, integrand, 1), 2);
    end

    % normalize to max 1 for convenience
    Rvec = real(Rvec);
    if max(Rvec)>0
        Rvec = Rvec ./ max(Rvec);
    end

    if P.verbose
        fprintf('integrate_R0: computed R(dli) on %d x %d grid.\n', length(theta_vec), length(lam_i_vec));
    end
end
