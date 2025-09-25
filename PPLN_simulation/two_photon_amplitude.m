% two_photon_amplitude.m
% Compute the (approximate) two-photon amplitude factor along energy-conserving slice.
% Returns amplitude weight W(omega_i, theta_i) including pump envelope and PM function.
%
% Usage:
% W = two_photon_amplitude(P, lambda_s_grid, lambda_i_grid, theta_i_grid)
% or for single points: W = two_photon_amplitude_point(...)

function W = two_photon_amplitude(P, lambda_s_vec, lambda_i_vec, theta_i_vec)
    % pump envelope: Gaussian around P.lambda_p with width implied by P.omega_rel_span
    % We'll use lambda-based Gaussian approx for simplicity
    sigma_frac = P.omega_rel_span; % fractional around center idler
    lam_p = P.lambda_p;
    % create 2D grid (ls x thetai)
    [LS, TI] = ndgrid(lambda_s_vec, theta_i_vec);
    LI = 1./(1/lam_p - 1./LS); % energy conservation
    % pump spectral envelope: use mean of ls and li
    Lmean = (LS + LI)/2;
    sigma_lam = sigma_frac * lam_p;
    pump_env = exp( - (Lmean - lam_p).^2 ./ (2 * sigma_lam^2) );

    % phase-matching amplitude along slice
    PM = zeros(size(pump_env));
    for ii = 1:numel(LS)
        ls = LS(ii); li = LI(ii); ti = TI(ii);
        [~, pm] = phase_match(lam_p, ls, li, 0, ti, P);
        PM(ii) = abs(pm).^2; % weight by sinc^2 (intensity)
    end

    % transverse (pump waist) factor approximate: exp(- (w0^2/2) * (k_i * theta_i)^2 )
    ni_vals = refractive_index(LI(:), P.Sellmeier);
    ki_vals = 2*pi*ni_vals ./ LI(:);
    KI = reshape(ki_vals, size(PM));
    trans_factor = exp( - (P.w0.^2 / 2) .* (KI .* TI).^2 );

    % final amplitude weight (unnormalized)
    W = pump_env .* PM .* trans_factor;
end
