% sample_model.m
% Provides sample transmission t(lambda) and phase shift phi(lambda) models.
%
% Usage:
% [T, phi] = sample_model(lambda, sample_params)
% sample_params: struct with fields center_nm, depth, width_nm, phi0

function [T, phi] = sample_model(lambda, sample_params)
    % lambda in meters
    lam_nm = lambda*1e9;
    center = sample_params.center_nm;
    width = sample_params.width_nm;
    depth = sample_params.depth; % amplitude transmission at center (0..1)
    phi0 = sample_params.phi0; % radians center
    % gaussian-shaped transmission profile (amplitude)
    T = depth + (1-depth)*exp(-0.5*((lam_nm - center)/width).^2);
    % phase shift gaussian-like
    phi = phi0 * exp(-0.5*((lam_nm - center)/width).^2);
end
