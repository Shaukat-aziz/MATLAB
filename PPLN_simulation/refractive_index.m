% refractive_index.m
% Sellmeier / refractive index wrapper.
% Usage: n = refractive_index(lambda, Sellmeier_struct)
% lambda in meters.

function n = refractive_index(lambda, Sell)
    % Simple placeholder Sellmeier - assumes coefficients given for microns.
    % Replace with accurate Sellmeier used in your paper (units matter).
    lam_um2 = (lambda*1e6).^2; % lambda in microns^2

    % handle vector input
    A = Sell.A;
    n2 = A * ones(size(lambda));
    if isfield(Sell,'B') && isfield(Sell,'C')
        n2 = n2 + Sell.B./(lam_um2 - Sell.C);
    end
    if isfield(Sell,'D') && isfield(Sell,'E') && Sell.D~=0
        n2 = n2 + Sell.D./(lam_um2 - Sell.E);
    end
    % ensure positive
    n2(n2<=0) = abs(n2(n2<=0)) + 1;
    n = sqrt(n2);
end
