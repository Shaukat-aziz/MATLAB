% helper_math.m
% Small helper functions used by modules.
function y = safe_sinc(x)
    y = ones(size(x));
    nz = (x~=0);
    y(nz) = sin(x(nz))./x(nz);
end

function v = trapz2D(X, Y)
    % simple 2D trapezoid flatten
    v = trapz(X, trapz(Y, 1));
end
