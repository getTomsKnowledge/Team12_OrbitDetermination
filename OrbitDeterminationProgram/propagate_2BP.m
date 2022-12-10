function out = propagate_2BP(t, X, mu)

    % magnitude of position vector, r:
    r = sqrt(X(1,1).^2 + X(2,1).^2 + X(3,1).^2);

    % SHM coefficient(s):
    k = (mu / r^3);

    % rDot in Cartestian coords:
    out(1,1) = X(4,1);
    out(2,1) = X(5,1);
    out(3,1) = X(6,1);

    % rDblDot in Cartesian coords
    out(4,1) = -k * X(1,1);
    out(5,1) = -k * X(2,1);
    out(6,1) = -k * X(3,1);

end