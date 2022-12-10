function lamOut = lambertSolver(rInit, rFin, TOF, DM, muBody)

    % Initial processing:
    ri = norm(rInit); 
    rf = norm(rFin);
    dNu = acos(dot(rInit, rFin) / ( norm(rInit) * norm(rFin)));
    A = DM * sqrt( (ri * rf) * (1 + cos(dNu)) );

    % Input validation:
    if ((dNu == 0) && (A == 0))
        error('Invalid vector inputs.');
        return;
    end

    % Initial Parameters:
    tolerance = 1e-6;
    psiGuess = 0;
    C2 = 0.5; % Hypergeometric Cos
    C3 = 1.0/6.0; % Hypergeometric Sin
    psiUp = 4 * (pi^2); % Max estimate
    psiLow = -4 * pi; %Min estimate
    dt = 0;
    dPsi = 0.0001;

    while (abs(TOF - dt) > tolerance)
        y = ri + rf + A * ((psiGuess * C3 - 1)/sqrt(C2));

        if ((A > 0) && (y < 0))
            psiLow = psiLow + dPsi;
        end

        x = sqrt(y / C2);
        dt = ((x^3)*C3 + A * sqrt(y))/sqrt(muBody);

        if (dt < TOF)
            psiLow = psiGuess;
        else
            psiUp = psiGuess;
        end

        psiGuess = (psiUp + psiLow)/2.0;

        % Check guess anomaly:
        if (psiGuess > tolerance)
            C2 = (1 - cos(sqrt(psiGuess)))/psiGuess;
            C3 = (sqrt(psiGuess) - sin(sqrt(psiGuess))) / sqrt(psiGuess^3);
        elseif (psiGuess < -tolerance)
                C2 = (1 - cosh(sqrt(-psiGuess)))/psiGuess;
                C3 = (sinh(sqrt(-psiGuess)) - sqrt(-psiGuess))/sqrt(-psiGuess^3);
        else
            C2 = 0.5;
            C3 = 1.0/6.0;
        end 
    end

    % Lagrange coefficients:
    f = 1 - y/ri;
    g = A*sqrt(y / muBody);
    gDot = 1 - y/rf;
    vInit = (rFin - f * rInit)/g;
    vFin = (gDot * rFin - rInit)/g;

    lamOut = [vInit, vFin];
end