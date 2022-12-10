function orbit = element_solver(R, V, muBody)

    % Make column vectors:
    if (isrow(R))
        R = R';
    end
    if (isrow(V))
        V = V';
    end

    % Create object for element storage:
    orbit = OrbitObject();

    % mu:
    orbit.muBody = muBody;

    %PV:
    orbit.position_m = zeros(3,1);
    orbit.position_m = R;
    orbit.velocity_ms = zeros(3,1);
    orbit.velocity_ms = V;

    % Magnitudes, r & v:
    orbit.rMag = norm(R);
    orbit.vMag = norm(V);

    % Angular momentum, H & h:
    orbit.hVec = zeros(3,1);
    orbit.hVec = cross(R,V);
    orbit.hMag = norm(orbit.hVec);

    % Line of nodes, n:
    kHat = [0 0 1]';
    N = cross(kHat, orbit.hVec);
    orbit.nHat = N ./ norm(N);

    % Eccentricity, eVec & e:
    eVec = (1 / muBody)*(((orbit.vMag^2) - muBody/orbit.rMag) * R - (dot(R,V)) * V);
    orbit.eHat = eVec / norm(eVec);
    orbit.e = norm(eVec);

    % Semi-latus rectum, p:
    orbit.p = (orbit.hMag^2)/muBody;

    % Semi-major axis, a:
    orbit.a = orbit.p / (1 - orbit.e^2);

    % Inclination, i:
    orbit.i = acos(orbit.hVec(3,1)/orbit.hMag);

end