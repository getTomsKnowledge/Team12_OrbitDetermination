function [V1, V2, V3] = gibbs_solver(R1, R2, R3, muBody)

    % magnitudes:
    r1 = norm(R1);
    r2 = norm(R2);
    r3 = norm(R3);

    alphaCoP = asin( ...
        dot(R1, cross(R2, R3))...
        /(r1 * norm(cross(R2, R3))) );

    % Coplanar test:
    if (abs(alphaCoP) > 0.5)
        error('Oops! Non-classical orbit. Vectors are not coplanar.');
        return ;
    end
    
    % cross the p's, dot the i-hats...
    % Cross products:
    iXii = cross(R1,R2);
    iXiii = cross(R1,R3);
    iiXi = cross(R2,R1);
    iiXiii = cross(R2,R3); 
    iiiXi = cross(R3,R1);
    iiiXii = cross(R3,R2);
    % Dot products:
    iOii = dot(R1,R2);
    iOiii = dot(R1,R3);
    iiOi = dot(R2,R1);
    iiOiii = dot(R2,R3);
    
    % Form N(umerator) and D(enominator) vectors:
    N = r3*iXii + r1*iiXiii + r2*iiiXi;      
    D = iXii + iiXiii + iiiXi;
    if ((norm(N) < 0.001) || (norm(D) < 0.001))
        error('Zero sum of scalar triples or cross products.');
        return;
    end
    
    % Semi-latus rectum, p:
    p = norm(N) / norm(D);

    % (p-q-r-)S vector:
    S = (r2 - r3)*R1 + (r3 - r1)*R2 + (r1 - r2)*R3;
   
    % Eccentricity, e:
    e = norm(S) / norm(D);
    
    % Perifocal unit vectors in Cartesian coords:
    qHat = S / norm(S);
    wHat = N / norm(N);
    pHat = cross(qHat, wHat);
    
    % Get B, L terms:
    B1 = cross(D,R1);
    B2 = cross(D,R2);
    B3 = cross(D,R3);
    L = sqrt(muBody / (norm(D) * norm(N)));

    % Calculate velocities:
    V1 = (L/r1) .* B1 + L .* S;
    V2 = (L/r2) .* B2 + L .* S;
    V3 = (L/r3) .* B3 + L .* S;

    % et voila!
end