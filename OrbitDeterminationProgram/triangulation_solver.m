function [rGuess, epsilon] = triangulation_solver(b1, c1, b2, c2)

    % Test for column inputs:
    if (isrow(c1))
        c1 = c1';
    end
    if (isrow(c2))
        c2 = c2';
    end
    % Set up linear system, solve:
    d = b1 - b2;
    A = [-(norm(c1)^2) dot(c1, c2); -dot(c1, c2) norm(c2)^2];
    D = [dot(c1, d), dot(c2, d)]';
    S = inv(A)*D;
    s1 = S(1);
    s2 = S(2);
    
    rGuess_1 = b1 + s1*c1;
    rGuess_2 = b2 + s2*c2;
    diff = rGuess_1 - rGuess_2;
    epsilon = norm(diff);
    
    midpoint_eci = 0.5*(rGuess_2 - rGuess_1);
    rGuess = rGuess_1 + midpoint_eci;


end