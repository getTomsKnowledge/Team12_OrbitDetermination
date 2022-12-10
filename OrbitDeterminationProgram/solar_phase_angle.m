function zeta = solar_phase_angle(rDiff, rSun)

    zeta = 0.0; % Assume 0
    r = norm(rDiff);
    rHat = rDiff / r;
    % Assume rSun >> rObj; may use rSun
    zeta = acosd(dot(rHat, rSun));
    
end