function [visBool, threshold_deg, angle_deg] = r1_r2_visibility(rOne, rTwo)

    visBool = 0; % assume not visible
    rEarth_m = OrbitConstants.R_earth_km*1e3;

    r1 = norm(rOne);
    r2 = norm(rTwo);

    theta1 = acosd(rEarth_m / r1);
    theta2 = acosd(rEarth_m / r2);
    threshold_deg = theta1 + theta2;
    
    % Calculate angle using inner product law:
    angle_deg = acosd(dot(rOne, rTwo) / (r1*r2));

    if (angle_deg <= threshold_deg)
        visBool = 1;
    end

end