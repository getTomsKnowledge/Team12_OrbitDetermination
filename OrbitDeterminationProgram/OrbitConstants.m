classdef  OrbitConstants
    properties ( Constant = true )
        % Globals:
        R_earth_km = 6378.1; % km
        mu_earth_km = 3.986004418e5; % earth grav parameter, km^3/s^2
        mu_sun_km = 1.32712440042e11; % sun grav param, km^3/s^2
        fineStruct = 0.007297352525693; % Fine-structure constant
        planck = 6.62607015e-34; % Planck's constant
        eFundamental = 1.602176634e-19; % Fundamental charge
        c = 299792458; % lightspeed m/s
        %{
        mu_0 = 2*fineStruct*planck / ((eFundamental^2)*c);
        epsilon_0 = 1/(mu_0*(c^2));
        %}
    end
end