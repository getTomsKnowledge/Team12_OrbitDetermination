%% Line of Sight Calculator
% Takes right ascensions, declinations; spits out Cartesian line-of-sight
% vectors

function [los1, los2, los3] = get_los(ra1, dec1, ra2, dec2, ra3, dec3)

    angles = zeros(3,2);
    angles = [ra1, dec1;
              ra2, dec2;
              ra3, dec3];
    
    cr1 = cosd(angles(1,1));
    cd1 = cosd(angles(1,2));
    sr1 = sind(angles(1,1));
    sd1 = sind(angles(1,2));
    cr2 = cosd(angles(2,1));
    cd2 = cosd(angles(2,2));
    sr2 = sind(angles(2,1));
    sd2 = sind(angles(2,2));
    cr3 = cosd(angles(3,1));
    cd3 = cosd(angles(3,2));
    sr3 = sind(angles(3,1));
    sd3 = sind(angles(3,2));
    
    los1 = [cd1*cr1;
            cd1*sr1;
            sd1];
    los2 = [cd2*cr2;
            cd2*sr2;
            sd2];
    los3 = [cd3*cr3;
            cd3*sr3;
            sd3];
    
end