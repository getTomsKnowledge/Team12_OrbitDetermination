classdef OrbitObject
    properties
        body = "Earth"
        muBody
        % PVT:
        position_m
        rMag
        velocity_ms
        vMag
        epoch
        % Energy:
        hVec
        hMag
        E
        % Perifocal plane, line of nodes:
        nHat
        eHat
        % Keplerian:
        a
        p
        e
        i
        o
        O
        nu
    end
end