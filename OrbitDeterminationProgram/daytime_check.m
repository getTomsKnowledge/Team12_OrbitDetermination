% Takes in ECI vectors for orbiting object and sun (unit):
function result = daytime_check(rObj, rSun)

    result = 0;

    % Take dot product of sun vector and object -- negative is dark side
    dotProduct = dot(rObj, rSun);

    if ( dotProduct >= 0.0)
        result = 1;
    end


end