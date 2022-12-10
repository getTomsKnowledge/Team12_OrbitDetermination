function [Xf, iterations] = trilateration_solver(...
    r1, rho1, r2, rho2, r3, rho3, tolerance)

    % Check if inputs row or column vectors:
    if (iscolumn(r1))
        r1 = r1';
    end
    if (iscolumn(r2))
        r2 = r2';
    end
    if (iscolumn(r3))
        r3 = r3';
    end
    
    % Prepare reference matrices:
    Xi = [r1; r2; r3]; % 3 X 3 Matrix of reference points
    rhoObs = [rho1, rho2, rho3]';

    % Algorithm:

    % 1.)  Start w/ estimation. x_0 = origin, always converges:
    Xn = [0,0,0]';
    diffA = zeros(3,3);
    rhoPred = zeros(3,1);
    A = zeros(3,3);
    counter = 0;
    epsilon = zeros(3,1);

    while (1)

        % Increment iteration counter:
        counter = counter + 1;

        % Get the variance norms:
        for n = 1:3
            for m = 1:3
                diffA(n,m) = -(Xi(n,m) - Xn(m));
            end
        end

        % Matrix of prediction norms:
        rhoPred = [norm(diffA(1,1:3)), ...
            norm(diffA(2,1:3)), ...
            norm(diffA(3,1:3))]';

        % 2.) Compute epsilon:
        for n = 1:length(rhoPred)
            epsilon(n) = rhoObs(n) - rhoPred(n); % compute the residual
        end

        % 3.) If the residual magnitude is small enough, quit w/ estimate:

        residual = norm(epsilon);

        if (residual < tolerance)

            Xf = Xn;
            iterations = counter;

            break;

        end
        
        % 4.)  Compute the Jacobian A: 

        % Jacobian of the predicted range:
        A = [ diffA(1,1)/rhoPred(1), diffA(1,2)/rhoPred(1), diffA(1,3)/rhoPred(1);
              diffA(2,1)/rhoPred(2), diffA(2,2)/rhoPred(2), diffA(2,3)/rhoPred(2);
              diffA(3,1)/rhoPred(3), diffA(3,2)/rhoPred(3), diffA(3,3)/rhoPred(3)];

        % 5.) Compute the Newtonian iteration step:
        Xn = Xn + inv(A)*epsilon;

        % 6.) Return to step 2.
    end

end