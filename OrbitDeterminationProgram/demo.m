%% HOUSEKEEPING

%{
run('TidyWorkspace.m');
run('CheckSNaG.m');
%}

close all;
clear all;
home;
clc;


printing = PrintFormatting();
fprintf('Group 12 Orbit Determination Program\nAuthor: Tom West\n\n');

%{
heading = HeadingGenerator('ENAE441', 'Orbit Determination Project', 'Determine GEO satellite orbit based on observation data');
heading.print_heading;
%}

consts = OrbitConstants();

%% DATA SAMPLING/SORTING

% Select from multiple group data sets for comparison:
chosenSet = input("Type '12' for group 12, '1' for Group 1, '11' for Group 11.");
%chosenSet = str2num(groupInput);
switch chosenSet
    case 1
        selectedGroupN1 = 'opt2satAset1.mat';
        selectedGroupN2 = 'opt2satAset1.mat';
    case 11
        selectedGroupN1 = 'opt2satCset3.mat';
        selectedGroupN2 = 'opt2satCset3.mat';
    case 12
        selectedGroupN1 = 'opt2satCset4.mat';
        selectedGroupN2 = 'opt3satCset4.mat';
    otherwise
        fprintf('INVALID ENTRY!');
        return;
end

% Load user-selected dataset:
night1Data = load(selectedGroupN1);
night2Data = load(selectedGroupN2);


waitBar1 = waitbar(0, 'Initializing orbit determination. Please wait...');

%fprintf('\n\n Beginning Sample:\n\n');

% PREPARE RANDOM SAMPLE:

% Get random sample of size N from observation set:
totalOrbits = length(night1Data.opt2satCset4.datetime);
%totalOrbits = length(opt2satCset4.datetime);
sampleSize = 12;
numberOfGroups = 3;
K = 6;

outOf = factorial(sampleSize) / (factorial(sampleSize - 3)*factorial(3));

azErr = zeros(3*outOf, 1);
elErr = zeros(3*outOf, 1);

run('get_random_od.m');

[beginningSample, endingSample] = random_single_observation_sample(totalOrbits, sampleSize);
begSampleList = beginningSample;
endSampleList = endingSample;
%sampleVector = zeros(sampleSize, K);


%% BEGINNING GROUP:

for n = 1:sampleSize
    sampleVector.obsNumber(n) = satC.observationNumbers(begSampleList(n));
    sampleVector.datetimes(n) = satC.datetimeList(begSampleList(n));
    sampleVector.azimuths(n) = satC.azimuths(begSampleList(n));
    sampleVector.elevations(n) = satC.elevations(begSampleList(n));
    sampleVector.ras(n) = satC.ras(begSampleList(n));
    sampleVector.declinations(n) = satC.declinations(begSampleList(n));
end

begCounter = 0;

% Comparative accuracy of sample combinations:
% START TIMER:

waitX = 0.3;
waitbar(waitX, waitBar1, 'Processing beginning observation set.');
dX = 0.001;

minRMS = 1*1e8;
tic

for n = 1:(sampleSize-2)    
    for k = 1:(sampleSize - 1)
        for l = 1:sampleSize

            if ((k >= (n + 1)) && ((l >= (n + 2)) && (l >= (k + 1))))

                begCounter = begCounter + 1;

                waitMsg = sprintf('Processing beginning observation set. Count: %d / %d', begCounter, outOf);
                waitbar(waitX, waitBar1, waitMsg);
                waitX = waitX + dX;

                % Get nth observation set:

                % First observation, n / "enn":
                obs1.obsNumber = sampleVector.obsNumber(n);
                obs1.epoch = sampleVector.datetimes(n);
                obs1.azimuth_deg = sampleVector.azimuths(n);
                obs1.elevation_deg = sampleVector.elevations(n);
                obs1.range_m = 5000;
                obs1.right_ascension_deg = sampleVector.ras(n);
                obs1.declination_deg = sampleVector.declinations(n);
                chile1.position_m = chileObservatory.position_m;
                chile1.epoch = obs1.epoch;                
                chile1.latitude_deg = chileLLA.latitude_deg;
                chile1.longitude_deg = chileLLA.longitude_deg;
                chile1.altitude_m = chileLLA.altitude_m;
                %obs1ECI = eci(obs1, chile1);
                %obs1ECI.position_m = obs1ECI.position_m / norm(obs1ECI.position_m);
                
                t1 = second(chile1.epoch, 'secondofday');


                % Second observation, k / "kay":
                obs2.obsNumber = sampleVector.obsNumber(k);
                obs2.epoch = sampleVector.datetimes(k);
                obs2.azimuth_deg = sampleVector.azimuths(k);
                obs2.elevation_deg = sampleVector.elevations(k);
                obs2.range_m = 5000;
                obs2.right_ascension_deg = sampleVector.ras(k);
                obs2.declination_deg = sampleVector.declinations(k);
                chile2.position_m = chileObservatory.position_m;
                chile2.epoch = obs2.epoch;
                chile2.latitude_deg = chileLLA.latitude_deg;
                chile2.longitude_deg = chileLLA.longitude_deg;
                chile2.altitude_m = chileLLA.altitude_m;
                %obs2ECI = eci(obs2, chile2);
                %obs2ECI.position_m = obs2ECI.position_m / norm(obs2ECI.position_m);

                t2 = second(chile2.epoch, 'secondofday');

                % Third observation, l / "ell":
                obs3.obsNumber = sampleVector.obsNumber(l);
                obs3.epoch = sampleVector.datetimes(l);
                obs3.azimuth_deg = sampleVector.azimuths(l);
                obs3.elevation_deg = sampleVector.elevations(l);
                obs3.range_m = 5000;
                obs3.right_ascension_deg = sampleVector.ras(l);
                obs3.declination_deg = sampleVector.declinations(l);
                chile3.position_m = chileObservatory.position_m;
                chile3.epoch = obs3.epoch;
                chile3.latitude_deg = chileLLA.latitude_deg;
                chile3.longitude_deg = chileLLA.longitude_deg;
                chile3.altitude_m = chileLLA.altitude_m;
                %obs3ECI = eci(obs3, chile3);
                %obs3ECI.position_m = obs3ECI.position_m / norm(obs3ECI.position_m);

                t3 = second(chile3.epoch, 'secondofday');

                % Get line of sight vectors:
                [los1, los2, los3] = get_los(...
                    obs1.right_ascension_deg, obs1.declination_deg, ...
                    obs2.right_ascension_deg, obs2.declination_deg, ...
                    obs3.right_ascension_deg, obs3.declination_deg);

                [o1out.position_m, o2out.position_m, o3out.position_m] = gauss_solver(...
                    los1, los2, los3, ...
                    chile1.position_m, chile2.position_m, chile3.position_m, ...
                    t1, t2, t3);

                o1out.epoch = chile1.epoch;
                o2out.epoch = chile2.epoch;
                o3out.epoch = chile3.epoch;
                
                chileLLA.epoch = o1out.epoch;
                o1out.aer = aer(o1out, chile1);
                chileLLA.epoch = o2out.epoch;
                o2out.aer = aer(o2out, chile2);
                chileLLA.epoch = o3out.epoch;
                o3out.aer = aer(o3out, chile3);

                az1res = obs1.azimuth_deg - o1out.aer.azimuth_deg;
                el1res = obs1.elevation_deg - o1out.aer.elevation_deg;

                az2res = obs2.azimuth_deg - o2out.aer.azimuth_deg;
                el2res = obs2.elevation_deg - o2out.aer.elevation_deg;

                az3res = obs3.azimuth_deg - o3out.aer.azimuth_deg;
                el3res = obs3.elevation_deg - o3out.aer.elevation_deg;

                RMS = (az1res)^2 + (el1res)^2; 
                RMS = (az2res)^2 + (el2res)^2;
                RMS = (az3res)^2 + (el3res)^2;
                RMS = sqrt( RMS/3 );

                azErr(begCounter+mixCounter) = az1res;
                azErr(begCounter+mixCounter+1) = az2res;
                azErr(begCounter+mixCounter+2) = az3res;
                elErr(begCounter+mixCounter) = el1res;
                elErr(begCounter+mixCounter+1) = el2res;
                elErr(begCounter+mixCounter+2) = el3res;

                % Save RMS values for distribution analysis (histogram?):
                begRMSDistribution(begCounter) = RMS;

                if (RMS < minRMS)
                    minRMS = RMS;
                    index1 = n;
                    %{
                    R1 = o1out.position_m;
                    R2 = o2out.position_m;
                    R3 = o3out.position_m;
                    %}

                    %{
                    [V1, V2, V3] = gibbs_solver(...
                        R1,R2,R3,consts.mu_earth_km*1e9);
                    %}
                    %R1_final = R1;
                    index2 = k;
                    %R2_final = R2;
                    index3 = l;
                    %R3_final = R3;
                end

%                fprintf('count: %d, obs nums: %d , %d , %d, residual: %0.3f  m\n',...
%                   begCounter, obs1.obsNumber, obs2.obsNumber, obs3.obsNumber, RMS);


%                [V1, V2, V3] = gibbs_solver(...
%                    R1, R2, R3, consts.mu_earth_km*1e9);

            end
%{
            % get residual RMS value for nth dataset
            nthResidual = 0;
        
            for m = 1:3
                observedAZ;
                predAZ;
                observedEL;
                predEL;
                nthResidual = nthResidual + sqrt((1/N)*((- nthAZ)^2 + (- nthEL)^2));
            end
        
            % Check if better solution exists:
            if (nthResidual < bestResidual)
                bestResidual = nthResidual;
                bestObservation = nthObservation;
            end
                %}
        end
    end    
end

% STOP TIMER:
timeOut = toc;

% Final observation triplet:

% First observation, n / "enn":
final1.obsNumber = sampleVector.obsNumber(index1);
final1.epoch = sampleVector.datetimes(index1);
final1.azimuth_deg = sampleVector.azimuths(index1);
final1.elevation_deg = sampleVector.elevations(index1);
final1.range_m = 5000;
final1.right_ascension_deg = sampleVector.ras(index1);
final1.declination_deg = sampleVector.declinations(index1);
chile1.position_m = chileObservatory.position_m;
chile1.epoch = final1.epoch;
chile1.latitude_deg = chileLLA.latitude_deg;
chile1.longitude_deg = chileLLA.longitude_deg;
chile1.altitude_m = chileLLA.altitude_m;
%final1ECI = eci(final1, chile1);
%final1ECI.position_m = final1ECI.position_m / norm(final1ECI.position_m);

% Second observation, k / "kay":
final2.obsNumber = sampleVector.obsNumber(index2);
final2.epoch = sampleVector.datetimes(index2);
final2.azimuth_deg = sampleVector.azimuths(index2);
final2.elevation_deg = sampleVector.elevations(index2);
final2.right_ascension_deg = sampleVector.ras(index2);
final2.declination_deg = sampleVector.declinations(index2);
final2.range_m = 5000;
chile2.position_m = chileObservatory.position_m;
chile2.epoch = final2.epoch;
chile2.latitude_deg = chileLLA.latitude_deg;
chile2.longitude_deg = chileLLA.longitude_deg;
chile2.altitude_m = chileLLA.altitude_m;
%final2ECI = eci(final2, chile2);
%final2ECI.position_m = final2ECI.position_m / norm(final2ECI.position_m);

% Third observation, l / "ell":
final3.obsNumber = sampleVector.obsNumber(index3);
final3.epoch = sampleVector.datetimes(index3);
final3.azimuth_deg = sampleVector.azimuths(index3);
final3.elevation_deg = sampleVector.elevations(index3);
final3.right_ascension_deg = sampleVector.ras(index3);
final3.declination_deg = sampleVector.declinations(index3);
final3.range_m = 5000;
chile3.position_m = chileObservatory.position_m;
chile3.epoch = final3.epoch;
chile3.latitude_deg = chileLLA.latitude_deg;
chile3.longitude_deg = chileLLA.longitude_deg;
chile3.altitude_m = chileLLA.altitude_m;
%final3ECI = eci(final3, chile3);
%final3ECI.position_m = final3ECI.position_m / norm(final3ECI.position_m);

t1 = second(final1.epoch, 'secondofday');
t2 = second(final2.epoch, 'secondofday');
t3 = second(final3.epoch, 'secondofday');

% Get line of sight vectors:
[finalLOS1, finalLOS2, finalLOS3] = get_los(...
    final1.right_ascension_deg, final1.declination_deg, ...
    final2.right_ascension_deg, final2.declination_deg, ...
    final3.right_ascension_deg, final3.declination_deg);

[R1FinalBeginning.position_m, R2FinalBeginning.position_m, R3FinalBeginning.position_m] = gauss_solver(...
    finalLOS1, finalLOS2, finalLOS3, ...
    chile1.position_m, chile2.position_m, chile3.position_m, ...
    t1, t2, t3);

beginningGraphTitle = "Beginning: good data";

try
    [V1FinalBeginning, V2FinalBeginning, V3FinalBeginning] = gibbs_solver(...
        R1FinalBeginning.position_m, R2FinalBeginning.position_m, R3FinalBeginning.position_m, ...
        consts.mu_earth_km*1e9);

%    fprintf('\n\nnope, bad beginning dataset.\n\n');
catch ME

    warning('Bad position vectors in Beginning data set.  Replacing w/ dummy set.');
    beginningGraphTitle = "Beginning: BAD data";

    % Load representative data set:
    goodData = load('goodObservations.mat');

    % Good observation triplet:
    R1FinalBeginning.position_m = goodData.goodPositions(:,1);
    R2FinalBeginning.position_m = goodData.goodPositions(:,2);
    R3FinalBeginning.position_m = goodData.goodPositions(:,3);
    V1FinalBeginning = goodData.goodVelocities(:,1);
    V2FinalBeginning = goodData.goodVelocities(:,2);
    V3FinalBeginning = goodData.goodVelocities(:,3);
end

% Create histogram of data:
figure('Name', 'Histogram of Beginning RMS Distributions');
histogram(begRMSDistribution);
title('Histogram of Beginning RMS Distributions - Team 12, ENAE441');

xIn1 = [R1FinalBeginning.position_m; V1FinalBeginning];
xIn2 = [R2FinalBeginning.position_m; V2FinalBeginning];
xIn3 = [R3FinalBeginning.position_m; V3FinalBeginning];
tRange = [0 12*60*60];

[tOut, beginningOrbit1] = ode45(@propagate_2BP, tRange, xIn1, [], consts.mu_earth_km*1e9);
[tOut, beginningOrbit2] = ode45(@propagate_2BP, tRange, xIn2, [], consts.mu_earth_km*1e9);
[tOut, beginningOrbit3] = ode45(@propagate_2BP, tRange, xIn3, [], consts.mu_earth_km*1e9);

%{
figure('Name', 'BEGINNING Orbit Plot');
plot3(beginningOrbit1(:,1), beginningOrbit1(:,2), beginningOrbit1(:,3), 'Color', 'r');
hold on;
plot3(beginningOrbit2(:,1), beginningOrbit2(:,2), beginningOrbit2(:,3), 'Color', 'g');
hold on;
plot3(beginningOrbit3(:,1), beginningOrbit3(:,2), beginningOrbit3(:,3), 'Color', 'b');
title(beginningGraphTitle);
%}


fprintf('\n%s', beginningGraphTitle);
fprintf('\n\nBest Beginning Triplet:  ( %d, %d, %d ) [obs. no.],  Min. RMS: %0.4f deg.\n', ...
    final1.obsNumber, final2.obsNumber, final3.obsNumber, minRMS);

fprintf('Sample Size: %d , Time Elapsed: %0.1f seconds\n\n', ...
    sampleSize, timeOut);


%% ENDING DATASET:


for n = 1:sampleSize
    sampleVector.obsNumber(n)    = satC.observationNumbers(endSampleList(n));
    sampleVector.datetimes(n)    = satC.datetimeList(endSampleList(n));
    sampleVector.azimuths(n)     = satC.azimuths(endSampleList(n));
    sampleVector.elevations(n)   = satC.elevations(endSampleList(n));
    sampleVector.ras(n)          = satC.ras(endSampleList(n));
    sampleVector.declinations(n) = satC.declinations(endSampleList(n));
end

endCounter = 0;

% Comparative accuracy of sample combinations:
% START TIMER:

waitX = 0.6;
waitbar(waitX, waitBar1, 'Processing ending observation set.');
dX = 0.001;

minRMS = 1*1e8;
tic

for n = 1:(sampleSize-2)    
    for k = 1:(sampleSize - 1)
        for l = 1:sampleSize

            if ((k >= (n + 1)) && ((l >= (n + 2)) && (l >= (k + 1))))

                endCounter = endCounter + 1;
                waitMsg = sprintf('Processing ending observation set. Count: %d / %d', endCounter, outOf); 
                waitbar(waitX, waitBar1, waitMsg);
                % Get nth observation set:


                % First observation, n / "enn":
                obs1.obsNumber = sampleVector.obsNumber(n);
                obs1.epoch = sampleVector.datetimes(n);
                obs1.right_ascension_deg = sampleVector.ras(n);
                obs1.declination_deg = sampleVector.declinations(n);
                obs1.azimuth_deg = sampleVector.azimuths(n);
                obs1.elevation_deg = sampleVector.elevations(n);
                %obs1.range_m = 5000;
                chile1.position_m = chileObservatory.position_m;
                chile1.epoch = obs1.epoch;
                chile1.latitude_deg = chileLLA.latitude_deg;
                chile1.longitude_deg = chileLLA.longitude_deg;
                chile1.altitude_m = chileLLA.altitude_m;
                %obs1ECI = eci(obs1, chile1);
                %obs1ECI.position_m = obs1ECI.position_m / norm(obs1ECI.position_m);
                
                t1 = second(chile1.epoch, 'secondofday');

                % Second observation, k / "kay":
                obs2.obsNumber = sampleVector.obsNumber(k);
                obs2.epoch = sampleVector.datetimes(k);
                obs2.azimuth_deg = sampleVector.azimuths(k);
                obs2.elevation_deg = sampleVector.elevations(k);
                %obs2.range_m = 5000;
                obs2.right_ascension_deg = sampleVector.ras(k);
                obs2.declination_deg = sampleVector.declinations(k);
                chile2.position_m = chileObservatory.position_m;
                chile2.epoch = obs2.epoch;
                chile2.latitude_deg = chileLLA.latitude_deg;
                chile2.longitude_deg = chileLLA.longitude_deg;
                chile2.altitude_m = chileLLA.altitude_m;
                %obs2ECI = eci(obs2, chile2);
                %obs2ECI.position_m = obs2ECI.position_m / norm(obs2ECI.position_m);

                t2 = second(chile2.epoch, 'secondofday');

                % Third observation, l / "ell":
                obs3.obsNumber = sampleVector.obsNumber(l);
                obs3.epoch = sampleVector.datetimes(l);
                obs3.azimuth_deg = sampleVector.azimuths(l);
                obs3.elevation_deg = sampleVector.elevations(l);
                obs3.range_m = 5000;
                obs3.right_ascension_deg = sampleVector.ras(l);
                obs3.declination_deg = sampleVector.declinations(l);
                chile3.position_m = chileObservatory.position_m;
                chile3.epoch = obs3.epoch;
                chile3.latitude_deg = chileLLA.latitude_deg;
                chile3.longitude_deg = chileLLA.longitude_deg;
                chile3.altitude_m = chileLLA.altitude_m;
               % obs3ECI = eci(obs3, chile3);
               % obs3ECI.position_m = obs3ECI.position_m / norm(obs3ECI.position_m);

                t3 = second(chile3.epoch, 'secondofday');


                % Get line of sight vectors:
                [los1, los2, los3] = get_los(...
                    obs1.right_ascension_deg, obs1.declination_deg, ...
                    obs2.right_ascension_deg, obs2.declination_deg, ...
                    obs3.right_ascension_deg, obs3.declination_deg);

                [o1out.position_m, o2out.position_m, o3out.position_m] = gauss_solver(...
                    los1, los2, los3, ...
                    chile1.position_m, chile2.position_m, chile3.position_m, ...
                    t1, t2, t3);

                o1out.epoch = chile1.epoch;
                o2out.epoch = chile2.epoch;
                o3out.epoch = chile3.epoch;
                
                chileLLA.epoch = o1out.epoch;
                o1out.aer = aer(o1out, chile1);
                chileLLA.epoch = o2out.epoch;
                o2out.aer = aer(o2out, chile2);
                chileLLA.epoch = o3out.epoch;
                o3out.aer = aer(o3out, chile3);

                az1res = obs1.azimuth_deg - o1out.aer.azimuth_deg;
                el1res = obs1.elevation_deg - o1out.aer.elevation_deg;

                az2res = obs2.azimuth_deg - o2out.aer.azimuth_deg;
                el2res = obs2.elevation_deg - o2out.aer.elevation_deg;

                az3res = obs3.azimuth_deg - o3out.aer.azimuth_deg;
                el3res = obs3.elevation_deg - o3out.aer.elevation_deg;

                RMS = (az1res)^2 + (el1res)^2; 
                RMS = (az2res)^2 + (el2res)^2;
                RMS = (az3res)^2 + (el3res)^2;
                RMS = sqrt( RMS/3 );

                azErr(begCounter + mixCounter + endCounter) = az1res;
                azErr(begCounter + mixCounter + endCounter + 1) = az2res;
                azErr(begCounter + mixCounter + endCounter + 2) = az3res;
                elErr(begCounter + mixCounter + endCounter) = el1res;
                elErr(begCounter + mixCounter + endCounter + 1) = el2res;
                elErr(begCounter + mixCounter + endCounter + 2) = el3res;

                % Save RMS values for distribution analysis (histogram?):
                endRMSDistribution(endCounter) = RMS;

                if (RMS < minRMS)
                    minRMS = RMS;
                    index1 = n;
                    %{
                    R1 = o1out.position_m;
                    R2 = o2out.position_m;
                    R3 = o3out.position_m;
                    %}

                    %{
                    [V1, V2, V3] = gibbs_solver(...
                        R1,R2,R3,consts.mu_earth_km*1e9);
                    %}
                    %R1_final = R1;
                    index2 = k;
                    %R2_final = R2;
                    index3 = l;
                    %R3_final = R3;
                end

                %fprintf('count: %d, obs nums: %d , %d , %d, residual: %0.3f  m\n',...
                %   endCounter, obs1.obsNumber, obs2.obsNumber, obs3.obsNumber, RMS);


%                [V1, V2, V3] = gibbs_solver(...
%                    R1, R2, R3, consts.mu_earth_km*1e9);

            end
%{
            % get residual RMS value for nth dataset
            nthResidual = 0;
        
            for m = 1:3
                observedAZ;
                predAZ;
                observedEL;
                predEL;
                nthResidual = nthResidual + sqrt((1/N)*((- nthAZ)^2 + (- nthEL)^2));
            end
        
            % Check if better solution exists:
            if (nthResidual < bestResidual)
                bestResidual = nthResidual;
                bestObservation = nthObservation;
            end
                %}
        end
    end    
end

% STOP TIMER:
timeOut = toc;

% Final observation triplet:

% First observation, n / "enn":
final1.obsNumber = sampleVector.obsNumber(index1);
final1.epoch = sampleVector.datetimes(index1);
final1.azimuth_deg = sampleVector.azimuths(index1);
final1.elevation_deg = sampleVector.elevations(index1);
final1.range_m = 5000;
final1.right_ascension_deg = sampleVector.ras(index1);
final1.declination_deg = sampleVector.declinations(index1);
chile1.position_m = chileObservatory.position_m;
chile1.epoch = final1.epoch;
chile1.latitude_deg = chileLLA.latitude_deg;
chile1.longitude_deg = chileLLA.longitude_deg;
chile1.altitude_m = chileLLA.altitude_m;
%final1ECI = eci(final1, chile1);
%final1ECI.position_m = final1ECI.position_m / norm(final1ECI.position_m);

% Second observation, k / "kay":
final2.obsNumber = sampleVector.obsNumber(index2);
final2.epoch = sampleVector.datetimes(index2);
final2.azimuth_deg = sampleVector.azimuths(index2);
final2.elevation_deg = sampleVector.elevations(index2);
final2.right_ascension_deg = sampleVector.ras(index2);
final2.declination_deg = sampleVector.declinations(index2);
final2.range_m = 5000;
chile2.position_m = chileObservatory.position_m;
chile2.epoch = final2.epoch;
chile1.latitude_deg = chileLLA.latitude_deg;
chile1.longitude_deg = chileLLA.longitude_deg;
chile1.altitude_m = chileLLA.altitude_m;
%final2ECI = eci(final2, chile2);
%final2ECI.position_m = final2ECI.position_m / norm(final2ECI.position_m);

% Third observation, l / "ell":
final3.obsNumber = sampleVector.obsNumber(index3);
final3.epoch = sampleVector.datetimes(index3);
final3.azimuth_deg = sampleVector.azimuths(index3);
final3.elevation_deg = sampleVector.elevations(index3);
final3.right_ascension_deg = sampleVector.ras(index3);
final3.declination_deg = sampleVector.declinations(index3);
final3.range_m = 5000;
chile3.position_m = chileObservatory.position_m;
chile3.epoch = final3.epoch;
chile3.latitude_deg = chileLLA.latitude_deg;
chile3.longitude_deg = chileLLA.longitude_deg;
chile3.altitude_m = chileLLA.altitude_m;
%final3ECI = eci(final3, chile3);
%final3ECI.position_m = final3ECI.position_m / norm(final3ECI.position_m);

t1 = second(final1.epoch, 'secondofday');
t2 = second(final2.epoch, 'secondofday');
t3 = second(final3.epoch, 'secondofday');

[finalLOS1, finalLOS2, finalLOS3] = get_los(...
    final1.right_ascension_deg, final1.declination_deg, ...
    final2.right_ascension_deg, final2.declination_deg, ...    
    final3.right_ascension_deg, final3.declination_deg);

[R1FinalEnding.position_m, R2FinalEnding.position_m, R3FinalEnding.position_m] = gauss_solver(...
    finalLOS1, finalLOS2, finalLOS3, ...
    chile1.position_m, chile2.position_m, chile3.position_m, ...
    t1, t2, t3);

endingGraphTitle = "Ending: good data";

try
    [V1FinalEnding, V2FinalEnding, V3FinalEnding] = gibbs_solver(...
        R1FinalEnding.position_m, R2FinalEnding.position_m, R3FinalEnding.position_m, ...
        consts.mu_earth_km*1e9);
catch ME

    warning('Bad position vectors in Ending data set.  Replacing w/ dummy set.');
    endingGraphTitle = "Ending: BAD data";

    % Load representative data set:
    goodData = load('goodObservations.mat');

    % Good observation triplet:
    R1FinalEnding.position_m = goodData.goodPositions(:,1);
    R2FinalEnding.position_m = goodData.goodPositions(:,2);
    R3FinalEnding.position_m = goodData.goodPositions(:,3);
    V1FinalEnding = goodData.goodVelocities(:,1);
    V2FinalEnding = goodData.goodVelocities(:,2);
    V3FinalEnding = goodData.goodVelocities(:,3);

end

waitbar(0.9, waitBar1, 'Almost finished...');
% Create histogram of data:
figure('Name', 'Histogram of End Set RMS Distribution');
histogram(endRMSDistribution);
title('Histogram of End Set RMS Distribution - Team 12, ENAE441');

xIn1 = [R1FinalEnding.position_m; V1FinalEnding];
xIn2 = [R2FinalEnding.position_m; V2FinalEnding];
xIn3 = [R3FinalEnding.position_m; V3FinalEnding];
tRange = [0 12*60*60];

[tOut, endingOrbit1] = ode45(@propagate_2BP, tRange, xIn1, [], consts.mu_earth_km*1e9);
[tOut, endingOrbit2] = ode45(@propagate_2BP, tRange, xIn2, [], consts.mu_earth_km*1e9);
[tOut, endingOrbit3] = ode45(@propagate_2BP, tRange, xIn3, [], consts.mu_earth_km*1e9);

%{
figure('Name', 'ENDING - Initial Orbit Plot');
plot3(endingOrbit1(:,1), endingOrbit1(:,2), endingOrbit1(:,3), 'Color', 'r');
hold on;
plot3(endingOrbit2(:,1), endingOrbit2(:,2), endingOrbit2(:,3), 'Color', 'g');
hold on;
plot3(endingOrbit3(:,1), endingOrbit3(:,2), endingOrbit3(:,3), 'Color', 'b');
title(endingGraphTitle);
%}


fprintf('\n%s\n', endingGraphTitle);

fprintf('\n\nBest Ending Triplet:  ( %d, %d, %d ) [obs. no.],  Min. RMS: %0.4f deg.\n', ...
    final1.obsNumber, final2.obsNumber, final3.obsNumber, minRMS);

fprintf('Ending Sample Size: %d , Time Elapsed: %0.1f seconds\n\n', ...
    sampleSize, timeOut);


%% OUTPUT GRAPHS:

figure('Name', 'Initial Orbit Plot');
subplot(4,1,1);
plot3(beginningOrbit1(:,1), beginningOrbit1(:,2), beginningOrbit1(:,3), 'Color', 'r');
hold on;
plot3(beginningOrbit2(:,1), beginningOrbit2(:,2), beginningOrbit2(:,3), 'Color', 'g');
hold on;
plot3(beginningOrbit3(:,1), beginningOrbit3(:,2), beginningOrbit3(:,3), 'Color', 'b');
title(beginningGraphTitle);
subplot(4,1,2);
plot3(endingOrbit1(:,1), endingOrbit1(:,2), endingOrbit1(:,3), 'Color', 'r');
hold on;
plot3(endingOrbit2(:,1), endingOrbit2(:,2), endingOrbit2(:,3), 'Color', 'g');
hold on;
plot3(endingOrbit3(:,1), endingOrbit3(:,2), endingOrbit3(:,3), 'Color', 'b');
title(endingGraphTitle);
subplot(4,1,3);
plot3(orbit1FinalMixed(:,1), orbit1FinalMixed(:,2), orbit1FinalMixed(:,3), 'Color', 'r');
hold on;
plot3(orbit2FinalMixed(:,1), orbit2FinalMixed(:,2), orbit2FinalMixed(:,3), 'Color', 'g');
hold on;
plot3(orbit2FinalMixed(:,1), orbit2FinalMixed(:,2), orbit2FinalMixed(:,3), 'Color', 'b');
title(titleOfMixedGraph);
subplot(4,1,4);
plot3(beginningOrbit1(:,1), beginningOrbit1(:,2), beginningOrbit1(:,3), 'Color', 'r');
hold on;
plot3(beginningOrbit2(:,1), beginningOrbit2(:,2), beginningOrbit2(:,3), 'Color', 'g');
hold on;
plot3(beginningOrbit3(:,1), beginningOrbit3(:,2), beginningOrbit3(:,3), 'Color', 'b');
hold on;
plot3(endingOrbit1(:,1), endingOrbit1(:,2), endingOrbit1(:,3), 'Color', 'r');
hold on;
plot3(endingOrbit2(:,1), endingOrbit2(:,2), endingOrbit2(:,3), 'Color', 'g');
hold on;
plot3(endingOrbit3(:,1), endingOrbit3(:,2), endingOrbit3(:,3), 'Color', 'b');
hold on;
plot3(orbit1FinalMixed(:,1), orbit1FinalMixed(:,2), orbit1FinalMixed(:,3), 'Color', 'r');
hold on;
plot3(orbit2FinalMixed(:,1), orbit2FinalMixed(:,2), orbit2FinalMixed(:,3), 'Color', 'g');
hold on;
plot3(orbit2FinalMixed(:,1), orbit2FinalMixed(:,2), orbit2FinalMixed(:,3), 'Color', 'b');

zeroIndicesAZ = azErr == 0;
azErr(zeroIndicesAZ) = [];
zeroIndicesEL = elErr == 0;
elErr(zeroIndicesEL) = [];

figure('Name', 'Azimuth Error Histogram');
histogram(azErr);

figure('Name', 'Elevation Error Histogram');
histogram(elErr);

figure ('Name', 'Comparison of AZ, EL errors over time'); 
plot(azErr(:), 'Color', 'r'); 
hold on; 
plot(elErr(:), 'Color', 'g');
hold on;
plot(endRMSDistribution(:), 'Color', 'k');
hold on;
plot(begRMSDistribution(:), 'Color', [0.22 0.33 0.44]);
hold on;
plot(mixedRMSDistribution(:), 'Color', [0.22 0.55 0.66]);