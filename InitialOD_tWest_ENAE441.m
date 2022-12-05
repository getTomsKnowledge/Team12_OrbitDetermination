%% Initial Orbit Determination - Tom West, Team 12, ENAE441

%{

Description:
Observations occurred over two nights: 10/29/2020, 10/30/2020.
Geosynchronous Equatorial Orbits.
4 sets of data, 6 satellites (24 in total)

Objective:
Determine best measurement set for OD by minimizing RMS value

%}

%% HOUSEKEEPING

run('TidyWorkspace.m');
run('CheckSNaG.m');
printing = PrintFormatting();
heading = HeadingGenerator('ENAE441', 'Orbit Determination Project', 'Determine GEO satellite orbit based on observation data');
heading.print_heading;
consts = OrbitConstants();

%% DATA SAMPLING/SORTING

% Get orbit data from datasets:
load('opt2satCset4.mat');
load('opt3satCset4.mat');

% Split data set into column vectors for speed/clarity:
satC.observationNumbers = opt2satCset4.observation_number;
satC.datetimeList = opt2satCset4.datetime;
satC.azimuths = opt2satCset4.azimuth_deg;
satC.elevations = opt2satCset4.elevation_deg;
satC.ras = opt2satCset4.right_ascension_deg;
satC.declinations = opt2satCset4.declination_deg;

% Set up initial epoch:
initialEpoch = opt2satCset4.datetime(1);
chileLLA.epoch = initialEpoch;

% Set up site LLA:
chileLLA.latitude_deg = opt2satCset4.site_latitude_deg(1);
chileLLA.longitude_deg = opt2satCset4.site_longitude_deg(1);
chileLLA.altitude_m = opt2satCset4.site_altitude_m(1);
chileObservatory = eci(chileLLA);



% PREPARE RANDOM SAMPLE:

% Get random sample of size N from observation set:
totalOrbits = length(opt2satCset4.datetime);
sampleSize = 15;
numberOfGroups = 3;
K = 6;

[begSamp, midSamp, endSamp] = random_observation_sample(totalOrbits, sampleSize, numberOfGroups);
sampleList = [begSamp, midSamp, endSamp];
%sampleVector = zeros(sampleSize, K);

for n = 1:sampleSize
    sampleVector.obsNumber(n) = satC.observationNumbers(sampleList(n));
    sampleVector.datetimes(n) = satC.datetimeList(sampleList(n));
    sampleVector.azimuths(n) = satC.azimuths(sampleList(n));
    sampleVector.elevations(n) = satC.elevations(sampleList(n));
    sampleVector.ras(n) = satC.ras(sampleList(n));
    sampleVector.declinations(n) = satC.declinations(sampleList(n));
end

counter = 0;

% Comparative accuracy of sample combinations:
% START TIMER:

minRMS = 1*1e8;
tic

for n = 1:(sampleSize-2)    
    for k = 1:(sampleSize - 1)
        for l = 1:sampleSize

            if ((k >= (n + 1)) && ((l >= (n + 2)) && (l >= (k + 1))))
                counter = counter + 1;

                % Get nth observation set:

                % First observation, n / "enn":
                obs1.obsNumber = sampleVector.obsNumber(n);
                obs1.epoch = sampleVector.datetimes(n);
                obs1.azimuth_deg = sampleVector.azimuths(n);
                obs1.elevation_deg = sampleVector.elevations(n);
                obs1.range_m = 5000;
                %obs1.right_ascension_deg = sampleVector.ras(n);
                %obs1.declination_deg = sampleVector.declinations(n);
                chile1.position_m = chileObservatory.position_m;
                chile1.epoch = obs1.epoch;
                obs1ECI = eci(obs1, chileLLA);
                obs1ECI.position_m = obs1ECI.position_m / norm(obs1ECI.position_m);
                
                t1 = second(chile1.epoch, 'secondofday');


                % Second observation, k / "kay":
                obs2.obsNumber = sampleVector.obsNumber(k);
                obs2.epoch = sampleVector.datetimes(k);
                obs2.azimuth_deg = sampleVector.azimuths(k);
                obs2.elevation_deg = sampleVector.elevations(k);
                obs2.range_m = 5000;
                %obs2.right_ascension_deg = sampleVector.ras(k);
                %obs2.declination_deg = sampleVector.declinations(k);
                chile2.position_m = chileObservatory.position_m;
                chile2.epoch = obs2.epoch;
                obs2ECI = eci(obs2, chileLLA);
                obs2ECI.position_m = obs2ECI.position_m / norm(obs2ECI.position_m);

                t2 = second(chile2.epoch, 'secondofday');

                % Third observation, l / "ell":
                obs3.obsNumber = sampleVector.obsNumber(l);
                obs3.epoch = sampleVector.datetimes(l);
                obs3.azimuth_deg = sampleVector.azimuths(l);
                obs3.elevation_deg = sampleVector.elevations(l);
                obs3.range_m = 5000;
                %obs3.right_ascension_deg = sampleVector.ras(l);
                %obs3.declination_deg = sampleVector.declinations(l);
                chile3.position_m = chileObservatory.position_m;
                chile3.epoch = obs3.epoch;
                obs3ECI = eci(obs3, chileLLA);
                obs3ECI.position_m = obs3ECI.position_m / norm(obs3ECI.position_m);

                t3 = second(chile3.epoch, 'secondofday');


                [o1out.position_m, o2out.position_m, o3out.position_m] = gauss_solver(...
                    obs1ECI.position_m, obs2ECI.position_m, obs3ECI.position_m, ...
                    chile1.position_m, chile2.position_m, chile3.position_m, ...
                    t1, t2, t3);

                o1out.epoch = chile1.epoch;
                o2out.epoch = chile2.epoch;
                o3out.epoch = chile3.epoch;
                
                chileLLA.epoch = o1out.epoch;
                o1out.aer = aer(o1out, chileLLA);
                chileLLA.epoch = o2out.epoch;
                o2out.aer = aer(o2out, chileLLA);
                chileLLA.epoch = o3out.epoch;
                o3out.aer = aer(o3out, chileLLA);

                RMS = (obs1.azimuth_deg - o1out.aer.azimuth_deg)^2 + (obs1.elevation_deg - o1out.aer.elevation_deg)^2; 
                RMS = (obs2.azimuth_deg - o2out.aer.azimuth_deg)^2 + (obs2.elevation_deg - o2out.aer.elevation_deg)^2;
                RMS = (obs3.azimuth_deg - o3out.aer.azimuth_deg)^2 + (obs3.elevation_deg - o3out.aer.elevation_deg)^2;
                RMS = sqrt( RMS/3 );

                if (RMS < minRMS)
                    minRMS = RMS;
                    index1 = n;
                    %R1_final = R1;
                    index2 = k;
                    %R2_final = R2;
                    index3 = l;
                    %R3_final = R3;
                end

                fprintf('count: %d, obs nums: %d , %d , %d, residual: %0.3f  m\n',...
                    counter, obs1.obsNumber, obs2.obsNumber, obs3.obsNumber, RMS);


%                [V1, V2, V3] = gibbs_solver(...
%                    R1, R2, R3, consts.mu_earth_km*1e3);

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
final1.right_ascension_deg = sampleVector.ras(index1);
final1.declination_deg = sampleVector.declinations(index1);
chile1.position_m = chileObservatory.position_m;
chile1.epoch = final1.epoch;
%final1ECI = eci(final1, chileLLA);
%final1ECI.position_m = final1ECI.position_m / norm(final1ECI.position_m);

% Second observation, k / "kay":
final2.obsNumber = sampleVector.obsNumber(index2);
final2.epoch = sampleVector.datetimes(index2);
final2.azimuth_deg = sampleVector.azimuths(index2);
final2.elevation_deg = sampleVector.elevations(index2);
final2.right_ascension_deg = sampleVector.ras(index2);
final2.declination_deg = sampleVector.declinations(index2);
chile2.position_m = chileObservatory.position_m;
chile2.epoch = final2.epoch;
%final2ECI = eci(final2, chileLLA);
%final2ECI.position_m = final2ECI.position_m / norm(final2ECI.position_m);

% Third observation, l / "ell":
final3.obsNumber = sampleVector.obsNumber(index3);
final3.epoch = sampleVector.datetimes(index3);
final3.azimuth_deg = sampleVector.azimuths(index3);
final3.elevation_deg = sampleVector.elevations(index3);
final3.right_ascension_deg = sampleVector.ras(index3);
final3.declination_deg = sampleVector.declinations(index3);
chile3.position_m = chileObservatory.position_m;
chile3.epoch = final3.epoch;
%final3ECI = eci(final3, chileLLA);
%final3ECI.position_m = final3ECI.position_m / norm(final3ECI.position_m);

printf('\n\n\nBest Observation Triplet:  ( %d, %d, %d ) [obs. no.],  Min. RMS: %0.4f meters\n', ...
    final1.obsNumber, final2.obsNumber, final3.obsNumber, minRMS);


fprintf('Sample Size: %d ,  %s\n\n', ...
    sampleSize, timeOut);
%% Angles-Only Vector Solution (Gauss' Method)


%% Orbit Determination (Gibbs' Method)