%% Random Split Observation Sample

% Get randomly-selected sample of N observations from
% beginning, middle, and end groups


function [begPartition, endPartition] = random_single_observation_sample(totalSize, sampleSize)

    % Set up partition limits:
    partitions = 20;
    
    % Beginning limits:
    partitionSize = floor(totalSize/partitions);
    begLim = [1 partitionSize];
    
    % Ending limits:
    endIndexOne = (partitions - 1)*floor(totalSize/partitions);
    endLim = [endIndexOne totalSize];        
    
    % Beginning group;
    for n = 1:sampleSize    
        collisionFlag = 1;
        while(collisionFlag)
            begNum = ceil(rand(1)*begLim(2));
            begPartition(n) = begNum;
            begCheck = unique(begPartition);
            if (length(begCheck) == length(begPartition))
                begPartition(n) = begNum;
                collisionFlag = 0;
            end
        end
    end
        
    begPartition = unique(begPartition);
    
    % End group;
    for n = 1:sampleSize    
        collisionFlag = 1;
        while(collisionFlag)
            endNum = endLim(1) + floor(rand(1)*partitionSize);
            endPartition(n) = endNum;
            endCheck = unique(endPartition);
            if (length(endCheck) == length(endPartition))
                endPartition(n) = endNum;
                collisionFlag = 0;
            end
        end
    end
        
    endPartition = unique(endPartition);
    
end