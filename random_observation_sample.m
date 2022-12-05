%% Random Observation Sample

% Get randomly-selected sample of N observations from
% beginning, middle, and end groups


function [begSample, midSample, endSample] = random_observation_sample(total, N, split)

    beginningSize = floor(total/split);
    middleSize = ceil(total/split);
    endingSize = floor(total/split);
    
    count = 0;
    group = 0;
    
    g1Num = 0;
    g1Count = 1;
    g1Check = 0;
    g2Num = 0;
    g2Count = 1;
    g3Num = 0;
    g3Count = 1;
    
    for n = 1:N
    
        g1Flag = 1;
        g2Flag = 1;
        g3Flag = 1;
    
        group = ceil(rand(1)*3);
    
        if (group == 1)
            while (g1Flag)
                g1Num = ceil(rand(1)*beginningSize);
                g1List(g1Count) = g1Num;
                g1Check = unique(g1List);
                if (length(g1Check) == length(g1List))
                    g1List(g1Count) = g1Num;
                    g1Count = g1Count + 1;
                    g1Flag = 0;
                end
            end        
          %  fprintf('%d\n', opt2satCset4.observation_number(g1Num));
        elseif (group == 2)
            while (g2Flag)
                g2Num = ceil(beginningSize + rand(1)*middleSize);
                g2List(g2Count) = g2Num;
                g2Check = unique(g2List);
                if (length(g2Check) == length(g2List))
                    g2List(g2Count) = g2Num;
                    g2Count = g2Count + 1;
                    g2Flag = 0;
                end
            end        
        else
            while (g3Flag)
                g3Num = ceil(beginningSize + middleSize + rand(1)*endingSize);
                g3List(g3Count) = g3Num;
                g3Check = unique(g3List);
                if (length(g3Check) == length(g3List))
                    g3List(g3Count) = g3Num;
                    g3Count = g3Count + 1;
                    g3Flag = 0;
                end
            end        
        end
    end
    
    g1Count = g1Count - 1;
    g1List = unique(g1List);
    g2Count = g2Count - 1;
    g2List = unique(g2List);
    g3Count = g3Count - 1;
    g3List = unique(g3List);
    
    if ((g1Count + g2Count + g3Count) == N)
        fprintf('\n\nIt worked!\n\n');
        begSample = g1List;
        midSample = g2List;
        endSample = g3List;
    end

end