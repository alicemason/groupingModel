clear all
close all

sigma=4; %free response stregth parameter
w= 0.05; % free weighting parameter


demoted=8; % demoted option points

%promoted options during training
prom(1,:)=[17 1]; % partial reinforcement
prom(2,:)=[9 9]; %full reinforce

prom_ext = [1 1]; % promoted points during extinction

nTrialsPerBlock=10; % number of trials per block
j=1;
k=2;
run=1;

nTrainBlocks = 10;
nTestBlocks = 10;

nBlocks = nTrainBlocks + nTestBlocks;
nTrials = nBlocks * nTrialsPerBlock;

Partial_Prom = zeros(1,nBlocks);
Full_Prom = zeros(1,nBlocks);
nReps=1000;

for reps=1:nReps
t = 1;
counter=1;

for p=1:2; %1=partial 2=full either running partial or full

    
    Q = zeros(2,1); % there are 2 options
    
    for block=1:nBlocks
        
        if block > nTrainBlocks
            curr_prom = prom_ext;
        else
            curr_prom = prom(p,:);
        end
        
        t_P = zeros(1,nTrialsPerBlock);
        
        for t=1:nTrialsPerBlock
            
            % v's need to be sampled--they are random
            v(1) = randsample(curr_prom,1);
            if v(1)==prom(1,1)
            counter=counter+1
            end
            if counter>51
                v(1)=prom(1,2)
            end
            
            v(2) = demoted;
            
            t_P(t)=1/(1+exp(sigma*(Q(2)-Q(1)))); % probability of choosing j over

            Q = (1-w).*Q + w.*v';
            
        end
        
        if p==1
            Partial_Prom(reps,block) = mean(t_P);
        else
            Full_Prom(reps,block) = mean(t_P);
        end
    end

end
end
 hold on
 plot(mean(Full_Prom))
 plot(mean(Partial_Prom))
 hold off
% 


