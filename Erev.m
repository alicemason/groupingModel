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

t = 1;

for p=1:2; %1=partial 2=full either running partial or full

    counter=1;
    
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
            v(2) = demoted;
            
            t_P(t)=1/(1+exp(sigma*(Q(2)-Q(1)))); % probability of choosing j over
            
            % propenitity to choose options k and j on next trial
%             Q(j,t+1)=(1-w)*Q(j,t)+(w*v(j,t)); % if it was j that was observed
%             Q(k,t+1)=(1-w)*Q(k,t)+(w*v(k,t)); % if it was k that was observed

            % people see all outcomes
            Q = (1-w).*Q + w.*v';
            
        end
        
        if p==1
            Partial_Prom(block) = mean(t_P);
        else
            Full_Prom(block) = mean(t_P);
        end
    end
%     if run==1
%         Partial_Prom=P_Prom;
%         run=run+1
%     elseif run==2
%         Full_Prom=P_Prom;
%     end
end
hold on
plot(Full_Prom)
plot(Partial_Prom)
hold off



