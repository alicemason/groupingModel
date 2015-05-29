
clear all; close all
% contingent-sampling model in which similarity was defined
%on the basis of the relative advantage of
%the promoted option in the m most recent trials
Exp=2;
demoted=[8 2]; % demoted option points

%promoted options during training
prom(1,:)=[17 1]; % partial reinforcement
prom(2,:)=[9 9]; %full reinforce

prom_ext = [1 1]; % promoted points during extinction

nTrialsPerBlock=10; % number of trials per block
j=1;
k=2;
run=1;

nTrainBlocks = 10;
nTestBlocks = 10; % these are really extinction blocks; should prob rename

nBlocks = nTrainBlocks + nTestBlocks;
nTrials = nBlocks * nTrialsPerBlock;

Partial_Prom = zeros(1,nBlocks);
Full_Prom = zeros(1,nBlocks);

groupSizes = [2 3 4 5];
recScale = 0; % set to 0 for no recency

nReps=200;

for p=1:2; %1=partial 2=full either running partial or full
    % nesting p inside reps didn't make sense to me, as full and
    % partial are different conditions (ppl don't see same thing)
    
    for reps=1:nReps
        
        t = 1;
        counter=1;
        absT=1;
        
        % this randomizes feedback while controlling overall
        % frequency
        vTrain = [repmat(prom(p,1),1,nTrialsPerBlock*nTrainBlocks/2) ...
            repmat(prom(p,2),1,nTrialsPerBlock*nTrainBlocks/2)];
        vExt = [repmat(prom_ext(1),1,nTrialsPerBlock*nTestBlocks/2) ...
            repmat(prom_ext(2),1,nTrialsPerBlock*nTestBlocks/2)];
        
        v = [randsample(vTrain,length(vTrain),false) ...
            randsample(vExt,length(vExt),false)];
        v = [v; repmat(demoted(Exp),1,nTrials)];
        [y, maxi] = max(v);
        v = [v; maxi==1];
        
        newGroup = 1;
        
        groupHist = NaN(nTrials,max(groupSizes));
        groupVals = NaN(nTrials,max(groupSizes));
        groupHistSize = zeros(nTrials,1);
        groupIndex = 0;
        
        for block=1:nBlocks
            
            if block > nTrainBlocks
                curr_prom = prom_ext;
            else
                curr_prom = prom(p,:);
            end
            
            %t_P = zeros(1,nTrialsPerBlock);
            
            for t=1:nTrialsPerBlock
               
                if newGroup==1
                    currGroupSize = randsample(groupSizes,1);
                    inpos = 1;
                    newGroup=0;
                    currGroup = NaN(1,currGroupSize);
                    currVals = NaN(1,currGroupSize);
                    groupIndex = groupIndex + 1;
                end
                
                if inpos>1
                    matches = sum(repmat(currGroup(1:(inpos-1)),nTrials,1)==...
                        groupHist(:,1:(inpos-1)),2)==(inpos-1) & groupHistSize>=inpos;
                else
                    matches = groupHistSize>=1;
                    %matches = 0; % forces guess
                end
                
                % 'strfind' is a good choice here (how I would have
                % done it) :)
                %matches=strfind(v(3,1:(absT-2)),sample); % find the start index of a sample
                % made this -2 to make sure it doesn't include current
                % choice in matches
                
                % We want the expected payoff given that sequence;
                % so average rewards following each match
                if sum(matches)>0
                    weights = exp(recScale.*((1:(groupIndex-1))-(groupIndex-1)));
                    weights(groupIndex:nTrials) = 0;
                    weights = weights(matches)./sum(weights(matches));
                    exVal = sum(groupVals(matches,inpos).*weights');
                    
                    % choose option with highest value
                    t_P(t,block) = exVal>demoted(Exp);
                    
                    % has trouble going to 0 for full condition as matches
                    % all previous trials, and expected value for promoted
                    % option in training phase often outweighs extinction
                    % phase
                    
                else % if we don't have any matches...
                    t_P(t,block) = rand>0.5; % ...just guess
                    %t_P(t,block) = 0; % pick demoted
                end
                
                % Note that this (the method in their paper) only
                % works when demoted option is const. Could be
                % extended to where demoted also varies
                
                
                %                         if demoted>score;
                %                             t_P(t,block)=0;
                %                         else
                %                             t_P(t,block)=1;
                %                         end
                %                         t_P(1:m,1)=0;
                
                currGroup(inpos) = v(3,absT);
                currVals(inpos) = v(1,absT);
                
                absT=absT+1;
                
                inpos=inpos+1;
                
                if inpos>currGroupSize
                    newGroup=1;
                    groupHist(groupIndex,1:(inpos-1)) = currGroup;
                    groupVals(groupIndex,1:(inpos-1)) = currVals;
                    groupHistSize(groupIndex) = currGroupSize;
                end
            end
            
            if p==1
                Partial_Prom(reps,block) = mean(t_P(:,block));
            else
                Full_Prom(reps,block) = mean(t_P(:,block));
            end
        end
    end
end

plot(mean(Partial_Prom))
hold on
%    subplot(1,5,m);
plot(mean(Full_Prom))
hold off
ylim([0 1])




