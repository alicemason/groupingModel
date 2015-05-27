
clear all
% contingent-sampling model in which similarity was defined
%on the basis of the relative advantage of
%the promoted option in the m most recent trials
Exp=2;
demoted=[8 2]; % demoted option points

%promoted options during training
prom(1,:)=[17 1]; % partial reinforcement
prom(2,:)=[9 9]; %full reinforce

prom_ext = [11]; % promoted points during extinction

nTrialsPerBlock=10; % number of trials per block
j=1;
k=2;
run=1;

nTrainBlocks = 10;
nTestBlocks = 20; % these are really extinction blocks; should prob rename

nBlocks = nTrainBlocks + nTestBlocks;
nTrials = nBlocks * nTrialsPerBlock;

Partial_Prom = zeros(1,nBlocks);
Full_Prom = zeros(1,nBlocks);

nReps=200;

for m=1:4
    
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
            
            for block=1:nBlocks
                
                if block > nTrainBlocks
                    curr_prom = prom_ext;
                else
                    curr_prom = prom(p,:);
                end
                
                %t_P = zeros(1,nTrialsPerBlock);
                
                for t=1:nTrialsPerBlock
                    % v's need to be sampled--they are random
                    
                    %                     v(1,absT) = randsample(curr_prom,1);
                    %
                    %                     if v(1,absT)==curr_prom(1,1);
                    %                         counter=counter+1;
                    %                     end
                    %                     if counter>51
                    %                         v(1,absT)=curr_prom(1,2);
                    %                     end
                    %
                    %                     v(2,absT) = demoted;
                    %
                    %                     if v(1,absT)>v(2,absT);
                    %                         v(3,absT)=1;% promoted is Higher
                    %                     else
                    %                         v(3,absT)=0;%promoted is Lower
                    %                     end
                    
                    start_pos=max(1,absT-m); % only use m-back if this is greater than one, else start at one
                    
                    sample=v(3,start_pos:absT-1);% all the times the promoted option was better
                    
                    
                    % 'strfind' is a good choice here (how I would have
                    % done it) :)
                    matches=strfind(v(3,1:(absT-2)),sample); % find the start index of a sample
                    % made this -2 to make sure it doesn't include current
                    % choice in matches
                    
                    % I don't think the next bit is correct, and have
                    % commented out; see my version below
                    %                         score=[];
                    %                         for i=1:length(matches);
                    %
                    %                             % add up the payoff for each instance of the
                    %                             % sequence - then get avergae
                    %                             % this is really nasty -sorry
                    %
                    %                             % gives you total payoff for instance i of
                    %                             % a sequence
                    %                             score(i)=sum(v(1,(matches(i):(matches(i)+(m-1)))));
                    %
                    %
                    %
                    %                         end
                    %
                    %                         % get average payoff - divide by i and by
                    %                         % number of items in sequence
                    %                         score=mean(score)/numel(sample);
                    
                    % We want the expected payoff given that sequence;
                    % so average rewards following each match
                    if length(matches)>0
                        exVal = mean(v(1,matches+m));
                    
                        % choose option with highest value
                        t_P(t,block) = exVal>demoted(Exp);
                    else % if we don't have any matches...
                        t_P(t,block) = rand>0.5; % ...just guess
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
                    
                    absT=absT+1;
                end
                
                
                if p==1
                    Partial_Prom(reps,block) = mean(t_P(:,block));
                else
                    Full_Prom(reps,block) = mean(t_P(:,block));
                end
            end
            
        end
    end
    
    Full_P_prom(m,1:nBlocks)=mean(Full_Prom(:,1:nBlocks));
    Partial_P_prom(m,1:nBlocks)=mean(Partial_Prom(:,1:nBlocks));
    
    subplot(1,5,m);
    
    plot(Full_P_prom(m,:))
    hold on
    %    subplot(1,5,m);
    plot(Partial_P_prom(m,:))
    hold off
    
end



