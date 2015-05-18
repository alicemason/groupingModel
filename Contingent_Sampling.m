    
clear all
 % contingent-sampling model in which similarity was defined
    %on the basis of the relative advantage of
    %the promoted option in the m most recent trials

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
    nReps=2000;
 for   m=1:5
        for reps=1:nReps
            t = 1;
            counter=1;
absT=1;
            for p=1:2; %1=partial 2=full either running partial or full

                for block=1:nBlocks

                    if block > nTrainBlocks
                        curr_prom = prom_ext;
                    else
                        curr_prom = prom(p,:);
                    end

                    %t_P = zeros(1,nTrialsPerBlock);

                    for t=1:nTrialsPerBlock
                            % v's need to be sampled--they are random
                            v(1,absT) = randsample(curr_prom,1);
                            if v(1,absT)==curr_prom(1,1);
                                counter=counter+1;
                            end
                            if counter>51
                                v(1,absT)=curr_prom(1,2);
                            end

                            v(2,absT) = demoted;

                            if v(1,absT)>v(2,absT);
                                v(3,absT)=1;% promoted is Higher
                            else
                                v(3,absT)=0;%promoted is Lower
                            end

                     if (absT-1) >=m;
                                start_pos=absT-m; % only use m-back if this is greater than one, else start at one
                            


                            sample=v(3,start_pos:absT-1);% all the times the promoted option was better
                            matches=strfind(v(3,1:(absT-1)),sample); % find the start index of a sample
score=[];
                            for i=1:length(matches);

                                % add up the payoff for each instance of the
                                % sequence - then get avergae 
                                % this is really nasty -sorry
                                
                                % gives you total payoff for instance i of
                                % a sequence
                               score(i)=sum(v(1,(matches(i):(matches(i)+(m-1)))));
 
                          
  
                            end
                            
 % get average payoff - divide by i and by
                                % number of items in sequence
                                score=mean(score)/numel(sample);
                                
                            if demoted>score;
                                t_P(t,block)=0;
                            else
                                t_P(t,block)=1;
                            end
                        t_P(1:m,1)=0;

                     end
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
   
   Full_P_prom(m,1:20)=mean(Full_Prom(:,1:20));
   Partial_P_prom(m,1:20)=mean(Partial_Prom(:,1:20));
   
   subplot(1,5,m);
   
plot(Full_P_prom(m,:))
 hold on
%    subplot(1,5,m);
   plot(Partial_P_prom(m,:))
   hold off

    end



