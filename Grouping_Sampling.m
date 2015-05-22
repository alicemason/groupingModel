
    Exp=1;
    demoted=[8 2]; % demoted option points

    %promoted options during training
    prom(1,:)=[17 1]; % partial reinforcement
    prom(2,:)=[9 9]; %full reinforce

    prom_ext = [1 1]; % promoted points during extinction

    nTrialsPerBlock=10; % number of trials per block
    listlength=nTrialsPerBlock;

    j=1;
    k=2;
    run=1;

    nTrainBlocks = 10;
    nTestBlocks = 10; % these are really extinction blocks; should prob rename

    nBlocks = nTrainBlocks + nTestBlocks;
    nTrials = nBlocks * nTrialsPerBlock;

    Partial_Prom = zeros(1,nBlocks);
    Full_Prom = zeros(1,nBlocks);

    possGroupSize=[2 3 4 5];

    nReps=1000;


    for p=1:2; %1=partial 2=full either running partial or full
        % nesting p inside reps didn't make sense to me, as full and
        % partial are different conditions (ppl don't see same thing)

        for reps=1:nReps

            t = 1;
            counter=1;
            absT=1;

            vTrain = [repmat(prom(p,1),1,nTrialsPerBlock*nTrainBlocks/2) ...
                repmat(prom(p,2),1,nTrialsPerBlock*nTrainBlocks/2)];
            vExt = [repmat(prom_ext(1),1,nTrialsPerBlock*nTestBlocks/2) ...
                repmat(prom_ext(2),1,nTrialsPerBlock*nTestBlocks/2)];

            v = [randsample(vTrain,length(vTrain),false) ...
                randsample(vExt,length(vExt),false)];
            v = [v; repmat(demoted(Exp),1,nTrials)];
            [y, maxi] = max(v);
            v = [v; maxi==1];
            start_pos=1;
            for block=1:nBlocks

                if block > nTrainBlocks
                    curr_prom = prom_ext;
                else
                    curr_prom = prom(p,:);
                end

                groupSize=randsample(possGroupSize,listlength,true);% random
                cumulz = cumsum(groupSize);
                numGroups = find(cumulz>=listlength, 1, 'first'); % finds first instance
                if numGroups>1
                    groupSize(numGroups) = listlength-cumulz(numGroups-1);
                    groupSize = groupSize(1:numGroups);
                else
                    groupSize = listlength;
                end


                Curr_List=1;
                absP = [];
                gContext=[];

                for gz=1:length(groupSize)
                    gContext = [gContext repmat(gz,1,groupSize(gz))];
                    absP = [absP 1:groupSize(gz)];
                end

    for t=1:nTrialsPerBlock

                    Curr_group=gContext(t);

                    % take the sequence so far from the current group - N.B if it's
                    % start of a new group - will sample position 1 from
                    % previous groups
                    Index_CG= strfind(gContext(1:(t-1)),Curr_group);
                    payoff_CG=v(3,Index_CG-1+start_pos); % find sequence for current group to t-1
                    sample_size=length(payoff_CG);
                    % compare current payoff to previous groups only if the
                    % previous group was as big as sample size
                   payoff=[];
                    if Curr_group>1
                        for g=1:(Curr_group-1);
                            groupStart=find(gContext==g,1,'first');
                            groupEnd=groupStart+(groupSize(g)-1);
                            %if the group size is bigger than current sample we
                            %can find payoff for sample+1
                          if groupSize(g)>sample_size 
                              % this is not very nice as I'm trying to index
                              % the trials in v
                              sequence=strfind(v(3,start_pos+(groupStart-1):(start_pos+(groupStart-1)+sample_size-1)),payoff_CG);
                         match=start_pos+(groupStart-1)+1;

                          end
                                              payoff=[payoff;match];    

                    end
     if absP(t)==1
     payoff=start_pos+(groupStart-1); 
                    end
                    end
                    % what to do for first item of a new group






                    if length(payoff)>0
                        exVal = mean(v(1,payoff'));
                        % choose option with highest value
                        t_P(t,block) = exVal>demoted(Exp);
                    else % if we don't have any matches default to "safe" option
                        t_P(t,block) = 0;% demoted(Exp); % ...just guess
                    end

                    absT=absT+1;
                end
                start_pos=(listlength*block)+1;

                Curr_List=Curr_List+1;

                if p==1
                    Partial_Prom(reps,block) = mean(t_P(:,block));
                else
                    Full_Prom(reps,block) = mean(t_P(:,block));
                end
            end

        end
    end

    Full_P_prom(1:20)=mean(Full_Prom)
    Partial_P_prom(1:20)=mean(Partial_Prom)


    plot(Full_P_prom)
    hold on
    plot(Partial_P_prom)
    hold off





