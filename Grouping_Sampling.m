
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
                listlength=nTrials;
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

                        absP = [];
                        gContext=[];
                        groups=1;
                        groupSize=randsample(possGroupSize,listlength,true);% random
                        cumulz = cumsum(groupSize);
                        numGroups = find(cumulz>=listlength, 1, 'first'); % finds first instance
                        if numGroups>1
                            groupSize(numGroups) = listlength-cumulz(numGroups-1);
                            groupSize = groupSize(1:numGroups);
                        else
                            groupSize = listlength;
                        end

                        for gz=1:length(groupSize)
                            gContext = [gContext repmat(gz,1,groupSize(gz))];
                            absP = [absP 1:groupSize(gz)];
                        end


                        for block=1:nBlocks

                            if block > nTrainBlocks
                                curr_prom = prom_ext;
                            else
                                curr_prom = prom(p,:);
                            end


                            for t=1:nTrialsPerBlock
                                %SIZE of current group
                                payoff=[];
                                %current group size
                                Curr_group=groupSize(gContext(absT));                        
                                %if this is the first item of a group then payoff will be 0 -no info available    
                                if absP(absT)==1
                                       payoff=[];
                                %if the groupsize of the current group is
                                %bigger than 1
                                elseif groupSize(gContext(absT))>1      
                                Start_CG=1+(find((gContext(1:absT)~=gContext(absT)),1,'last'));
                                Index_CG= Start_CG:(absT-1);
                                payoff_CG=v(3,Index_CG); % find sequence for current group to t-1
                                sample_size=length(payoff_CG);
                                    % find all other groups that have been seen that are smaller than
                                    % or equal to size of current group

                                    for g=1:(gContext(absT-1))
                                        groupStart=find(gContext==g,1,'first');
                                        if groupSize(gContext(groupStart))>=sample_size
                                            sequence=strfind(v(3,groupStart:(groupStart+sample_size-1)),payoff_CG);
                                            if absP(sequence)==2 % if it's the second item of a group just compare to first items
                                                match=groupStart; % otherswise compare to sequecne length
                                            else
                                            match=groupStart:(groupStart+sample_size-1);
                                            end
                                        payoff=[payoff;match]; % build up payoff index for each time the sequence has been sampled
                                    end

                                    end
                                end

                                   if length(payoff)>0
                                        exVal = mean(v(1,payoff'));
                                        % choose option with highest value
                                        t_P(t,block) = exVal>demoted(Exp);
                                    else % if we don't have any matches default to "safe" option
                                        t_P(t,block) = 0;% demoted(Exp); % ...just guess
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

                    Full_P_prom(1:20)=mean(Full_Prom)
                    Partial_P_prom(1:20)=mean(Partial_Prom)


                    plot(Full_P_prom)
                    hold on
                    plot(Partial_P_prom)
                    hold off





