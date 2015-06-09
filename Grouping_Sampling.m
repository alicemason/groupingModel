close all; clear all
w_CL=1;
w_PL=1;
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
listlength=nTrials/2;
Partial_Prom = zeros(1,nBlocks);
Full_Prom = zeros(1,nBlocks);

possGroupSize=[2 3 4 5];

nReps=200;


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
        groupSize=[];
        
        for list=1:2
            make_absP=[];
            make_gContext=[];
            GS=[];
            GS=randsample(possGroupSize,listlength,true);% random
            cumulz = cumsum(GS);
            numGroups = find(cumulz>=listlength, 1, 'first'); % finds first instance
            if numGroups>1
                GS(numGroups) = listlength-cumulz(numGroups-1);
                GS = GS(1:numGroups);
            else
                GS = listlength;
            end
            
            
            %             for gz=1:length(GS)
            %                 make_gContext = [make_gContext repmat(gz,1,GS(gz))];
            %                 make_absP= [make_absP 1:GS(gz)];
            %             end
            %             gContext=[gContext make_gContext];
            %             absP=[absP make_absP];
            %                        groupSize=[groupSize GS];
            
            
            groupSize=[groupSize GS];
            
            
            
        end
        for gz=1:length(groupSize)
            gContext = [gContext repmat(gz,1,groupSize(gz))];
            absP= [absP 1:groupSize(gz)];
        end
        
        
        for block=1:nBlocks
            
            if block > nTrainBlocks
                curr_prom = prom_ext;
                list=2;
                Start_CL=101;
            else
                curr_prom = prom(p,:);
                list=1;
                Start_CL=1;
            end
            
            
            
            for t=1:nTrialsPerBlock
                Curr_group=groupSize(gContext(absT));
                payoff=[];  
                List_ID(absT)=list;
                weight=1;
                payoff_L=[];

                if absP(absT)==1
                    payoff=rand>0.5;
                    
                elseif groupSize(gContext(absT))>1
                    Start_CG=1+(find((gContext(1:absT)~=gContext(absT)),1,'last'));
                    Index_CG= Start_CG:(absT-1);
                    payoff_CG=v(3,Index_CG);
                    sample_size=length(payoff_CG);
                    
                    
                    for g=1:(gContext(absT-1))
                        groupStart=find(gContext==g,1,'first');
                        if  groupSize(gContext(groupStart))>=sample_size
                            sequence=strfind(v(3,groupStart:(groupStart+sample_size-1)),payoff_CG);
                            if sequence==1
                                match=groupStart+sample_size;
                                payoff=[payoff;match];
                                payoff_L=List_ID(payoff);% build up payoff index for each time the sequence has been sampled
                            end
                        end
                        
                    end
                end
                
               for i=1:length(payoff_L)
                    if payoff_L(i)==list
                        weight(i)=w_CL;
                    else 
                        weigth(i)=w_PL;
                    end
               end
                
                    
                if length(payoff)>0
                    Val =v(1,payoff');
                    exVal=mean(Val.*weight);
                    % choose ption with highest value
                    t_P(t,block) = exVal>demoted(Exp);
                else % if we don't have any matches default to "safe" option
                    %t_P(t,block) = 0;% demoted(Exp); % ...just guess
                    t_P(t,block) = rand>.5;
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


plot(mean(Partial_Prom))
hold on
%    subplot(1,5,m);
plot(mean(Full_Prom))
hold off
ylim([0 1])
