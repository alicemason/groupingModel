clear all
close all


theta=4; %free response stregth parameter
w= 0.05% free weighting parameter


demoted=8;
prom(1,:)=[17 1]; % partial
prom(2,:)=[9 9]; %full


nTrials=10;
j=1;
k=2;
run=1;


     
blocks=20;      
for p=1:2; %1=partial 2=full either running partial or full
v(j,:)=repmat((prom(p,:)),1,(nTrials*blocks/2));
v(k,:)=repmat(demoted,1,nTrials*blocks);
block=1;
counter=1;
while block<21
    
    for t=1:nTrials

        Q(j,1)=0;
        Q(k,1)=0;
        P(j,t)=1/(1+exp(theta*(Q(k,t)-Q(j,t)))); % probability of choosing j over   

        % propenitity to choose options k and j on next trial
        Q(j,t+1)=(1-w)*Q(j,t)+(w*v(j,t)); % if it was j that was observed
        Q(k,t+1)=(1-w)*Q(k,t)+(w*v(k,t)); % if it was k that was observed
        
        
    end
    % start extinction
    P_Prom(block)=mean(P)
    block=block+1;
    if block==11;
v(j,:)=repmat(1,1,nTrials*blocks);  
    end
    
    
end
if run==1
Partial_Prom=P_Prom;
run=run+1
elseif run==2
    Full_Prom=P_Prom;
end
end
hold on
plot(Full_Prom)
plot(Partial_Prom)
hold off



