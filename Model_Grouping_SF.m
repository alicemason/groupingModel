
clear all
%Model Parameters
phi_p = 0.1;
phi_g=0.3;
gamma=0.97% rate of change of primacy gradient across groups
sigma_gp=0.02%
rho=0.3;

% experimental details
nTrials=1000;
listlength=12;
possGroupSize=1:5;

v=zeros(nTrials, listlength); % this can be easily pre-allocated, so should be

for t=1:nTrials
    % Set the group sizes
    % moved following here
    %numGroups = 4; % number of groups - did we decide that there needed to be 4
% or should this be random too?
%    groupSize=randsample(possGroupSize,numGroups,true);
%     while (sum(groupSize)>listlength || sum(groupSize)<listlength)
%         groupSize=randsample(possGroupSize,numGroups,true);
%     end

    
    % Here is what happens in model
    % We generated a longish vector of group sizes, then find first
    % group that takes us beyond list length. We truncate that group
    % and use only the groups we need
    
    % completely random
    groupSize=randsample(possGroupSize,listlength,true);
    
    % or constant within list, but varies across lists
    %groupSize=repmat(randsample(possGroupSize,1,true),1,listlength);
    
    cumulz = cumsum(groupSize);
    numGroups = find(cumulz>=listlength, 1, 'first'); % finds first instance
    groupSize(numGroups) = listlength-cumulz(numGroups-1);
    groupSize = groupSize(1:numGroups);
    
    % make group markers
    
    gContext = [];
    pContext = [];
    
    for gz=1:length(groupSize)
        gContext = [gContext repmat(gz,1,groupSize(gz))];
        pContext = [pContext linspace(0,1,groupSize(gz))];
    end
    item_probe=listlength-(groupSize(numGroups))+1 % want to probe for first item of last group
                                                % (what you had was close)
    
    v_GV = phi_g.^abs(gContext(item_probe)-gContext);
    v_PV = phi_p.^abs(pContext(item_probe)-pContext);

    t_v = rho*v_GV + (1-rho)*v_PV;
    v(t,:)=t_v;
end
Av_v=mean(v);
plot(Av_v);

% simulation with variable group sizes
% simulate individual trials (list=12) 100 trials
% on each trial you will have a different v (put into maxtrix)
% within the trial (ie 12 items) have different size groups
% group size vector uniform distrubution from 1-5
% look at paper to see how done in paper
% known list length


% noisy retrieval