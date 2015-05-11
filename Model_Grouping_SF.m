
clear all
%Model Parameters
phi_p = 0.1;
phi_g=0.3;
phi_l=0.35;
gamma=0.97;% rate of change of primacy gradient across groups
sigma_gp=0.02;%
sigma_L=0.02;
sigma_v=0.005;
rho=0.3;
theta=0.003;
eta_NC=0.15;% learning rate for assoication between context and group

% experimental details
nTrials=1000;
listlength=12;
possGroupSize=1:5;

v=zeros(nTrials, listlength); % this can be easily pre-allocated, so should be
recalled_item=zeros(nTrials,1);

N_O=0; % counter for non recalled items
for t=1:nTrials
    
    r=zeros(1,listlength);
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
    
    % control element for list context
    lContext = ones(1,numGroups);
    %Ignore timing (assume, e.g., p_j^{CG} = 1)
    eta_LC=1+randn(numGroups)*sigma_L; % Eq A3
    
    List_cue = 1;
    % Simulate recall of current (last presented) list
    C_LC=eta_LC.*phi_l.^abs(List_cue-lContext); %Eq A7
    
    % Length; % association of each item to the group context
    C_NC=zeros(1,numGroups);
    
    % if group cue elemment was previously associated with control elemnent
    % takes eta_NC else will be 0
    % Assume people "tag" the first group - don't we want to assume they
    %tag last group?
    C_NC(Group_cue)=eta_NC;
    % I'm using last group as context here:
%     if C(Group_cue)==1
%         C_NC(Group_cue)=eta_NC;
%     else C(Group_cue)=0;
%     end
    
    % select most actiavted control element
    C=C_NC+C_LC;
    [max_value,max_idx]=max(C)
    C(max_idx)=1; % set control element for most activared to 1
    
    % make group markers
    gContext = [];
    pContext = [];
    absP = [];
    % set control elements
    C=zeros(1,numGroups);
    
    
    for gz=1:length(groupSize)
        gContext = [gContext repmat(gz,1,groupSize(gz))];
        pContext = [pContext linspace(0,1,groupSize(gz))];
        absP = [absP 1:groupSize(gz)];
    end
    
    item_probe=listlength-(groupSize(numGroups))+1; % want to probe for first item of last group
    Group_cue=gContext(item_probe);
    

    
    %scale by noisy learning paramter
    noise=randn(1,12)*sigma_gp; % watch out, this needs to be multiplication
    eta_gv=gamma.^(absP-1)+noise; % absP is the within-group ordinal position of
    v_GV = P_CG*eta_gv.*phi_g.^abs(gContext(item_probe)-gContext);
    
    v_PV = phi_p.^abs(pContext(item_probe)-pContext);
    
    % implemented primacy gradient
    v_PV=eta_gv.*v_PV;
    
    % sum group and item vectors to get activation of each item
    t_v = rho*v_GV + (1-rho)*v_PV;
    v(t,:)=t_v;
    %Next step: implement item recall;
    %for the moment, just assume recall of a single item using first item/last group as cue.
    %Take average of the recalls, which will effectively be a first recall probability function.
    
    
    noise=randn(1,listlength)*sigma_v;
    a=(t_v+noise).*(1-r);
    %a(t,:)=(v(t,:)+N)*(1-r(t,:));
    
    % activation of two highest items
    % largest
    [max_value,max_idx] = max(a);
    a(max_idx) = NaN;
    second_max = max(a);
    a = max_value; % not sure what this line does
    
    
    % rewrote the following slightly to save space
    if (max_value-second_max)>theta
        recalled_item(t)=max_idx;
    else
        % should put something like this in so you know what is going in to
        % vector
        recalled_item(t)=0;
    end
end

% I tried to vectorise this but I'm not sure it's possible with "find"
for i=1:listlength;
    prop(i)=numel(find(recalled_item==i))/nTrials;
end


plot(1:i,prop)

Av_v=mean(v);
Av_a=mean(a);
