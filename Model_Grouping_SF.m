
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
    
    % Here is what happens in model
    % We generated a longish vector of group sizes, then find first
    % group that takes us beyond list length. We truncate that group
    % and use only the groups we need
    
    groupSize=randsample(possGroupSize,listlength,true);% random
    
    % or constant within list, but varies across lists
    %groupSize=repmat(randsample(possGroupSize,1,true),1,listlength);
    
    cumulz = cumsum(groupSize);
    numGroups = find(cumulz>=listlength, 1, 'first'); % finds first instance
    groupSize(numGroups) = listlength-cumulz(numGroups-1);
    groupSize = groupSize(1:numGroups);
    
    lContext = ones(1,numGroups);     % control element for list context (each group has a list context)
    List_cue = 1; % cues the list - at moment set to last LIST
    eta_LC=1+randn(1,numGroups)*sigma_L; % Eq A3
    
    C_LC=eta_LC.*phi_l.^abs(List_cue-lContext); %Eq A7 - control element for list
    C_NC=zeros(1,numGroups); % control element for group cue
    C_NC(1)=eta_NC; % this shouldn't be list cue, but rather the
                            % cue to a particular group (let's say the
                            % first, but try setting it to last group)
    %C=zeros(1,numGroups); % don't need this, as it is overwritten in next line
    C=C_NC+C_LC; % list and group control elements added
    [max_value,max_idx]=max(C); % select most activated
    C(max_idx)=1; % set control element for most activared to 1
    % not sure what effect previous line has, as it isn't used below
    
    % make group markers
    gContext = [];
    pContext = [];
    absP = [];
    % set control elements
    
    for gz=1:length(groupSize)
        gContext = [gContext repmat(gz,1,groupSize(gz))];
        pContext = [pContext linspace(0,1,groupSize(gz))];
        absP = [absP 1:groupSize(gz)];
    end
    
    % commented out following line as it didn't make sense
    %item_probe=listlength-(sum(groupSize(max_idx:end)))+1; % want to probe for first item of cued group
    % instead, let's just use 1 for the moment, and soon we will
    % implement sequential retrieval
    item_probe = 1;
    
    Group_cue=max_idx;% conext cued of the most activated group
    P_CG=1; % assume no effect of time
    
    eta_gv=gamma.^(absP-1)+randn(1,listlength)*sigma_gp; %Eq A10
    v_GV = P_CG*eta_gv.*phi_g.^abs(Group_cue-gContext); %Eq A11
    v_PV = phi_p.^abs(pContext(item_probe)-pContext); %Eq A14
    
    % implemented primacy gradient
    v_PV=eta_gv.*v_PV;
    
    % sum group and item vectors to get activation of each item
    t_v = rho*v_GV + (1-rho)*v_PV; % Eq A15
    v(t,:)=t_v;
    % noisy retrieval Eq A16
    noise=randn(1,listlength)*sigma_v;
    a=(t_v+noise).*(1-r);
    % activation of two highest items
    [max_value,max_idx] = max(a);
    a(max_idx) = NaN;
    second_max = max(a);
    a = max_value; % retruns max_value into a
    
    if (max_value-second_max)>theta
        recalled_item(t)=max_idx;
    else
        recalled_item(t)=0;
    end
end

for i=1:listlength;
    prop(i)=numel(find(recalled_item==i))/nTrials;
end

plot(1:i,prop)

Av_v=mean(v);
Av_a=mean(a);
