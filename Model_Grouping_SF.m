
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

%rec=zeros(nTrials, listlength); % this can be easily pre-allocated, so should be
recalled_item=zeros(nTrials,listlength);

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
    listCue = 1; % cues the list - at moment set to last LIST
    
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
    
    groupCue = 1; % start with first group; this will get incremented below
    needNewGroup = 1; % a flag to indicate whether we need a new group
        
    for outpos=1:3
        
        if needNewGroup
            
            eta_LC=1+randn(1,numGroups)*sigma_L; % Eq A3
            
            C_LC=eta_LC.*phi_l.^abs(listCue-lContext); %Eq A7 - control element for list
            C_NC=zeros(1,numGroups); % control element for group cue
            C_NC(groupCue)=eta_NC; % this shouldn't be list cue, but rather the
            % cue to a particular group (let's say the
            % first, but try setting it to last group)
            %C=zeros(1,numGroups); % don't need this, as it is overwritten in next line
            C=C_NC+C_LC; % list and group control elements added
            [max_value,currGroup]=max(C); % select most activated
            %currGroup refers to group we are currently trying to recall from
            
            % cue sequentially from first item
            % we find out how long the current group is, and then
            % construct scaled [0-1] within-grpup markers based on this
            % information
            currpContext = linspace(0,1,groupSize(currGroup));
            withinPos = 1;

            P_CG=1; % assume no effect of time
            
            needNewGroup = 0;
        end
        
        eta_gv=gamma.^(absP-1)+randn(1,listlength)*sigma_gp; %Eq A10
        v_GV = P_CG*eta_gv.*phi_g.^abs(currGroup-gContext); %Eq A11
        v_PV = phi_p.^abs(currpContext(withinPos)-pContext); %Eq A14
        
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
            recalled_item(t,outpos)=max_idx;
        else
            recalled_item(t,outpos)=0;
        end
        
        withinPos = withinPos + 1;
        
        if withinPos>groupSize(currGroup) % have we gone beyond the end of
                                            % current group?
            needNewGroup = 1;
            groupCue = groupCue + 1;
        end
    end
end

for i=1:listlength;
    prop(i)=numel(find(recalled_item==i))/nTrials;
end

plot(1:i,prop)

Av_v=mean(v);
Av_a=mean(a);
