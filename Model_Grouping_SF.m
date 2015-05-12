
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


            for List_cue =1:3 % cues the list - at moment set to last LIST


            % list cue and group cue to context
            lContext = ones(1,numGroups);     % control element for list context (each group has a list context)
            eta_LC=1+randn(1,numGroups)*sigma_L; % Eq A3
            C_LC=eta_LC.*phi_l.^abs(List_cue-lContext); %Eq A7 - control element for list
            C_NC=zeros(1,numGroups); % control element for group cue
            C_NC(List_cue)=eta_NC; %cue to a particular group                       
            C=C_NC+C_LC; % list and group control elements added
            [max_value,max_idx]=max(C); % select most activated
            C(max_idx)=1; % set control element for most activared to 1

    Group_cue=max_idx;
    gContext=[];
    gContext = repmat(Group_cue,1,groupSize(Group_cue));
    % retrieve N items from current group assoicated with context
        for attempt=1:groupSize(Group_cue)
        % make group markers
            pContext = [];
            absP = [];
            pContext = linspace(0,1,groupSize(Group_cue));
            absP = [absP 1:groupSize(Group_cue)];
            item_probe = absP(attempt);

           % conext cued of the most activated group
            P_CG=1; % assume no effect of time

            eta_gv=gamma.^(absP-1)+randn(1,groupSize(Group_cue))*sigma_gp; %Eq A10
            v_GV = P_CG*eta_gv.*phi_g.^abs(Group_cue-gContext); %Eq A11
            v_PV = phi_p.^abs(pContext(item_probe)-pContext); %Eq A14

            % implemented primacy gradient
            v_PV=eta_gv.*v_PV;

            % sum group and item vectors to get activation of each item
            t_v = rho*v_GV + (1-rho)*v_PV; % Eq A15
            i=length(t_v);
            s= listlength-(sum(groupSize(Group_cue:end)))+attempt; % which item overall is being acyivated 
            v(t,s:(s+i-1))=t_v;
            % noisy retrieval Eq A16
            noise=randn(1,i)*sigma_v;
            a=(t_v+noise).*(1-r(length(t_v)));
            % activation of two highest items
            [max_value,max_idx] = max(a);
            a(max_idx) = NaN;
            second_max = max(a);
            a = max_value; % retruns max_value into a
            item_id= listlength-(sum(groupSize(Group_cue:end)))+max_idx;% which item across whole list has been retrived
            
            
            if (max_value-second_max)>theta
                recalled_item(t,s)=item_id;
            else
                recalled_item(t,s)=0;
            end

            end
            end
        end

    for i=1:listlength;
        prop(i)=numel(find(recalled_item==i))/nTrials;
    end

    plot(1:i,prop)

    Av_v=mean(v);
    Av_a=mean(a);
