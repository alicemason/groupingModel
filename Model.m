clear all
%Model Parameters
phi_p = 0.1;
phi_g=0.3;
gamma=0.97% rate of change of primacy gradient across groups
sigma_gp=0.02% 
rho=0.3;

% Experiment Setting

listlength=12;
j=1:listlength; % position of each item in the list
p=linspace(0,1,listlength);% information about within group position
groups=3;
g=linspace(0,1,groups); % information about group position
k=1:listlength/groups; % position of item within a group
item_cue=p(9);
group_cue=g;

% control elelements to lists (A1-3)
% control elelements to groups (A4-6)
% selecting control element to 1 (A7)

% set group context (g) - vector of pattern of scaled activations 1/sqr(numbergroups)
% group strengths

group_strength=phi_g.^abs(g-groups);
% group cues associated with current group context by updating
% to cue an item retrieve group context+position of item
%  % Associating items with group conetxt at encoding
%  % define strength of learning - primacy gradient in each group
 %eta_gv=gamma.^(k-1)+(rand*sigma_gp); %(A10)
%  % associate an item with it's group context
 %W_GV= eta_gv*(g*v)'  %(A9)
 
 % extent to which item is activated by group context
 %P_CG=1 %if group context is carried over from encoding
 % upsilon_GV=P_CG*sum(eta_gv*phi_g)
 % k indexes groups with which item been associated
 % i indexes context of group placed across all group context units -what's
 % this

 % within list position
item_strength= phi_p.^abs(p-item_cue)
  % associated item to within-group context
% W_PV=eta_g*(p*v)'; (A13)
% then scaled by noisy primacy gradient (A14) - 

% define v (A15)
% select item (A16) - competion between items

%% TASK
group_strength=phi_g.^abs(g-groups); % group strength relative to group n
item_strength= phi_p.^abs(p-item_cue) % item strength relative to first item in last group
% need to weightt each item relative to group strength 
% must be better way to do this?
group_1=item_strength(1:4)*group_strength(1)
group_2=item_strength(5:8)*group_strength(2)
group_3=item_strength(9:12)*group_strength(3)

% create a vector that tells you the group strength of each item
item_group_strength=[group_1,group_2,group_3]

plot(1:4,group_1,'+',5:8,group_2,'+',9:12,group_3,'+')

%% Simon's Version
% make group markers

groupSize = [3 3 3 3];

gContext = [];
pContext = [];

for gz=1:length(groupSize)
    gContext = [gContext repmat(gz,1,groupSize(gz))];
    pContext = [pContext linspace(0,1,groupSize(gz))];
end

v_GV = phi_g.^abs(gContext(10)-gContext);
v_PV = phi_p.^abs(pContext(10)-pContext);

v = rho*v_GV + (1-rho)*v_PV;

plot(v);


% multiply items by eta_gv to get primacy gradient?
% group_1=(item_strength(1:4)*group_strength(1)).*eta_gv
% group_2=(item_strength(5:8)*group_strength(2)).*eta_gv
% group_3=(item_strength(9:12)*group_strength(3)).*eta_gv


% item=1:listlength;
% plot(item,item_strength)
% xlabel('Item Strength')
% ylabel('P(recall)')




