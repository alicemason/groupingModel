clear all
%Model Parameters
phi_p = 0.1;
phi_g=0.3;
gamma=0.97% rate of change of primacy gradient across groups
sigma_gp=0.02% 
rho=0.3;

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