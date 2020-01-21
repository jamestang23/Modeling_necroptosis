function [simData,BASAL,v] = getSimData(id,Type)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Construct model information structure and simulate TNF stimulus
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Get default rate parameters, initial values, model time/timestep, and stimulus
[v.NP,v.IP] = getRateParams();
v.INIT = getInit(Type);
v.DT = id.DT;
v.SIM_TIME     = id.sim_time; % min of stimulation phase (phase 2)
v.TNF_DOSE = id.dose*(1.96e-4); % 1.96e-4uM= 1ng/mL TNF;
v.D_FLAG = 0;
v.Perturb=id.Perturb;
v.k=id.k;v.kd=id.kd;v.K=id.K;v.n=id.n;v.KO=id.KO;


% Update parameters using input modifications
if isfield(id,'inputPid')
    v.IP(id.inputPid) =id.inputP;
end
if isfield(id,'inputvPid')
    v.NP(id.inputvPid) =id.inputvP;
end

if isfield(id,'SET_IKK')
    if id.SET_IKK  
        v.IKK_input = id.IKK_input;
        v = generate_ikk(v); 
        v.SET_IKK = id.SET_IKK;
    end
else 
    v.SET_IKK = 0;
end


% Run ODE with no stimulus to calculate basal values
v.T_EQUILIBRATE   = -4000; 
%v.TNF_TIME = v.SIM_TIME+1; % Total length of stimulus
v.IKK_TOTAL  = sum(v.INIT.VALS(ismember(v.INIT.NAMES,{'IKK', 'IKK_off', 'IKK_i'})));
v = nfkbBasal(v);

% Run ODE with stimulus
v = nfkbSim(v,Type);


%% Return requested species dynamics
if numel(id.output)==1
    index =find(strcmp(v.INIT.NAMES,id.output));
    simData.induced =v.OUTPUT(:,index)';
elseif numel(id.output)>1
    for i=1:numel(id.output)
        index(i) =find(strcmp(v.INIT.NAMES,id.output{i}));
    end
    simData =v.OUTPUT(:,index)';
    BASAL =v.BASAL_VALUES(:,index)';
else
    error('func:getSim','Wrong output id number');
end