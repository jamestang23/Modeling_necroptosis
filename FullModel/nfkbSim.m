function v = nfkbSim(v,Type)
%% intro
% The simulation is run in 2 phases
%   1. An equilibrium phase.  This runs in increments of START_TIME
%       (by default 4000min) until there is no change > 1%
%   2. A stimulation phase. This runs from 0 to SIM_TIME minutes

%% Call the ODE function to reset persistent variables
nfkbOde([],[],[],v);

%% Run Phase 2
if v.SIM_TIME > 0
    v.PHASE         = 2;
    initial_values  = v.BASAL_VALUES;
    initial_values(strcmp('TNF',v.INIT.NAMES))  = v.TNF_DOSE; % TNF dose set
    
    
%     disp(initial_values);dd
    %initial_values(20)=0;initial_values(16)=0;%IkBdt,IkBd in order for threshold
    %initial_values(37)=0.05;%a20
    %Initial values for RIP1, RIP3, pMLKL
    initial_values(39)=0.03;initial_values(40)=0;initial_values(41)=0;
   if Type==3
    disp('2KO');
    %initial_values(6:9)=0.02*ones(1,4);%IkBb
    initial_values(16:19)=0.005*ones(1,4);%IkBd
   end
    if Type==1
    disp('WT');
    if v.Perturb(19)==1 % synthetic NFkB, make A20mRNA and A20 initial zero
        disp('synthetic')
        % initial_values(38)=0.5e-4;%A20mRNA
        % initial_values(37)=0.06;%a20
    end
    if v.Perturb(14)~=.6 % A20 constitutive scan
        disp('A20 constitutive scan');
        %initial_values(37)=0.06;%a20
    end
    
    if v.Perturb(18)==0 % A20-KO
%        initial_values(37)=0.06;%a20
    end
    
     %%initial_values(26)=0;%IKK % as suggested that it should be prior than NFkB
%     initial_values(38)=initial_values(38)*0.3;%A20mRNA
%     initial_values(37)=0.05;%a20
%     %initial_values(6:9)=0.000*ones(1,4);%IkBb
%     %initial_values(11:14)=0.001*ones(1,4);%IkBb
%     initial_values(16:19)=0.000*ones(1,4);%IkBd
    end
    
     if Type==2
    disp('RelA');
    initial_values(22)=0;%3e-3;%NFkBn
    %initial_values(21:22)=0.000*ones(1,2);%NFkB
    %initial_values(37)=0.0;%a20
    %initial_values(38)=0.0;%a20
    end
    
    if v.Perturb(19)==1
     initial_values(22)=v.Perturb(21);%delta_NFkBn*v.KO(2);
        %r2(Duration:end,22)=0;%delta_NFkBn*v.KO(2);
    end

    % TNF chronic
    [t2, r2] = ode15s('nfkbOde', [0 v.SIM_TIME],[initial_values],[],v);%Initial for RIP3 and pMLKL is 0
    
%     if v.Perturb(19)==1
%     Duration=round(size(r2,1)*v.Perturb(20)/24);  %synthetic NFkB
%         r2(1:Duration,22)=v.Perturb(21);%delta_NFkBn*v.KO(2);
%         r2(Duration:end,22)=0;%delta_NFkBn*v.KO(2);
%     end

% disp(v.SIM_TIME);
%     disp(size(t2));
   % disp(size(r2));%ddd
    
    % Interpolate Phase 2 Results
    %   This makes an array of length SIM_TIME+1 with per min values
    v.OUTPUT = interp1(t2,r2(:,:),0:v.DT:v.SIM_TIME);
end
end