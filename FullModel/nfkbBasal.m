function v=nfkbBasal(v,flag)
  
%% Call the ODE function to reset persistent variables
nfkbOde([],[],[],v);

%% Run Phase 1 (equilibrium phase)
v.PHASE     = 1;
static      = false;    % Set to 'true' after reaching equilibrium
count       = 1;        % Iteration Counter
threshold   = 1;        % Max % difference used by evaluate_phase1()
initial_values = v.INIT.VALS;
v.BASAL_VALUES = [];

while ~static % Iterate through Phase 1 until equilibrium is reached
    options = odeset('RelTol', 1e-10,'AbsTol',1e-10);
    

    [t1, r1] = ode15s('nfkbOde', [v.T_EQUILIBRATE 0], [initial_values],options,v);
%     plot(t1,r1(:,44));hold on ; 
    % Evaluate results and return true if at equilibrium
    
    if (evaluatePhase1(r1, threshold)) % values have converged
        v.T_EQUILIBRATE    = (count * v.T_EQUILIBRATE) - (count-1);
        static          = true;
        v.BASAL_VALUES  = r1(end,:);
        if v.D_FLAG; disp(['Equilibrium met at ' num2str(v.T_EQUILIBRATE), ' min']); end

    elseif count > 100 % values are not converging so stop running
        if v.D_FLAG; disp('===> max phase 1 steps reached (100 steps)'); end
        %disp(v.IP)
        static          = true;
        v.T_EQUILIBRATE    = (count * v.T_EQUILIBRATE) - (count-1);
    else % equilibrium not met, so run another round of Phase 1
        count           = count + 1;
        initial_values  = r1(end,:);
    end
end

%% subsubroutine static
    function static = evaluatePhase1(results,threshold)
        % threshold= percentage
        
        prev_values     = abs(results(1,:)'); %init value
        current_values  = abs(results(end,:)'); %end value
        
        max = current_values * (1 + (threshold/100));
        min = current_values * (1 - (threshold/100));
        
        
        for i = 1:length(prev_values) % examine all the state valuables (SVs).
            if ( (prev_values(i) >= max(i)) || (prev_values(i) <= min(i)) )
                if (current_values(i) > 1e-21)  % Resolves noise errors
                    static = false; % make sure all SVs statisfy the condition
                    return;
                end
            end
        end
        static = true;
    end
end