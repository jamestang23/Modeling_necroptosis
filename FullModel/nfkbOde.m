function delta = nfkbOde(t,x,ode_options,v)
% Full TLR4-NFkB ODE model - stiff system, designed to be solved using ode15 (see this function's help file to set
% options). Simulation is run until convergence (with no stimulus)
%
% (earlier version history unknown)
% version 3.0: BT, 9.7.2014 - adds explicit CD14 binding reaction, tracks it as species.
%

% Initialize the iteration

% Set the persistant variables that handle the delay reactions
persistent DELAY_NFKB;
persistent DELAY_A20;

if isempty(t)
    % The sim function calls the ode function once without arguments to
    %   reset the persistent variables. Once done, it returns.
    DELAY_NFKB        = cell(3,1);       % Time, concentration, index
    DELAY_NFKB([1 2]) = {zeros(1000,1)}; % increase if getting errors
    DELAY_NFKB(3)     = {1}; % index (starts at 1)


    DELAY_A20         = cell(3,1);       % Time, concentration, index
    DELAY_A20([1 2])  = {zeros(1000,1)}; % increase if getting errors
    DELAY_A20(3)      = {1}; % index (starts at 1)
    
    return;
end




if v.Perturb(19)==1 %synthetic NFkB
if t/60<=v.Perturb(23)
x(22)=v.Perturb(22);
elseif v.Perturb(23)<t/60<v.Perturb(20)  
x(22)=v.Perturb(21);%delta_NFkBn*v.KO(2);
else
x(22)=0;
end
end
    
% Load Previous Concentrations
[a,b] = size(x);
delta = zeros(a,b);

% Species variables (can copy these from spreadsheet)
% NFkB module

IkBa	 = x(1);
IkBan	 = x(2);
IkBaNFkB	 = x(3);
IkBaNFkBn	 = x(4);
IkBat	 = x(5);
IkBb	 = x(6);
IkBbn	 = x(7);
IkBbNFkB	 = x(8);
IkBbNFkBn	 = x(9);
IkBbt	 = x(10);
IkBe	 = x(11);
IkBen	 = x(12);
IkBeNFkB	 = x(13);
IkBeNFkBn	 = x(14);
IkBet	 = x(15);
IkBd	 = x(16);
IkBdn	 = x(17);
IkBdNFkB	 = x(18);
IkBdNFkBn	 = x(19);
IkBdt	 = x(20);
NFkB	 = x(21);
NFkBn	 = x(22);
IKKK_off	 = x(23);
IKKK	 = x(24);
IKK_off	 = x(25);
IKK	 = x(26);
IKK_i	 = x(27);

% TNF module
TNF           = x(28);
tnfrm          = x(29);
TNFR           = x(30);
TNFRtnf        = x(31);
C1             = x(32);
C1_off         = x(33);
C1tnf          = x(34);
C1tnf_off      = x(35);
TTR            = x(36);
a20            = x(37);
a20t           = x(38);


RIP3           = x(39);
pMLKL           = x(40);
RIP1           = x(41);
%% Set iteration changes (delta) to zero
delta_IkBa        = 0;
delta_IkBan       = 0;
delta_IkBaNFkB    = 0;
delta_IkBaNFkBn   = 0;
delta_IkBat       = 0;
delta_IkBb        = 0;
delta_IkBbn       = 0;
delta_IkBbNFkB    = 0;
delta_IkBbNFkBn   = 0;
delta_IkBbt       = 0;
delta_IkBe        = 0;
delta_IkBen       = 0;
delta_IkBeNFkB    = 0;
delta_IkBeNFkBn   = 0;
delta_IkBet       = 0;
delta_IkBd        = 0;
delta_IkBdn       = 0;
delta_IkBdNFkB    = 0;
delta_IkBdNFkBn   = 0;
delta_IkBdt       = 0;
delta_NFkB        = 0;
delta_NFkBn       = 0;

%IKK module
delta_IKKK        = 0;
delta_IKKK_off    = 0;
delta_IKK         = 0;
delta_IKK_off     = 0;
delta_IKK_i       = 0;

% TNFR module
delta_TNF  = 0 ;
delta_tnfrm       = 0;
delta_TNFR        = 0;
delta_TNFRtnf     = 0;
delta_C1          = 0;
delta_C1_off      = 0;
delta_C1tnf       = 0;
delta_C1tnf_off   = 0;
delta_TTR         = 0;
delta_a20         = 0;
delta_a20t        = 0;


%% init delay
% Calculate Inducible Transcription and A20 Translation Delays
% By default, set the delays to the current values (no delay)
if v.KO(2)==0
    NFkBn=0;%3e-3;
    NFkB=0;
    x(21)=0;
    x(22)=0;%3e-3;
end
delayed_nfkbn_a   = NFkBn;
delayed_nfkbn_b   = NFkBn;
delayed_nfkbn_e   = NFkBn;
delayed_nfkbn_d = NFkBn;
delayed_nfkbn_a20 = NFkBn;
delayed_a20t      = a20t ;


if v.PHASE == 2 % no delay in phase 1

    % Create a cache of previous NFkBn values for txn delays.
    %   The find command is needed to prevent duplicate values that
    %   break the interpolation function
    if (isempty(find(DELAY_NFKB{1}(1:DELAY_NFKB{3})== t,1)))
        DELAY_NFKB{1}(DELAY_NFKB{3}) = t;
        DELAY_NFKB{2}(DELAY_NFKB{3}) = NFkBn;
        DELAY_NFKB{3}= DELAY_NFKB{3} + 1; % iterate the index + 1
    end

    if v.NP(10) > 0  % IkBa inducible txn delay
        if t > v.NP(10)
            delayed_nfkbn_a = ...
                interp1(DELAY_NFKB{1}(1:DELAY_NFKB{3}-1),...
                DELAY_NFKB{2}(1:DELAY_NFKB{3}-1),t);%-v.NP(10));
            % because iterate inde
        else
            delayed_nfkbn_a = ...
                DELAY_NFKB{2}(1); % 1st index is the basal state value
        end
    end

    if v.NP(11) > 0  % IkBb inducible txn delay
        if t > v.NP(11)
            delayed_nfkbn_b = ...
                interp1(DELAY_NFKB{1}(1:DELAY_NFKB{3}-1),...
                DELAY_NFKB{2}(1:DELAY_NFKB{3}-1),t);%-v.NP(11));
        else
            delayed_nfkbn_b = ...
                DELAY_NFKB{2}(1); % 1st index is the basal state value
        end
    end

    if v.NP(12) > 0  % IkBe inducible txn delay
        if t > v.NP(12)
            if v.NP(12) == v.NP(11)  % saves cpu time if IkBb=IkBe delay
                delayed_nfkbn_e = delayed_nfkbn_b;
            else
                delayed_nfkbn_e = ...
                    interp1(DELAY_NFKB{1}(1:DELAY_NFKB{3}-1),...
                    DELAY_NFKB{2}(1:DELAY_NFKB{3}-1),t);%-v.NP(12));
            end
        else
            delayed_nfkbn_e = ...
                DELAY_NFKB{2}(1); % 1st index is the basal state value
        end
    end

    if v.NP(75) >0 % have tried to delete IkBd delay %> 0  % IkBd inducible txn delay
        if t > v.NP(75)
            delayed_nfkbn_d = ...
               ... interp1(DELAY_NFKB{1}(1:DELAY_NFKB{3}-1),...
                ...        DELAY_NFKB{2}(1:DELAY_NFKB{3}-1),t-v.NP(75));
                ...        %change delay mode: previouis delay by following earlier NFkB
                interp1(DELAY_NFKB{1}(1:DELAY_NFKB{3}-1),...
                        DELAY_NFKB{2}(1:DELAY_NFKB{3}-1),t);%-v.NP(75));
        else
            delayed_nfkbn_d = ...
                DELAY_NFKB{2}(1); % 1st index is the basal state value
        end
    end
    
    
    if v.NP(66) > 0  % A20 inducible txn delay
        if t > v.NP(66) 
            delayed_nfkbn_a20 = ...
                interp1(DELAY_NFKB{1}(1:DELAY_NFKB{3}-1),...
                        DELAY_NFKB{2}(1:DELAY_NFKB{3}-1),t);%-v.NP(66));        
        else
            delayed_nfkbn_a20 = ...
                DELAY_NFKB{2}(1); % 1st index is the basal state value  
        end
        disp('A20 delay');
    end
         
    % This code enables the A20 translational delay
    %   It does so by caching A20 mRNA levels
    if v.NP([64,68,69]) > 0
        % Create a cache of previous A20 mRNA values. 
        if(isempty(find(DELAY_A20{1}(1:DELAY_A20{3}) == t,1)))
            DELAY_A20{1}(DELAY_A20{3}) = t;
            DELAY_A20{2}(DELAY_A20{3}) = delayed_a20t;
            DELAY_A20{3}= DELAY_A20{3} + 1; % iterate the index + 1
        end
        
        % Calculate the A20 mRNA value at time = t-delay
        if t > v.NP(69)
            delayed_a20t = ...
                interp1(DELAY_A20{1}(1:DELAY_A20{3}-1),...
                        DELAY_A20{2}(1:DELAY_A20{3}-1),t);%-v.NP(69));
        else
            delayed_a20t = DELAY_A20{2}(1);
        end
    end
        
end

%% set IKK flux
if v.PHASE == 1
    IKK_flux = IKK/v.IKK_TOTAL;
    IKK1_flux = .01; % basal
elseif v.SET_IKK ==0
    IKK_flux = IKK/v.IKK_TOTAL;
    IKK1_flux = .01; % basal
else
    [~,b] = size(v.IKK_CURVE);
    t_range = 0:1:(b-1);
    IKK_flux = interp1(t_range, v.IKK_CURVE, t);
    
    [~,b] = size(v.IKK_CURVE1);
    t_range = 0:1:(b-1);
    IKK1_flux = interp1(t_range, v.IKK_CURVE1, t);
end




%% Calculate Iteration Fluxes
%-- NFkB Activation Module
flux_rsu_a         = v.NP(1) ;
flux_rsu_b         = v.NP(2) *v.Perturb(12);
flux_rsu_e         = v.NP(3) ;
flux_rsu_d         = v.NP(72)*v.Perturb(13) ;

flux_rsr_an        = v.NP(4)   * (delayed_nfkbn_a ^ v.NP(7)) ;
flux_rsr_bn        = v.NP(5)   * (delayed_nfkbn_b ^ v.NP(8))*v.Perturb(9) ;
flux_rsr_en        = v.NP(6)   * (delayed_nfkbn_e ^ v.NP(9)) ;
flux_rsr_dn        = v.NP(73)  * (delayed_nfkbn_d ^ v.NP(74))*v.Perturb(7) ;

flux_rd_a          = v.NP(13)  * IkBat	;
flux_rd_b          = v.NP(14)  * IkBbt	;
flux_rd_e          = v.NP(15)  * IkBet	;
flux_rd_d          = v.NP(76)  * IkBdt	;

flux_ps_c_a        = v.NP(16)	* IkBat	;
flux_ps_c_b        = v.NP(17)	* IkBbt	;
flux_ps_c_e        = v.NP(18)	* IkBet	;        
flux_ps_c_d        = v.NP(77)	* IkBdt	 ;        

flux_in_a          = v.NP(19)	* IkBa	;
flux_in_b          = v.NP(20)	* IkBb	;
flux_in_e          = v.NP(21)	* IkBe	;
flux_in_d          = v.NP(78)	* IkBd;
flux_in_n          = v.NP(22)	* NFkB	;

flux_ex_a          = v.NP(23)	* IkBan	;
flux_ex_b          = v.NP(24)	* IkBbn	;
flux_ex_e          = v.NP(25)	* IkBen	;
flux_ex_d          = v.NP(79)	* IkBdn	;
flux_ex_n          = v.NP(26)	* NFkBn	;

flux_in_2an        = v.NP(27)	* IkBaNFkB	;
flux_in_2bn        = v.NP(28)	* IkBbNFkB	;
flux_in_2en        = v.NP(29)	* IkBeNFkB	;
flux_in_2dn        = v.NP(80)	* IkBdNFkB	;

flux_ex_2an        = v.NP(30)	* IkBaNFkBn	;
flux_ex_2bn        = v.NP(31)	* IkBbNFkBn	;
flux_ex_2en        = v.NP(32)	* IkBeNFkBn	;
flux_ex_2dn        = v.NP(81)	* IkBdNFkBn	;

flux_pd_c_a        = v.NP(33)	* IkBa	;
flux_pd_c_b        = v.NP(34)	* IkBb	;
flux_pd_c_e        = v.NP(35)	* IkBe	;
flux_pd_c_d        = v.NP(82)	* IkBd;

flux_pd_n_a        = v.NP(36)	* IkBan	;
flux_pd_n_b        = v.NP(37)	* IkBbn	;
flux_pd_n_e        = v.NP(38)	* IkBen	;
flux_pd_n_d        = v.NP(83)	* IkBdn	;

flux_pd_c_2an      = v.NP(39)	* IkBaNFkB	;
flux_pd_c_2bn      = v.NP(40)	* IkBbNFkB	;
flux_pd_c_2en      = v.NP(41)	* IkBeNFkB	;
flux_pd_c_2dn      = v.NP(84)	* IkBdNFkB	;

flux_pd_n_2an      = v.NP(42)	* IkBaNFkBn	;
flux_pd_n_2bn      = v.NP(43)	* IkBbNFkBn	;
flux_pd_n_2en      = v.NP(44)	* IkBeNFkBn	;
flux_pd_n_2dn      = v.NP(85)	* IkBdNFkBn	;

flux_a_c_an        = v.NP(45)  * IkBa	* NFkB	;
flux_a_c_bn        = v.NP(46)  * IkBb	* NFkB	;
flux_a_c_en        = v.NP(47)  * IkBe	* NFkB	;
flux_a_c_dn        = v.NP(86)  * IkBd	* NFkB	;

flux_a_n_an        = v.NP(48)	* IkBan	* NFkBn	;
flux_a_n_bn        = v.NP(49)	* IkBbn	* NFkBn	;
flux_a_n_en        = v.NP(50)  * IkBen	* NFkBn	;
flux_a_n_dn        = v.NP(87)  * IkBdn	* NFkBn	;

flux_d_c_an        = v.NP(51)	* IkBaNFkB	;
flux_d_c_bn        = v.NP(52)	* IkBbNFkB	;
flux_d_c_en        = v.NP(53)	* IkBeNFkB	;
flux_d_c_dn        = v.NP(88)	* IkBdNFkB	;

flux_d_n_an        = v.NP(54)	* IkBaNFkBn	;
flux_d_n_bn        = v.NP(55)	* IkBbNFkBn	;
flux_d_n_en        = v.NP(56)	* IkBeNFkBn	;
flux_d_n_dn        = v.NP(89)	* IkBdNFkBn	;

% IKK Mediated IkB Degradation (free and bound)
flux_ph_c_a        = v.NP(57)   * IkBa * IKK_flux ;
flux_ph_c_b        = v.NP(58)	* IkBb * IKK_flux ;
flux_ph_c_e        = v.NP(59)	* IkBe * IKK_flux ;
flux_ph_c_d        = v.NP(90)	* IkBd * IKK1_flux ;

flux_ph_c_an       = v.NP(60)   * IkBaNFkB * IKK_flux ;
flux_ph_c_bn       = v.NP(61)	* IkBbNFkB * IKK_flux ;
flux_ph_c_en       = v.NP(62)	* IkBeNFkB * IKK_flux ;
flux_ph_c_dn       = v.NP(91)	* IkBdNFkB * IKK1_flux ;


% A20 Fluxes
flux_rsu_a20       = v.NP(63)*v.Perturb(14) ; %constituive %/4 because to shift the consistuitive scan for WT 2^0 fold
flux_rsr_a20       = v.NP(64)  * (delayed_nfkbn_a20 ^ v.NP(65))  *v.Perturb(5).*v.Perturb(10);%%%%%distribute a20 systhesis %Induction
flux_rd_a20        = v.NP(67)  * a20t.*v.Perturb(11);%%%%%%%%%distribute a20 mRNA degradation 
flux_ps_c_a20      = v.NP(68)  * delayed_a20t*v.Perturb(6) ;
flux_pd_c_a20      = v.NP(70)  * a20*v.Perturb(4)  ;    

if v.Perturb(19)==1 % synthetic NFkB disables inducible A20 at different time
    if( t >(v.Perturb(20))*60)  % disables inducible A20 txn
    flux_rsr_a20 = 0; %min
    end

else
    if( t > v.NP(71))  % disables inducible A20 txn to match exp data
    flux_rsr_a20 = 0;
    end
end




% TNF module 
% Adapted from Werner et al. 2008
%TNFR part
flux_pd_m_tnf      = v.IP(1)    * TNF*v.Perturb(1); %Perturb % pd_m_tnf 45' half life of exogenous TNF
if v.Perturb(16)==1 %if TNF pulse experiment
if t/60>v.Perturb(17)  %pulse-TNF
    flux_pd_m_tnf      = v.IP(1)    * TNF*v.Perturb(1)*1500;% *150; %Perturb % pd_m_tnf 45' half life of exogenous TNF    
end
end

flux_syn_tnfrm     = v.IP(2);  % tnfrm synthesis (txn, tsl, localization)
flux_pd_tnfrm      = v.IP(3)    * tnfrm;% tnfrm --> deg  -- 120' halflife
flux_a_tnfrm       = v.IP(4)    * tnfrm;% 3tnfrm --> TNFR
flux_d_TNFR        = v.IP(5)    * TNFR;% TNFR   --> 3tnfrm
flux_i_TNFR        = v.IP(6)    * TNFR;% TNFR internalization -- 30' halflife
flux_a_C1_off      = v.IP(7)    * TNFR* TTR;% TNFR + TTR --> C1_off
flux_d_C1_off      = v.IP(8)    * C1_off;% C1_off --> TNFR + TTR
flux_a_C1          = v.IP(9)    * C1_off*v.Perturb(2);% C1_off --> C1
flux_C1_off        = v.IP(10)   * C1;% C1     --> C1_off
flux_C1_A20        = v.IP(11)   * C1 * a20*v.Perturb(15);% C1     --> C1_off (A20 Mediated)
flux_d_C1          = v.IP(12)   * C1; %UPDATE % C1     --> TNFR + TTR
flux_i_C1_off      = v.IP(13)   * C1_off*v.Perturb(3); % C1_off internalization
flux_i_C1          = v.IP(14)   * C1; % C1 internalization
flux_a_tnfrmtnf    = v.IP(15)   * tnfrm * TNF; % 3tnfrm + tnf --> TNFRtnf
flux_a_TNFRtnf     = v.IP(16)   * TNFR * TNF; % TNFR + tnf --> TNFRtnf
flux_d_TNFRtnf     = v.IP(17)   * TNFRtnf; % TNFRtnf   --> TNFR + tnf
flux_i_TNFRtnf     = v.IP(18)   * TNFRtnf;  % TNFRtnf internalization -- 30' halflife
flux_a_C1tnf_off   = v.IP(19)   * TNFRtnf * TTR; % TNFRtnf + TTR --> C1tnf_off
flux_d_C1tnf_off   = v.IP(20)   * C1tnf_off; % C1tnf_off --> TNFRtnf + TTR
flux_a_C1tnf       = v.IP(21)   * C1tnf_off; % C1tnf_off --> C1tnf
flux_C1tnf_off     = v.IP(22)   * C1tnf; % C1tnf     --> C1tnf_off
flux_C1tnf_A20     = v.IP(23)   * C1tnf * a20*v.Perturb(15);  % C1tnf     --> C1tnf_off (A20 Mediated)
flux_d_C1tnf       = v.IP(24)   * C1tnf; % C1tnf --> TNFRtnf + TTR
flux_i_C1tnf_off   = v.IP(25)   * C1tnf_off;  % C1tnf_off internalization
flux_i_C1tnf       = v.IP(26)   * C1tnf; % C1tnf internalization
flux_d_tnf_C1_off  = v.IP(27)   * C1tnf_off; % C1tnf_off --> C1_off + tnf
flux_a_tnf_C1_off  = v.IP(28)   * C1_off * TNF; % C1_off + tnf --> C1tnf_off
flux_d_tnf_C1      = v.IP(29)   * C1tnf; % C1tnf    --> C1 + tnf
flux_a_tnf_C1      = v.IP(30)   * C1 * TNF; % C1 + tnf --> C1tnf
flux_IKKK_on_C1    = v.IP(31)   * IKKK_off * C1; % IKKK_off --> IKKK (C1 mediated)?500
flux_IKKK_on_C1tnf = v.IP(32)   * IKKK_off * C1tnf;% IKKK_off --> IKKK (C1tnf mediated)
flux_IKKK_on       = v.IP(33)   * IKKK_off;
flux_IKKK_off      = v.IP(34)   * IKKK;


% IKK cycle 
flux_IKK_on_IKKK   = v.IP(35)   * IKK_off * IKKK;
flux_IKK_on        = v.IP(36)   * IKK_off;
flux_IKK_off       = v.IP(37)   * IKK;
flux_IKK_off_i     = v.IP(38)   * IKK;
flux_IKKi          = v.IP(39)   * IKK_i;


%% Add Fluxes to component concentrations, then Save

% IkB Transcription
delta_IkBat     = delta_IkBat + flux_rsu_a;
delta_IkBbt     = delta_IkBbt + flux_rsu_b;
delta_IkBet     = delta_IkBet + flux_rsu_e;
delta_IkBdt     = delta_IkBdt + flux_rsu_d;

delta_IkBat     = delta_IkBat + flux_rsr_an;
delta_IkBbt     = delta_IkBbt + flux_rsr_bn;
delta_IkBet     = delta_IkBet + flux_rsr_en;
delta_IkBdt     = delta_IkBdt + flux_rsr_dn;

% IkB Transcrv.IPt Degradation
delta_IkBat     = delta_IkBat - flux_rd_a;
delta_IkBbt     = delta_IkBbt - flux_rd_b;
delta_IkBet     = delta_IkBet - flux_rd_e;
delta_IkBdt     = (delta_IkBdt - flux_rd_d);%*(v.Perturb(8)-IkBdt); %Perturb: set a threshold on IkBd mRNA synthesis


% IkB Translation
delta_IkBa      = delta_IkBa  + flux_ps_c_a;
delta_IkBb      = delta_IkBb  + flux_ps_c_b;
delta_IkBe      = delta_IkBe  + flux_ps_c_e;
delta_IkBd      = delta_IkBd  + flux_ps_c_d;%perturb

% IkB:NFkB Shuttling (Free and Bound)
delta_IkBa      = delta_IkBa  - flux_in_a;
delta_IkBan     = delta_IkBan + flux_in_a;

delta_IkBb      = delta_IkBb  - flux_in_b;
delta_IkBbn     = delta_IkBbn + flux_in_b;

delta_IkBe      = delta_IkBe  - flux_in_e;
delta_IkBen     = delta_IkBen + flux_in_e;

delta_IkBd      = delta_IkBd  - flux_in_d;
delta_IkBdn     = delta_IkBdn + flux_in_d;

delta_NFkB      = delta_NFkB  - flux_in_n;
delta_NFkBn     = delta_NFkBn + flux_in_n;

delta_IkBan     = delta_IkBan - flux_ex_a;
delta_IkBa      = delta_IkBa  + flux_ex_a;

delta_IkBbn     = delta_IkBbn - flux_ex_b;
delta_IkBb      = delta_IkBb  + flux_ex_b;

delta_IkBen     = delta_IkBen - flux_ex_e;
delta_IkBe      = delta_IkBe  + flux_ex_e;

delta_IkBdn     = delta_IkBdn - flux_ex_d;
delta_IkBd      = delta_IkBd  + flux_ex_d;

delta_NFkBn     = delta_NFkBn - flux_ex_n;
delta_NFkB      = delta_NFkB  + flux_ex_n;

delta_IkBaNFkB  = delta_IkBaNFkB  - flux_in_2an;
delta_IkBaNFkBn = delta_IkBaNFkBn + flux_in_2an;

delta_IkBbNFkB  = delta_IkBbNFkB  - flux_in_2bn;
delta_IkBbNFkBn = delta_IkBbNFkBn + flux_in_2bn;

delta_IkBeNFkB  = delta_IkBeNFkB  - flux_in_2en;
delta_IkBeNFkBn = delta_IkBeNFkBn + flux_in_2en;

delta_IkBdNFkB  = delta_IkBdNFkB  - flux_in_2dn;
delta_IkBdNFkBn = delta_IkBdNFkBn + flux_in_2dn;

delta_IkBaNFkBn = delta_IkBaNFkBn - flux_ex_2an;
delta_IkBaNFkB  = delta_IkBaNFkB  + flux_ex_2an;

delta_IkBbNFkBn = delta_IkBbNFkBn - flux_ex_2bn;
delta_IkBbNFkB  = delta_IkBbNFkB  + flux_ex_2bn;

delta_IkBeNFkBn = delta_IkBeNFkBn - flux_ex_2en;
delta_IkBeNFkB  = delta_IkBeNFkB  + flux_ex_2en;

delta_IkBdNFkBn = delta_IkBdNFkBn - flux_ex_2dn;
delta_IkBdNFkB  = delta_IkBdNFkB  + flux_ex_2dn;

% IkB:NFkB Association (Cytoplasm and Nucleus)
delta_IkBa      = delta_IkBa - flux_a_c_an;
delta_NFkB      = delta_NFkB - flux_a_c_an;
delta_IkBaNFkB  = delta_IkBaNFkB + flux_a_c_an;

delta_IkBb      = delta_IkBb - flux_a_c_bn;
delta_NFkB      = delta_NFkB - flux_a_c_bn;
delta_IkBbNFkB  = delta_IkBbNFkB + flux_a_c_bn;

delta_IkBe      = delta_IkBe - flux_a_c_en;
delta_NFkB      = delta_NFkB - flux_a_c_en;
delta_IkBeNFkB  = delta_IkBeNFkB + flux_a_c_en;

delta_IkBd      = delta_IkBd - flux_a_c_dn;
delta_NFkB      = delta_NFkB - flux_a_c_dn;
delta_IkBdNFkB  = delta_IkBdNFkB + flux_a_c_dn;

delta_IkBan     = delta_IkBan - flux_a_n_an;
delta_NFkBn     = delta_NFkBn - flux_a_n_an;
delta_IkBaNFkBn = delta_IkBaNFkBn + flux_a_n_an;

delta_IkBbn     = delta_IkBbn - flux_a_n_bn;
delta_NFkBn     = delta_NFkBn - flux_a_n_bn;
delta_IkBbNFkBn = delta_IkBbNFkBn + flux_a_n_bn;

delta_IkBen     = delta_IkBen - flux_a_n_en;
delta_NFkBn     = delta_NFkBn - flux_a_n_en;
delta_IkBeNFkBn = delta_IkBeNFkBn + flux_a_n_en;

delta_IkBdn     = delta_IkBdn - flux_a_n_dn;
delta_NFkBn     = delta_NFkBn - flux_a_n_dn;
delta_IkBdNFkBn = delta_IkBdNFkBn + flux_a_n_dn;

% IkB:NFkB Dissociation (Cytoplasm and Nucleus)
delta_IkBaNFkB  = delta_IkBaNFkB - flux_d_c_an;
delta_IkBa      = delta_IkBa + flux_d_c_an;
delta_NFkB      = delta_NFkB + flux_d_c_an;

delta_IkBbNFkB  = delta_IkBbNFkB - flux_d_c_bn;
delta_IkBb      = delta_IkBb + flux_d_c_bn;
delta_NFkB      = delta_NFkB + flux_d_c_bn;

delta_IkBeNFkB  = delta_IkBeNFkB - flux_d_c_en;
delta_IkBe      = delta_IkBe + flux_d_c_en;
delta_NFkB      = delta_NFkB + flux_d_c_en;

delta_IkBdNFkB  = delta_IkBdNFkB - flux_d_c_dn;
delta_IkBd      = delta_IkBd + flux_d_c_dn;
delta_NFkB      = delta_NFkB + flux_d_c_dn;

delta_IkBaNFkBn = delta_IkBaNFkBn - flux_d_n_an;
delta_IkBan     = delta_IkBan + flux_d_n_an;
delta_NFkBn     = delta_NFkBn + flux_d_n_an;

delta_IkBbNFkBn = delta_IkBbNFkBn - flux_d_n_bn;
delta_IkBbn     = delta_IkBbn + flux_d_n_bn;
delta_NFkBn     = delta_NFkBn + flux_d_n_bn;

delta_IkBeNFkBn = delta_IkBeNFkBn - flux_d_n_en;
delta_IkBen     = delta_IkBen + flux_d_n_en;
delta_NFkBn     = delta_NFkBn + flux_d_n_en;

delta_IkBdNFkBn = delta_IkBdNFkBn - flux_d_n_dn;
delta_IkBdn     = delta_IkBdn + flux_d_n_dn;
delta_NFkBn     = delta_NFkBn + flux_d_n_dn;

% Free IkB Degradation (Cytoplasm and Nucleus)
delta_IkBa      = delta_IkBa  - flux_pd_c_a;
delta_IkBb      = delta_IkBb  - flux_pd_c_b;
delta_IkBe      = delta_IkBe  - flux_pd_c_e;
delta_IkBd      = delta_IkBd  - flux_pd_c_d;

delta_IkBan     = delta_IkBan - flux_pd_n_a;
delta_IkBbn     = delta_IkBbn - flux_pd_n_b;
delta_IkBen     = delta_IkBen - flux_pd_n_e;
delta_IkBdn     = delta_IkBdn - flux_pd_n_d;

% IkB:NFkB Degradation (Cytoplasm and Nucleus)
delta_IkBaNFkB  = delta_IkBaNFkB - flux_pd_c_2an;
delta_NFkB      = delta_NFkB + flux_pd_c_2an;

delta_IkBbNFkB  = delta_IkBbNFkB - flux_pd_c_2bn;
delta_NFkB      = delta_NFkB + flux_pd_c_2bn;

delta_IkBeNFkB  = delta_IkBeNFkB - flux_pd_c_2en;
delta_NFkB      = delta_NFkB + flux_pd_c_2en;

delta_IkBdNFkB  = delta_IkBdNFkB - flux_pd_c_2dn;
delta_NFkB      = delta_NFkB + flux_pd_c_2dn;

delta_IkBaNFkBn = delta_IkBaNFkBn - flux_pd_n_2an;
delta_NFkBn     = delta_NFkBn + flux_pd_n_2an;

delta_IkBbNFkBn = delta_IkBbNFkBn - flux_pd_n_2bn;
delta_NFkBn     = delta_NFkBn + flux_pd_n_2bn;

delta_IkBeNFkBn = delta_IkBeNFkBn - flux_pd_n_2en;
delta_NFkBn     = delta_NFkBn + flux_pd_n_2en;

delta_IkBdNFkBn = delta_IkBdNFkBn - flux_pd_n_2dn;
delta_NFkBn     = delta_NFkBn + flux_pd_n_2dn;

% IKK Mediated IkB Degradation
delta_IkBa      = delta_IkBa - flux_ph_c_a;
delta_IkBb      = delta_IkBb - flux_ph_c_b;
delta_IkBe      = delta_IkBe - flux_ph_c_e;
delta_IkBd      = (delta_IkBd - flux_ph_c_d);%*(v.Perturb(9)-IkBd);

% IKK Mediated IkB:NFkB Degradation
delta_IkBaNFkB  = delta_IkBaNFkB - flux_ph_c_an;
delta_NFkB      = delta_NFkB + flux_ph_c_an;

delta_IkBbNFkB  = delta_IkBbNFkB - flux_ph_c_bn;
delta_NFkB      = delta_NFkB + flux_ph_c_bn;

delta_IkBeNFkB  = delta_IkBeNFkB - flux_ph_c_en;
delta_NFkB      = delta_NFkB + flux_ph_c_en;

delta_IkBdNFkB  = delta_IkBdNFkB - flux_ph_c_dn;
delta_NFkB      = delta_NFkB + flux_ph_c_dn;

%---------------TNR module---------------
% A20 Expression 
% Transcription
delta_a20t      = delta_a20t + flux_rsu_a20;
delta_a20t      = delta_a20t + flux_rsr_a20;
delta_a20t      = delta_a20t - flux_rd_a20;

% Translation
delta_a20       = delta_a20  + flux_ps_c_a20;

% Degradation
delta_a20       = delta_a20 - flux_pd_c_a20;            


% Upstream Pathway  [TNF] --> IKK activation
% TNF-independent Activation
% tnfrm trimerization 3tnfrm<--> TNFR
delta_tnfrm     = delta_tnfrm- 3*flux_a_tnfrm;
delta_TNFR      = delta_TNFR +   flux_a_tnfrm;

delta_tnfrm     = delta_tnfrm+ 3*flux_d_TNFR;
delta_TNFR      = delta_TNFR -   flux_d_TNFR;

% TNF Receptor Metabolism
delta_TNFR     = delta_TNFR     - flux_i_TNFR;
delta_TNFRtnf  = delta_TNFRtnf  - flux_i_TNFRtnf;

delta_tnfrm    = delta_tnfrm    - flux_pd_tnfrm;
delta_tnfrm    = delta_tnfrm    + flux_syn_tnfrm;

% Complex I Generation  TTR + TNFR <--> C1_off, C1
delta_C1_off    = delta_C1_off  + flux_a_C1_off;
delta_TTR       = delta_TTR     - flux_a_C1_off;
delta_TNFR      = delta_TNFR    - flux_a_C1_off;

delta_C1_off    = delta_C1_off  - flux_d_C1_off;
delta_TTR       = delta_TTR     + flux_d_C1_off;
delta_TNFR      = delta_TNFR    + flux_d_C1_off;

delta_C1        = delta_C1      - flux_d_C1;
delta_TTR       = delta_TTR     + flux_d_C1;
delta_TNFR      = delta_TNFR    + flux_d_C1;

% C1_off  <--> C1
delta_C1_off    = delta_C1_off  - flux_a_C1;
delta_C1        = delta_C1      + flux_a_C1;

delta_C1_off    = delta_C1_off  + flux_C1_A20;
delta_C1        = delta_C1      - flux_C1_A20;

delta_C1_off    = delta_C1_off  + flux_C1_off;
delta_C1        = delta_C1      - flux_C1_off;

% C1 and C1_off internalization (treated as deg)
delta_C1        = delta_C1      - flux_i_C1;
delta_C1_off    = delta_C1_off  - flux_i_C1_off;

% TNF-Dependent Activation
% TNF degradation
delta_TNF       = delta_TNF     - flux_pd_m_tnf;

% TNF:TNFR binding
delta_TNFRtnf      = delta_TNFRtnf    +   flux_a_tnfrmtnf;
delta_tnfrm     = delta_tnfrm   - 3*flux_a_tnfrmtnf;
delta_TNF       = delta_TNF     -   flux_a_tnfrmtnf;

delta_TNFR      = delta_TNFR    - flux_a_TNFRtnf;
delta_TNF       = delta_TNF     - flux_a_TNFRtnf;
delta_TNFRtnf   = delta_TNFRtnf + flux_a_TNFRtnf;

delta_TNFR      = delta_TNFR    + flux_d_TNFRtnf;
delta_TNF       = delta_TNF     + flux_d_TNFRtnf;
delta_TNFRtnf   = delta_TNFRtnf - flux_d_TNFRtnf;

% Complex I Generation  TTR + TNFRtnf <--> C1tnf_off, C1tnf
delta_C1tnf_off = delta_C1tnf_off   + flux_a_C1tnf_off;
delta_TTR       = delta_TTR         - flux_a_C1tnf_off;
delta_TNFRtnf   = delta_TNFRtnf     - flux_a_C1tnf_off;

delta_C1tnf_off = delta_C1tnf_off   - flux_d_C1tnf_off;
delta_TTR       = delta_TTR         + flux_d_C1tnf_off;
delta_TNFRtnf   = delta_TNFRtnf     + flux_d_C1tnf_off;

delta_C1tnf     = delta_C1tnf       - flux_d_C1tnf;
delta_TTR       = delta_TTR         + flux_d_C1tnf;
delta_TNFRtnf   = delta_TNFRtnf     + flux_d_C1tnf;

% C1tnf_off <--> C1tnf
delta_C1tnf_off = delta_C1tnf_off   - flux_a_C1tnf;
delta_C1tnf     = delta_C1tnf       + flux_a_C1tnf;

delta_C1tnf_off = delta_C1tnf_off   + flux_C1tnf_off;
delta_C1tnf     = delta_C1tnf       - flux_C1tnf_off;

delta_C1tnf_off = delta_C1tnf_off   + flux_C1tnf_A20;
delta_C1tnf     = delta_C1tnf       - flux_C1tnf_A20;

% C1tnf and C1tnf_off internalization (treated as deg)
delta_C1tnf_off = delta_C1tnf_off   - flux_i_C1tnf_off;
delta_C1tnf     = delta_C1tnf       - flux_i_C1tnf;

% C1tnf_off,C1tnf <--> C1_off,C1 + TNF
delta_C1tnf_off = delta_C1tnf_off   - flux_d_tnf_C1_off;
delta_C1_off    = delta_C1_off      + flux_d_tnf_C1_off;
delta_TNF       = delta_TNF         + flux_d_tnf_C1_off;

delta_C1tnf_off = delta_C1tnf_off   + flux_a_tnf_C1_off;
delta_C1_off    = delta_C1_off      - flux_a_tnf_C1_off;
delta_TNF       = delta_TNF         - flux_a_tnf_C1_off;

delta_C1tnf     = delta_C1tnf       - flux_d_tnf_C1;
delta_C1        = delta_C1          + flux_d_tnf_C1;
delta_TNF       = delta_TNF         + flux_d_tnf_C1;

delta_C1tnf     = delta_C1tnf       + flux_a_tnf_C1;
delta_C1        = delta_C1          - flux_a_tnf_C1;
delta_TNF       = delta_TNF         - flux_a_tnf_C1;

% IKKK Regulation
delta_IKKK      = delta_IKKK     + flux_IKKK_on;
delta_IKKK_off  = delta_IKKK_off - flux_IKKK_on;

delta_IKKK      = delta_IKKK     + flux_IKKK_on_C1tnf;
delta_IKKK_off  = delta_IKKK_off - flux_IKKK_on_C1tnf;

delta_IKKK      = delta_IKKK     + flux_IKKK_on_C1;
delta_IKKK_off  = delta_IKKK_off - flux_IKKK_on_C1;

delta_IKKK_off  = delta_IKKK_off + flux_IKKK_off;
delta_IKKK      = delta_IKKK     - flux_IKKK_off;

% IKK Regulation
delta_IKK       = delta_IKK     + flux_IKK_on;
delta_IKK_off   = delta_IKK_off - flux_IKK_on;

delta_IKK       = delta_IKK     + flux_IKK_on_IKKK;
delta_IKK_off   = delta_IKK_off - flux_IKK_on_IKKK;

delta_IKK       = delta_IKK     - flux_IKK_off;
delta_IKK_off   = delta_IKK_off + flux_IKK_off;

delta_IKK       = delta_IKK     - flux_IKK_off_i;
delta_IKK_i     = delta_IKK_i   + flux_IKK_off_i;

delta_IKK_off   = delta_IKK_off + flux_IKKi;
delta_IKK_i     = delta_IKK_i   - flux_IKKi;



%%
% delta_RIP3 and delta_pMLKL
% label the flux as 92, 93, and degradation flux 94
delta_RIP1=v.k(3)*((C1+C1tnf)/v.K(4)).^v.n(4)./(((C1+C1tnf)/v.K(4)).^v.n(4)+1)-v.kd(3)*RIP1;
% label the flux as 95, inhibition flux 96, and degradation flux 97
delta_RIP3=v.k(1)*(RIP1/v.K(1)).^v.n(1)./((RIP1/v.K(1)).^v.n(1)+1).*1./((a20/v.K(2)).^v.n(2)+1)-v.kd(1)*RIP3;
% label the flux as 98, and degradation flux 99
delta_pMLKL=v.k(2)*(RIP3/v.K(3)).^v.n(3)./((RIP3/v.K(3)).^v.n(3)+1)-v.kd(2)*pMLKL;

%==========================================================
%% Save concentrations for next time step

delta(1,1)=delta_IkBa*v.KO(1);
delta(2,1)=delta_IkBan*v.KO(1);
delta(3,1)=delta_IkBaNFkB*v.KO(1);
delta(4,1)=delta_IkBaNFkBn*v.KO(1);
delta(5,1)=delta_IkBat*v.KO(1);
delta(6,1)=delta_IkBb;
delta(7,1)=delta_IkBbn;
delta(8,1)=delta_IkBbNFkB;
delta(9,1)=delta_IkBbNFkBn;
delta(10,1)=delta_IkBbt;
delta(11,1)=delta_IkBe*v.KO(3);
delta(12,1)=delta_IkBen*v.KO(3);
delta(13,1)=delta_IkBeNFkB*v.KO(3);
delta(14,1)=delta_IkBeNFkBn*v.KO(3);
delta(15,1)=delta_IkBet*v.KO(3);
delta(16,1)=delta_IkBd;%*(v.Perturb(9)-IkBd);
delta(17,1)=delta_IkBdn;
delta(18,1)=delta_IkBdNFkB;
delta(19,1)=delta_IkBdNFkBn;
delta(20,1)=delta_IkBdt;
delta(21,1)=delta_NFkB*v.KO(2);
if v.Perturb(19)==1
%     if t/60<v.Perturb(20)  %synthetic NFkB
%         delta(22,1)=v.Perturb(21);%delta_NFkBn*v.KO(2);
%     else
        delta(22,1)=0;
%     end
else
    delta(22,1)=delta_NFkBn*v.KO(2);
end
delta(23,1)=delta_IKKK_off;
delta(24,1)=delta_IKKK;
delta(25,1)     =delta_IKK_off;
delta(26,1)     =delta_IKK;
delta(27,1)     =delta_IKK_i;
delta(28,1)     = delta_TNF  ;
delta(29,1)     = delta_tnfrm;
delta(30,1)     = delta_TNFR;
delta(31,1)     = delta_TNFRtnf;
delta(32,1)     = delta_C1;
delta(33,1)     = delta_C1_off  ;
delta(34,1)     = delta_C1tnf  ;
delta(35,1)     = delta_C1tnf_off   ;
delta(36,1)     = delta_TTR ;
delta(37,1)     = delta_a20;%*v.Perturb(18) ;% Previoulsy keep the A20 protein sustained...
delta(38,1)     = delta_a20t*v.Perturb(18)  ;
delta(39,1)     = delta_RIP3 ;
delta(40,1)     = delta_pMLKL;
delta(41,1)     = delta_RIP1 ;
end 