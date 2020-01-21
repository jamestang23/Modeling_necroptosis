function [n,i] = getRateParams()
%% intro
% The model uses 2 parameter lists: one upstream of IKK activation, and one downstream.
% for the 
%   This allows for separate simultion of an individual module
%     if you want to do so.
% return: np and ip

    n = zeros(91,1); % NFkB Activation Module
    i = zeros(39,1); % IKK Activation Module


    %% ----- IkB:NFkB Module (including A20 synthesis/deg) -----
    %           ikba       ikbb        ikbe       ikbd
    n([1:3 72])  = [7e-5;      1e-5;       1e-6;   1e-7]; %1e-7 constitutive txn
    n([4:6 73])  = [8;         0.02;       0.3;    .025]; % inducible txn, 
    n([7:9 74])  = [3;         3;          3 ;      3]; % Hill Coefficient
    n([10:12 75])= [10;         37;         37;     240];%90]; % %Change delay here for IkBs ...%inducible txn delay ## updated,10
    n([13:15 76])= [.035;      3e-3;       4e-3;    2e-3]; % mRNA Degradation
    n([16:18 77])= [0.25;      0.25;       0.25;    .25]; % translation rate

    n([19:21 78])= [0.09;      9e-3;       0.045;   .045]; % Free IkB Import
    n([23:25 79])= [0.012;     0.012;      0.012;   .012]; % Free IkB Export
    n([27:29 80])= [0.276;     0.0276;     0.138;   .276]; % IkB:NFkB Import
    n([30:32 81])= [0.828;     0.414;      0.414;   .414]; % IkB:NFkB Export
    n(22)   = 5.4;   % Free NFkB import
    n(26)   = 0.0048;% Free NFkB export

    n([33:35 82])= [0.12;      0.18;       0.18;    1.4e-3]; % IkB deg cytoplasm
    n([36:38 83])= [0.12;      0.18;       0.18;    1.4e-3]; % IkB deg nucleus
    n([39:41 84])= [6e-5;      6e-5;       6e-5;    6e-5]; % Bound IkB deg cyt
    n([42:44 85])= [6e-5;      6e-5;       6e-5;    6e-5]; % Bound IkB deg nuc

    n([45:47 86])= [30;        30;         30   ;   30]; % IkB:NFkB Asn cyt
    n([48:50 87])= [30;        30;         30   ;   30]; % IkB:NFkB Asn nuc
    n([51:53 88])= [6e-5;      6e-5;       6e-5 ;   6e-5]; % IkB:NFkB Dsn cyt
    n([54:56 89])= [6e-5;      6e-5;       6e-5 ;   6e-5]; % IkB:NFkB Dsn nuc

    n([57:59 90])= [0.36;      0.12;       0.18 ;   0.18]; % Free IkB + IKK
                                                           % deg 102 modified
    n([60:62 91])= [0.36;      0.12;       0.18 ;   0.18 ] ; % IkB:NFkB +
                                                             % IKK deg 106

    % A20 Parameters
    n(63)  = 2e-06;   % constitutive txn
    n(64)  = 0.4;     % inducible txn
    n(65)  = 3;       % Hill Coefficient
    n(66)  = 0;       % inducible txn delay
    n(67)  = 0.035;   % mRNA Degradation
    n(68)  = 0.25;    % translation rate
    n(69)  = 30;      % tsl delay
    n(70)  = 0.0029;  % protein degradation
    n(71)  = 120;     % promoter shutdown (experiments show ~120min)

    %--------------- TNF & IKK module ---------------%
    % tnfrm metabolism (synthesis and degradation)
    i(1)   = 0.0154; % pd_m_tnf 45' half life of exogenous TNF, 49    
    i(2)   = 2e-7;   % tnfrm synthesis (txn, tsl, localization), 36
    i(3)   = 0.0058; % tnfrm --> deg  -- 120' halflife, 37

    % TNF-Independent C1 Activation
    i(4)   = 1e-5;   % 3tnfrm --> TNFR, 38
    i(5)   = 0.1;    % TNFR   --> 3tnfrm, 39
    i(6)   = 0.023;  % TNFR internalization -- 30' halflife, 40 

    i(7)   = 100;    % TNFR + TTR --> C1_off, 41
    i(8)   = 0.1;    % C1_off --> TNFR + TTR, 42
    i(9)   = 30;     % C1_off --> C1, 43
    i(10)  = 2;      % C1     --> C1_off, 44
    i(11)  = 1000;   % C1     --> C1_off (A20 Mediated), 45
    i(12)  = i(8);   % C1     --> TNFR + TTR, 46
    i(13)  = i(6);   % C1_off internalization, 47
    i(14)  = i(6);   % C1 internalization, 48

    % TNF-dependent C1 Activation
    i(15)  = 1100;   % 3tnfrm + tnf --> TNFRtnf, 50
    i(16)  = i(15);  % TNFR + tnf --> TNFRtnf, 51
    i(17)  = 0.021;  % TNFRtnf   --> TNFR + tnf, 52
    i(18)  = i(6);   % TNFRtnf internalization -- 30' halflife, 53
    i(19)  = i(7);   % TNFRtnf + TTR --> C1tnf_off, 54
    i(20)  = i(8);   % C1tnf_off --> TNFRtnf + TTR, 55
    i(21)  = i(9);   % C1tnf_off --> C1tnf, 56
    i(22)  = i(10);  % C1tnf     --> C1tnf_off, 57
    i(23)  = i(11);  % C1tnf     --> C1tnf_off (A20 Mediated), 58 
    i(24)  = i(8);   % C1tnf     --> TNFRtnf + TTR, 59
    i(25)  = i(6);   % C1tnf_off internalization, 60 
    i(26)  = i(6);   % C1tnf internalization, 61
    i(27)  = i(17);  % C1tnf_off --> C1_off + tnf, 62
    i(28)  = i(15);  % C1_off + tnf --> C1tnf_off, 63
    i(29)  = i(17);  % C1tnf    --> C1 + tnf, 64
    i(30)  = i(15);  % C1 + tnf --> C1tnf, 65

    % IKKK
    i(31)  = 500;    % IKKK_off --> IKKK (C1 mediated)?500, 69
    i(32)  = i(31);  % IKKK_off --> IKKK (C1tnf mediated), 70 
    i(33)  = 5e-7;%    5e-7;   % IKKK_off --> IKKK (constitutive),30
    i(34)  = 0.25;    % IKKK     --> IKKK_off (constitutive), 31

    % IKK
    i(35)  = 520;    %  IKK_off  --> IKK (IKKK mediated),32
    i(36)  = 5e-5 ;    % IKK_off  --> IKK (constitutive),33
    i(37)  = 0.02  ;    % IKK      --> IKK_off (constitutive), 34
    i(38)  = 0.15 ;    % IKK      --> IKK_i (constitutive)  , 35
    i(39)  = 0.02;%0.02  ;    % IKKi     --> IKK_off (constitutive) , 36
end
