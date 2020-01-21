function init = getInit(Type)

init     = struct; % initial concentrations and names
% (can copy this array from spreadsheet)
init.NAMES  = {...
 'IkBa'	...	1
    'IkBan'	...	2
    'IkBaNFkB'	...	3
    'IkBaNFkBn'	...	4
    'IkBat'	...	5
    'IkBb'	...	6
    'IkBbn'	...	7
    'IkBbNFkB'	...	8
    'IkBbNFkBn'	...	9
    'IkBbt'	...	10
    'IkBe'	...	11
    'IkBen'	...	12
    'IkBeNFkB'	...	13
    'IkBeNFkBn'	...	14
    'IkBet'	...	15
    'IkBd'	...	16
    'IkBdn'	...	17
    'IkBdNFkB'	...	18
    'IkBdNFkBn'	...	19
    'IkBdt'	...	20
    'NFkB'	...	21
    'NFkBn'	...	22
    'IKKK_off'	...	23
    'IKKK'	...	24
    'IKK_off'	...	25
    'IKK'	...	26
    'IKK_i'	...	27
    'TNF'	...	28
    'tnfrm'	...	29
    'TNFR'	...	30
    'TNFRtnf'	...	31
    'C1'	...	31
    'C1_off'	...	33
    'C1tnf'	...	34
    'C1tnf_off'	...	35
    'TTR'	...	36
    'a20'	...	37
    'a20t'	...	38
    'RIP3'	...	39
    'pMLKL'	...	40
'RIP1'	...	41
                 };

init.VALS  = zeros(length(init.NAMES),1); % all species start at zero by default

         
% Set initial values for fixed species (can cycle between forms, but are neither created nor destroyed)
init.VALS(strcmp('NFkBn',init.NAMES))   = 0.125;
init.VALS(strcmp('IKKK_off',init.NAMES))   = 0.1;    % IKKK_off
init.VALS(strcmp('IKK_off',init.NAMES))   = 0.1;    % IKK_off
init.VALS(strcmp('TTR',init.NAMES)) = 8.3e-4;


% Notes on steady state
% - [TLR4] is not held constant, but should equilibrate to uM-range based on endocytosis/recycling
% - All IKK should be in IKK_off form before stimulus
