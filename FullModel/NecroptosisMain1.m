%(c) 2020 Signaling Systems Lab UCLA
%All rights reserved. 
%This MATLAB code implements the mathematical modeling on the signaling dynamics of 
%necroptosis. It plots the time course of all the signaling molecules 
%and the death time distribution.

%If the simulation was run and the data files were alredy saved, the code 
%directly load the data to plot the girues

%A detailed description on the model is given in the main text. 



%% Initialize settings on the simulation 
clear;

sizeScan=300;%The number of cells to simulate
Threshold=0.4;%The threshold of cell death marker protein

plotAllspecies=0;%1 is to plot all species and 0 is not.
Type=1;%WT or 2 is RelKO 
id = struct;

id.Perturb(16)=0; %%TNF pulse experiment: 1 is on and 0 is off.
id.Perturb(17)=12;%% TNF pulse time length. 3, 6, 12, 24

scanA20Constituive=0;% 0 is WT, 1 is scan, 2 is A20-KO
PlotHeatMap=0;%To plot heatmap of constitutive A20 scan: 1 is on and 0 is off.

PlotHeatMapWT=1;%Plot NFkB, etc heatmap: 1 is on and 0 is off.
id.Perturb(19)=0;%Synthetic NFkB as input: 1 is on and 0 is off.
PlotSynTogether=0;%Plot synthetic colorlines together: should be 1 one when id.Perturb(19)=1;


TNFConcentration=0;% Scan TNF concentration: 1 is on and 0 is off. Should turn off when scan A20

%death time distribution by distributing A20
deathTime=1;%plot death time distribution.
NoiseStrength=1;%Noise is not added in the current model

%%
% %%
fig.PlotHeatMapWT=PlotHeatMapWT;
fig.scanA20Constituive=scanA20Constituive;
%% Load the experimental data

[DeathRates, txt]= xlsread('datasummary.xlsx',1,'B5:C779');
% %[data.NuclearNFkB, txt]= xlsread('datasummary.xlsx',2,'C3:D9');
% [A20mRNA, txt]= xlsread('datasummary.xlsx',3,'B3:C10');
% [A20ProteinLong, txt]= xlsread('datasummary.xlsx',3,'B16:C22');
% [A20ProteinShort, txt]= xlsread('datasummary.xlsx',3,'F16:G21');

% %Death rates for 3IKBKO and RelAKO

[DeathRates2, txt]= xlsread('Mutant_data.xlsx',1,'B5:D779'); % WT & RelAKO
[DeathRates3, txt]= xlsread('Mutant_data.xlsx',2,'B5:D779'); % WT & 3IkBKO

%[data.NuclearNFkB, txt]= xlsread('NFkB_mutant_data_all_mod2.xlsx',4,'J12:N18');% WT &3KO& A20-KO &RelAKO
[data.NuclearNFkB, txt]= xlsread('RelARelB_new.xlsx',1,'B3:D9');% WT &3KO only
[NuclearNFkB2, txt]= xlsread('RelARelB_new.xlsx',1,'B10:D16');% WT &3KO only


% %New A20 mRNA and Protein data

[data.A20mRNA_new, txt]= xlsread('total_A20_new2.xlsx',2,'B5:D12');
[A20mRNA3, txt]= xlsread('total_A20_new2.xlsx',2,'B21:D28');% 3IkBKO

[data.A20Protein_new, txt]= xlsread('proteinA20final.xlsx',1,'C4:E11');%xlsread('total_A20_new2.xlsx',1,'B4:D11');
%[data.A20Protein_new, txt]= xlsread('2KO data.xlsx',3,'C5:E12');
[A20Protein3, txt]= xlsread('total_A20_new2.xlsx',1,'B20:D27');% 3IkBKO
% % %  [data.A20Protein_new, txt]= xlsread('NFkBmutantdata_A20.xlsx',2,'E4:F10');
% % % [A20Protein2, txt]= xlsread('NFkBmutantdata_A20.xlsx',2,'E11:F14');% RelAKO
% % % [A20Protein3, txt]= xlsread('NFkBmutantdata_A20.xlsx',2,'E15:F21');% 3IkBKO

% [A20Protein_10_insoluble, txt]= xlsread('NFkBmutantdata_NonSoluableA20.xlsx',3,'F5:G10');
% [A20Protein3_insoluble, txt]= xlsread('NFkBmutantdata_NonSoluableA20.xlsx',3,'H5:I10');% 3IkBKO
% [A20Protein_11_insoluble, txt]= xlsread('NFkBmutantdata_NonSoluableA20.xlsx',3,'F17:G20');
% [A20Protein2_insoluble, txt]= xlsread('NFkBmutantdata_NonSoluableA20.xlsx',3,'H17:I20');% RelAKO
% 
% [A20Protein_new_insoluble_WT, txt]= xlsread('total_A20_new.xlsx',1,'B3:C10');
% [A20Protein_new_insoluble_RelAKO, txt]= xlsread('total_A20_new.xlsx',1,'B11:C18');
% [A20Protein_new_insoluble_3IkBKO, txt]= xlsread('total_A20_new.xlsx',1,'B19:C26');

[data.IkBaProtein, txt]= xlsread('IkBa_data_new.xlsx',1,'B5:D15');%
[IkBamRNA, txt]= xlsread('IkBa_data_new.xlsx',1,'F5:G10');% 


[data.IkBTime, txt]= xlsread('proteinIkBsfinal.xlsx',1,'C7:C17');%
[data.IkBeProtein, txt]= xlsread('proteinIkBsfinal.xlsx',1,'F7:G17');%
[data.IkBdProtein, txt]= xlsread('proteinIkBsfinal.xlsx',1,'H7:I17');%
[data.IkBbProtein, txt]= xlsread('proteinIkBsfinal.xlsx',1,'J7:K17');%



[data.RIP3_insoluble_3IkBKO, txt]= xlsread('NFkB_mutant_data_A20+RIP3.xlsx',3,'H28:J33');% WT & 3IkBKO
%[data.RIP3_insoluble, txt]= xlsread('activeRIP3_pMLKL_mod.xlsx',1,'B18:C25');% WT & 3IkBKO
[data.RIP3_insoluble, txt]= xlsread('RIP3_MLKL final.xlsx',1,'E5:G15');% WT 


%pMLKL
% [pMLKL_soluble, txt]= xlsread('pMLKL_sol_vs_insol.xlsx',1,'A3:C8');% WT & 3IkBKO
% %[data.pMLKL_insoluble, txt]= xlsread('pMLKL_sol_vs_insol.xlsx',1,'A14:C19');% WT & 3IkBKO
% %[data.pMLKL_insoluble, txt]= xlsread('activeRIP3_pMLKL.xlsx',1,'B3:C13');% WT 
[data.pMLKL_insoluble, txt]= xlsread('RIP3_MLKL final.xlsx',1,'B5:D15');% WT 


[C1, txt]= xlsread('C1data.xlsx',1,'B4:C10');
[C1_off, txt]= xlsread('C1data.xlsx',2,'B4:C9');


id.KO=ones(1,3);

%RelA KO
if Type==2
    [data.NuclearNFkB, txt]= xlsread('Mutant_data.xlsx',3,'C3:E9');% WT & 3IkBKO
    [data.A20mRNA_new, txt]= xlsread('total_A20_new2.xlsx',2,'B13:D20');% RelAKO
    [data.A20Protein_new, txt]= xlsread('total_A20_new2.xlsx',1,'B12:D19');% RelAKO
    id.KO(2)=0;
end


%IkBa, IkBe KO
if Type==3
    [data.NuclearNFkB, txt]= xlsread('2KO data.xlsx',1,'C11:E17');
    [data.A20mRNA_new, txt]= xlsread('2KO data.xlsx',2,'C14:E22');
    [data.A20Protein_new, txt]= xlsread('2KO data.xlsx',3,'C13:E20');
    [data.data.IkBdProtein_2KO, txt]= xlsread('proteinp100andIkBbin2KO.xlsx',1,'C3:E11');
    [data.IkBbProtein_2KO, txt]= xlsread('proteinp100andIkBbin2KO.xlsx',1,'F3:G11');
    id.KO(1)=0;id.KO(3)=0;%IkBa, IkBan, IkBat, IkBe, IkBen, IkBet
end


%IkBa KO
if Type==4
    id.KO(1)=0;%IkBa, IkBan, IkBat
end


%% Initiatlize parameters

% BMDM model 

%  Vary i parameters
alter =  [...
 ];
if ~isempty(alter)
    id.inputPid = alter(:,1)';
    id.inputP  = alter(:,2)';
end

% Vary n parameters
alter =       [...
    ];
if ~isempty(alter)
    id.inputvPid = alter(:,1)';
    id.inputvP  = alter(:,2)';
end




%TNF concentration 
if TNFConcentration==1
doses = [10 100] ;%[ .1, 1, 10, 100] ;%[10];%[ .1, 1, 10, 100] ; %ng/ml
else
doses =[10];%[ .1, 1, 10, 100] ; %ng/ml
end 
mod_colormap = jet(1);
%divergingmap(0:1/(length(doses)-1):1,[12 12 77]/255,[158 ...
%                    4 0]/255);

id.Perturb(1)=0.1;%0.01;%0.01 is good one%TNF decay
id.Perturb(2)=1;% C1_off --> C1
id.Perturb(3)=1;% C1_off internalization %TNFR + TTR --> C1_off/// % C1tnf_off --> C1_off + tnf
id.Perturb(4)=1;%A20 decay
%id.Perturb(5)=.08;%A20 mRNA transcription
id.Perturb(6)=20;%A20 translation
id.Perturb(12)=.85;%IkBb transcription basal
id.Perturb(13)=.3;%IkBd transcription basal
id.Perturb(14)=.6;%A20 mRNA transcription basal
%id.Perturb(15)=1;%A20 function to CI 

%For A20 difference in WT and A20 
id.Perturb(5)=.15;%A20 mRNA transcription
%id.Perturb(15)=3;%A20 function to CI 
id.Perturb(15)=1;%A20 function to CI 


id.Perturb(7)=4;%4: %Now, IkBd induced transcription
% 4 is called 2-fold in slide, 2 is called 1-fold, which is bettern for IkBs



if Type==3
%id.Perturb(7)=1;%Previous IkBd protein synthesis rate=0.15 
id.Perturb(8)=1;%Now, threshold of IkBd mRNA syntheis. Indeed, mRNA is already too high in KO, and can set a threshold there.  
%id.Perturb(9)=0.018;%Now, threshold of IkBd protein syntheis
else

end
id.Perturb(9)=20;%Now, IkBb induced transcription 


%The parameters for added necroptosis module:
id.k=[.00015,4,.015];
%id.K=[1*1e-11,1e-1,0.08,1*1e-11];
%id.K=[1*1e-4,1e-1,0.08,1*1e-6];%1*1e-11];
id.K=[1*1e-4,1e-1*1,0.08,1*1e-6];%1*1e-11];
id.n=[3,3,5,3];
id.kd=[.0005,3,.5];
k1=id.k(1);

% id.kd(3)=id.kd(3)/5;
id.K(4)=id.K(4)/10;
id.kd(3)=id.kd(3)/20;
id.k(3)=id.k(3)/20;
id.K(1)=id.K(1)*3;
id.kd(2)=id.kd(2)/10*2.1;
id.k(2)=id.k(2)/10;
id.K(3)=id.K(3)*0.8;

id.Perturb(10)=1;%distribute a20 systhesis :set below
id.Perturb(11)=1;%l

%% Adjust the parameter for the A20 KO
 id.Perturb(18)=1;
A20KO=0;
if A20KO==1 %KO A20 role to left
 id.Perturb(15)=0;%A20 function to CI 
elseif A20KO==2 %KO A20 role to right RIP3
 id.K(2)=1e10;  
elseif A20KO==3 %KO A20 
    id.Perturb(18)=0;%A20KO
end


%% Rename foldername for TNF decay simulation


if id.Perturb(16)==1
filename=['TNFDecay'];
else
filename=[''];
end

%%  Adjust the parameter for synthetic NFkB
if id.Perturb(19)==1
       
    id.Perturb(21)=1/12;%1/12/2;%Conversion factor is 12; 0.0035;
    id.Perturb(22)=3e-3;%basal 3 nM
    id.Perturb(23)=0;%1/3600;%initial time lag in hour
    SyntheticPerb=[0.5,1,2,4,8,16]+id.Perturb(23);
    filename=['Synthetic'];
else
    filename=[''];
    SyntheticPerb=[24];
end

%% Adjust the parameter for the A20 constitutive scan


A20Constituive=id.Perturb(14);
FoldExponent=-10:0.5:10;
Fold=2.^FoldExponent;%[1/16 1/8 1/4 1/2 1 2 4 8 16 2^5 2^10];
if scanA20Constituive==1
    scanA20ConstNum=length(Fold);
    FirstPassageSummary=cell(scanA20ConstNum,1);
   
else
    scanA20ConstNum=1; filename=[filename,'Induced'];
end

%% Adjust the parameter for TNF concentration change
if TNFConcentration==1
    scanTNFConcentration=length(doses);
    %FirstPassageSummary=cell(scanTNFConcentration,1); 
    filename=[filename,'TNFConc']; 
else
    scanTNFConcentration=1; 
end


%% Start the simulation

for nnn=1:size(SyntheticPerb,2)
    id.Perturb(20)=SyntheticPerb(nnn);
    if id.Perturb(19)==1%run synthetic
        if nnn==1
        filename2=filename;
        end
        filename=[filename2,num2str(SyntheticPerb(nnn))];
    end
    
    
for kkk=1:scanTNFConcentration
    if TNFConcentration==1
        if kkk==1
        filename2=filename;
        end
        filename=[filename2,num2str(doses(kkk))];
    end
    
for iii=1:scanA20ConstNum
    
    if scanA20Constituive==1
    id.Perturb(14)=A20Constituive*Fold(iii);
    elseif scanA20Constituive==2
       % id.Perturb(14)=0;
        id.Perturb(18)=0;
    end
%% More parameters settings: species name, time window, time step length 
%id.output = {'TNF','IKKK','IkBd','IKK','C1','NFkBn','a20t','a20'}; % output names are in getInit.m
output2 = {'TNF','C1_off','IkBat','IkBa','NFkBn','a20t','a20','RIP3','pMLKL','C1tnf','RIP1','IkBbt','IkBb','IkBdt','IkBd','IkBet','IkBe','IKK','IKKK'}; % output names are in getInit.m
IndexToShow=[6,8,9,5,4,13,15,17,11,7,1,30];%shoulld be the index in id.output
id.output = {'TNF','C1_off','IkBat','IkBa','NFkBn','a20t','a20','RIP3','pMLKL','C1tnf','RIP1','IkBbt','IkBb','IkBdt','IkBd','IkBet','IkBe',...
    'IkBan','IkBaNFkB','IkBaNFkBn','IkBbn','IkBbNFkB','IkBbNFkBn','IkBdn','IkBdNFkB','IkBdNFkBn','IkBen','IkBeNFkB','IkBeNFkBn','IKK','IKKK',...
    'NFkB','IKKK_off','IKK_off','IKK_i','tnfrm','TNFR','TNFRtnf','C1','C1tnf_off','TTR','IkBa','IkBb','IkBe','IkBd'}; % output names are in getInit.m

id.DT = 5;%0.5; 
id.sim_time = 60*24;
[n,i] = getRateParams(); % Vary "i" params with id.inputPid/inputP and "n" params with id.inputvPid/inputvP
% Simulate
wt_sim = zeros(length(id.output),round(id.sim_time/id.DT)+1, ...
               length(doses(kkk)));
ttt=0:id.DT:id.sim_time;       

%%


for index2=1:1%3% 1 distribute syntheis, 2 distribute degradation, 3 distribute syntheis and degradation

    % Figure names:
fig.figurename2=[pwd,'\DeathConcept2\',filename,'Full','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename3=[pwd,'\DeathConcept2\',filename,'A20Frac','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename4=[pwd,'\DeathConcept2\',filename,'A20Frac2','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename5=[pwd,'\DeathConcept2\',filename,'deathTime','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename6=[pwd,'\DeathConcept2\',filename,'A20mRNA','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename7=[pwd,'\DeathConcept2\',filename,'pMLKL','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename8=[pwd,'\DeathConcept2\',filename,'NFkB','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename9=[pwd,'\DeathConcept2\',filename,'Data','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.mat'];
fig.figurename10=[pwd,'\DeathConcept2\',filename,'deathRates','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename11=[pwd,'\DeathConcept2\',filename,'HeatMap1','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename12=[pwd,'\DeathConcept2\',filename,'HeatMap2','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename121=[pwd,'\DeathConcept2\',filename,'HeatMap31','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename122=[pwd,'\DeathConcept2\',filename,'HeatMap32','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename123=[pwd,'\DeathConcept2\',filename,'HeatMap33','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename124=[pwd,'\DeathConcept2\',filename,'HeatMap34','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];

fig.figurename13=[pwd,'\DeathConcept2\',filename,'Mean_A20mRNA','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename14=[pwd,'\DeathConcept2\',filename,'Mean_pMLKL','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename15=[pwd,'\DeathConcept2\',filename,'Mean_NFkB','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];

fig.figurename21=[pwd,'\DeathConcept2\',filename,'Mean_IkBa','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename22=[pwd,'\DeathConcept2\',filename,'Mean_IkBb','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename23=[pwd,'\DeathConcept2\',filename,'Mean_IkBd','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename24=[pwd,'\DeathConcept2\',filename,'Mean_IkBe','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename25=[pwd,'\DeathConcept2\',filename,'Mean_RIP1','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename26=[pwd,'\DeathConcept2\',filename,'Mean_A20Protein','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename27=[pwd,'\DeathConcept2\',filename,'Mean_TNF','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename28=[pwd,'\DeathConcept2\',filename,'Mean_RIP3','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename29=[pwd,'\DeathConcept2\',filename,'Mean_IKK','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];


%% Plot heatmaps of NFkB and A20 
if PlotHeatMap==1
    
    load(fig.figurename9);
    fig.figurename121=[pwd,'\DeathConcept2\',filename,'HeatMap31','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename122=[pwd,'\DeathConcept2\',filename,'HeatMap32','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename123=[pwd,'\DeathConcept2\',filename,'HeatMap33','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];
fig.figurename124=[pwd,'\DeathConcept2\',filename,'HeatMap34','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.jpg'];

    FirstPassageSummary{iii,1}=FirstPassage;
    NFkBSummary(iii,:)=DataSave{1,1}(4,:); %NFkB;DataSave{:,1}(4,:);
    if iii==scanA20ConstNum
        edges = [0:2:24];
        for kkk=1:scanA20ConstNum
        hh = histcounts(FirstPassageSummary{kkk,1},edges);
        if ~isempty(FirstPassageSummary{kkk,1})
            DeathTimeSummary(kkk,:)=hh*300/sizeScan;
            %DeathTimeSummary(kkk,:)=DeathTimeSummary(kkk,:)/sum(DeathTimeSummary(kkk,:));
        else
            DeathTimeSummary(kkk,:)=zeros(size(hh));
        end
        end
        
        figure('position', [-1600, 10, 660, 600])
        surf(ttt/60,FoldExponent,NFkBSummary);
        h=colorbar; 
        ylabel(h, 'NF\kappaB (A.U.)')
        h.Ticks = linspace(0, 0.08, 5) ; %Create 8 ticks from zero to 1
        h.TickLabels = num2cell([0:0.3:1.2]) ;
        shading interp;
        view(2);
        ylim([min(FoldExponent), max(FoldExponent)]);
        xlim([min(ttt/60), max(ttt/60)]);
         colormap jet;
        caxis([0 0.08]);
        title('NF\kappaB');
        ylabel('constitutive A20 (2^X)')
        xlabel('time (h)')
         set(gca,'FontSize',20);
        saveas(gcf,fig.figurename11);
        
        figure('position', [-1600, 10, 660, 600])
        surf(edges(1:end),FoldExponent,[DeathTimeSummary,zeros(size(DeathTimeSummary,1),1)]);
        h=colorbar; 
        ylabel(h, 'cell count') 
        %ylabel(h, 'Frequency') 
        shading interp;
        %shading flat;
        view(2);
         colormap jet;
        caxis([0 50]);
        %title('death event frequency')
        title('death time distribution')
        ylabel('constitutive A20 (2^X)')
         xlabel('time (h)')
        ylim([min(FoldExponent), max(FoldExponent)]);
        xlim([min(edges(1:end-1)), max(edges(1:end-1))]);
         set(gca,'FontSize',20);
        saveas(gcf,fig.figurename12);
        
     
        figure('position', [-1600, 10, 800, 200])
        surf(edges(1:end),[0 1],[[DeathTimeSummary(1,:),0];[DeathTimeSummary(1,:),0]]);
%         h=colorbar; 
%         ylabel(h, 'cell count') 
        %ylabel(h, 'Frequency') 
        shading interp;
        %shading flat;
        view(2);
        title('A20 KO')
        h=colorbar; 
         colormap jet;
        caxis([0 50]);
        ylabel(h, 'cell count') 
        %ylabel('constitutive A20 (2^X)')
        set(gca,'ytick',[]);
         xlabel('time (h)')
        %ylim([min(FoldExponent), max(FoldExponent)]);
        xlim([min(edges(1:end)), max(edges(1:end))]);
         set(gca,'FontSize',20);
        saveas(gcf,fig.figurename121);
        
        figure('position', [-1600, 10, 800, 200])
        %surf(ttt(1:100:end)/60,[0 1],[NFkBSummary(1,1:100:end);NFkBSummary(1,1:100:end)]);
        surf(ttt/60,[0 1],[NFkBSummary(1,:);NFkBSummary(1,:)]);
        %surf(edges(1:end-1),[0 1],[DeathTimeSummary(end,:);DeathTimeSummary(end,:)]);
%         h=colorbar; 
%         ylabel(h, 'cell count') 
        %ylabel(h, 'Frequency') 
        shading interp;
       % shading flat;
       h=colorbar; 
        colormap jet;
        caxis([0 0.08]);
        ylabel(h, 'NF\kappaB');
        %h.Ticks = linspace(0, max(NFkBSummary2), 3) ; %Create 8 ticks from zero to 1
        h.Ticks = linspace(0, 0.08, 5) ; %Create 8 ticks from zero to 1
        h.TickLabels = num2cell([0:0.3:1.2]) ;
        view(2);
        title('A20 KO')
        %ylabel('constitutive A20 (2^X)')
        set(gca,'ytick',[]);
         xlabel('time (h)')
        %ylim([min(FoldExponent), max(FoldExponent)]);
        xlim([min(edges(1:end)), max(edges(1:end))]);
         set(gca,'FontSize',20);
        saveas(gcf,fig.figurename122);
        
        
        
        
        close all;
%         dddd
    end
    continue;
end

%% Generate simulations on the time course of all species, with distributing the parameters
if deathTime==1 && ~exist(fig.figurename9)
    ToSave=1;
DataSave=cell(sizeScan,1);

IDs=cell(1,size(IndexToShow,2));
tt=0:id.DT:id.sim_time;

%disp(figurename2)
% if exist(figurename2)
%     display('Have it');
%     return;
% end
% dd
    ParaDistribute=zeros(1,sizeScan);
    for kscan=1:sizeScan
       
        solution=[];
        
        if index2==1
%             id.Perturb(10)=lognrnd(0,NoiseStrength);%distribute a20 mRNA systhesis
            %id.Perturb(10)=lognrnd(0,NoiseStrength);%distribute a20 mRNA systhesis
%             if rand(1)<0.3
%                 %id.Perturb(10)=abs(randn(1)*4-3);%lognrnd(0,NoiseStrength)+5;%distribute a20 mRNA systhesis
%                 id.Perturb(10)=abs(randn(1)*1-3);%lognrnd(0,NoiseStrength)+5;%distribute a20 mRNA systhesis
%             elseif rand(1)<0.6
%                 %id.Perturb(10)=abs(randn(1)*14-16);%lognrnd(0,NoiseStrength);%distribute a20 mRNA systhesis
%                 id.Perturb(10)=abs(randn(1)*6-18);%lognrnd(0,NoiseStrength);%distribute a20 mRNA systhesis
%             else
%                 id.Perturb(10)=abs(randn(1)*10-60);
%             end
%             if rand(1)<0.3
%                 %id.Perturb(10)=abs(randn(1)*4-3);%lognrnd(0,NoiseStrength)+5;%distribute a20 mRNA systhesis
%                 id.Perturb(10)=abs(randn(1)*2-3);%lognrnd(0,NoiseStrength)+5;%distribute a20 mRNA systhesis
%             elseif rand(1)<0.55
%                 %id.Perturb(10)=abs(randn(1)*14-16);%lognrnd(0,NoiseStrength);%distribute a20 mRNA systhesis
%                 id.Perturb(10)=abs(randn(1)*8-18);%lognrnd(0,NoiseStrength);%distribute a20 mRNA systhesis
%             else
%                 id.Perturb(10)=abs(randn(1)*10-60);
%             end
%             
%             
%             id.Perturb(10)=id.Perturb(10)/25;
%             ParaDistribute(kscan)=id.Perturb(10);            
%             id.k(1)=k1*abs(randn(1)*0.2+0.7);
            scale2=3;
            randTemp=rand(1);
%             if randTemp<0.45
%                 %id.Perturb(10)=abs(randn(1)*4-3);%lognrnd(0,NoiseStrength)+5;%distribute a20 mRNA systhesis
%                 id.Perturb(10)=abs(randn(1)*10*scale2-3);%lognrnd(0,NoiseStrength)+5;%distribute a20 mRNA systhesis
%             elseif randTemp<0.65
%                 %id.Perturb(10)=abs(randn(1)*14-16);%lognrnd(0,NoiseStrength);%distribute a20 mRNA systhesis
%                 id.Perturb(10)=abs(randn(1)*10*scale2-400);%lognrnd(0,NoiseStrength);%distribute a20 mRNA systhesis
%             else
%                 id.Perturb(10)=abs(randn(1)*30*scale2-500);
%             end
            if randTemp<0.6
                %id.Perturb(10)=abs(randn(1)*4-3);%lognrnd(0,NoiseStrength)+5;%distribute a20 mRNA systhesis
                id.Perturb(10)=abs(randn(1)*10*scale2-3);%lognrnd(0,NoiseStrength)+5;%distribute a20 mRNA systhesis
            elseif randTemp<0.5
                %id.Perturb(10)=abs(randn(1)*14-16);%lognrnd(0,NoiseStrength);%distribute a20 mRNA systhesis
                id.Perturb(10)=abs(randn(1)*10*scale2-400);%lognrnd(0,NoiseStrength);%distribute a20 mRNA systhesis
            else
                id.Perturb(10)=abs(randn(1)*45*scale2-600);
            end
            
            id.Perturb(10)=id.Perturb(10)/25/scale2/1.5;
            ParaDistribute(kscan)=id.Perturb(10);            
            %id.k(1)=k1*abs(randn(1)*0.2+0.7);
%              id.k(1)=k1*abs(randn(1)*1+0.6);
             id.k(1)=k1*lognrnd(0,0.5)*0.7;%k1*lognrnd(0,1)*0.8;%k1*abs(randn(1)*1.8);%k1*lognrnd(0,1)
            
           %id.k(1)=k1*0.7;
            
            
%             %%%%%test original parameters:
%             id.Perturb(10)=1;
%             id.k(1)=k1;
%             %%%%%
            
        elseif  index2==2
            id.Perturb(11)=lognrnd(0,NoiseStrength);%lognrnd(0,.5);%distribute a20 mRNA degradation
            id.Perturb(14)=id.Perturb(11);%hange basal TF rate  to ensure uniform basal value
        else  
            id.Perturb(10)=lognrnd(0,NoiseStrength);%distribute a20  mRNA systhesis
            id.Perturb(11)=lognrnd(0,NoiseStrength);%lognrnd(0,.5);%distribute a20  mRNA degradation
            id.Perturb(14)=id.Perturb(11);%hange basal TF rate  to ensure uniform basal value
        end

if scanA20Constituive==1
    id.Perturb(10)=0;
end
        
        
for idx_d = 1:length(doses(kkk))
    run_id = id;
    run_id.dose = doses(kkk);%doses(idx_d);
    [wt_sim(:,:,idx_d),BASAL,v] = getSimData(run_id,Type);
end

n_output2 = 17;%length(output2);
for idx_o =IndexToShow
    
    if idx_o==4 %IkBa    
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+1:n_output2+3,:,idx_d),1);
        elseif idx_o==13%IkBb
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+4:n_output2+6,:,idx_d),1);
        elseif idx_o==15%IkBd
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+7:n_output2+9,:,idx_d),1);
        elseif idx_o==17%IkBe
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+10:n_output2+12,:,idx_d),1);
        else
            dataplot=wt_sim(idx_o,:,idx_d);
    end
   solution=[solution',dataplot']';
end
%IDs={'A20mRNA','RIP3','pMLKL','NFkBn','IkBa','IkBb','IkBd','IkBe','RIP1','A20Protein','TNF'};
IDs={'A20mRNA','RIP3','pMLKL','NF\kappaBn','I\kappaB\alpha','I\kappaB\beta','I\kappaB\delta','I\kappaB\epsilon','RIP1','A20Protein','TNF','IKK'};
DataSave{kscan,1}=solution;

if v.Perturb(19)==1
        Duration1=round(size(DataSave{kscan,1},2)*(1/2)/24);  %synthetic NFkB
        Duration=round(size(DataSave{kscan,1},2)*v.Perturb(20)/24);  %synthetic NFkB
        DataSave{kscan,1}(4,1:Duration1)=v.Perturb(22);%delta_NFkBn*v.KO(2);
        DataSave{kscan,1}(4,Duration1:Duration)=v.Perturb(21);%delta_NFkBn*v.KO(2);
        DataSave{kscan,1}(4,Duration:end)=0;
end

    if  plotAllspecies==1
        break;
     end
    end
    
   

else
    load(fig.figurename9);
    ToSave=0;fig.PlotHeatMapWT=PlotHeatMapWT;fig.scanA20Constituive=scanA20Constituive;
    display('Have it');
    
    
    fig.figurename2=[pwd,'\DeathConcept2\',filename,'Full','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename3=[pwd,'\DeathConcept2\',filename,'A20Frac','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename4=[pwd,'\DeathConcept2\',filename,'A20Frac2','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename5=[pwd,'\DeathConcept2\',filename,'deathTime','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename6=[pwd,'\DeathConcept2\',filename,'A20mRNA','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename7=[pwd,'\DeathConcept2\',filename,'pMLKL','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename8=[pwd,'\DeathConcept2\',filename,'NFkB','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename9=[pwd,'\DeathConcept2\',filename,'Data','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.mat'];
fig.figurename10=[pwd,'\DeathConcept2\',filename,'deathRates','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename11=[pwd,'\DeathConcept2\',filename,'HeatMap1','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename12=[pwd,'\DeathConcept2\',filename,'HeatMap2','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename121=[pwd,'\DeathConcept2\',filename,'HeatMap31','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename122=[pwd,'\DeathConcept2\',filename,'HeatMap32','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename123=[pwd,'\DeathConcept2\',filename,'HeatMap33','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename124=[pwd,'\DeathConcept2\',filename,'HeatMap34','_A20Const_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];

fig.figurename13=[pwd,'\DeathConcept2\',filename,'Mean_A20mRNA','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename14=[pwd,'\DeathConcept2\',filename,'Mean_pMLKL','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename15=[pwd,'\DeathConcept2\',filename,'Mean_NFkB','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];

fig.figurename21=[pwd,'\DeathConcept2\',filename,'Mean_IkBa','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename22=[pwd,'\DeathConcept2\',filename,'Mean_IkBb','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename23=[pwd,'\DeathConcept2\',filename,'Mean_IkBd','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename24=[pwd,'\DeathConcept2\',filename,'Mean_IkBe','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename25=[pwd,'\DeathConcept2\',filename,'Mean_RIP1','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename26=[pwd,'\DeathConcept2\',filename,'Mean_A20Protein','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename27=[pwd,'\DeathConcept2\',filename,'Mean_TNF','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename28=[pwd,'\DeathConcept2\',filename,'Mean_RIP3','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];
fig.figurename29=[pwd,'\DeathConcept2\',filename,'Mean_IKK','_A20Const_',num2str(iii),'_A20KO',num2str(A20KO),'_',num2str(Threshold),'_',num2str(index2),'_',num2str(sizeScan),'.fig'];



    if PlotSynTogether==1
        for kscan2=1:sizeScan 
         DataSynthetic{nnn}(kscan2,:)=DataSave{kscan2}(1,:);%A20 mRNA;
        end
        if nnn<size(SyntheticPerb,2)      
            continue;
        else
            cc=jet(size(SyntheticPerb,2));
            for ii=1:size(SyntheticPerb,2) 
                MeanTemp=[];
                for kscan2=1:sizeScan 
                    MeanTemp(kscan2,:)=DataSynthetic{ii}(kscan2,:);
                end
                plot(tt,max(mean(MeanTemp,1),0),'-','Color',cc(ii,:),'linewidth',2);hold on;
            
                x2 = [tt, fliplr(tt)];
                MeanTempSort=sort(MeanTemp,1);
                Percentile=0.25;
                inBetween = [MeanTempSort(round(sizeScan*Percentile),:), fliplr(MeanTempSort(round(sizeScan*(1-Percentile)),:))];
                %inBetween = [MeanTempSort(300*Percentile,:), fliplr(MeanTempSort(300*(1-Percentile),:))];

                h2=fill(x2, inBetween, [0.5 0.5 0.5]);
                alpha(0.25);  
            end
             %set(gca,'YLim',[0 200]/1e4/8/3,'Ytick',[0:40:200]/1e4/8/3,'Yticklabel',0:4:20);
             ylabel('\muM');
                ylim([0 3e-3]);
             xlabel3=xlabel('time (h)','FontSize',20);
                xlim([0 tf/3600])
            set(gca,'XLim',[0 24],'Xtick',0:4:24,'Xticklabel',0:4:24);    
            set(gca,'linewidth',2);set(gca,'fontweight','bold','FontSize',20);
            if id.Perturb(19)==1        
                      ylabel(sprintf('%s ',char(Name(1)),'(\muM)'),'FontSize',20)
            end
            
          dddd
        end
    end
end
if plotAllspecies~=1
 [Threshold,FirstPassage]=DeathTime(DataSave,IndexToShow,sizeScan,IDs,Threshold,tt,index2,fig, ToSave,data,Type,id);
end
end




close all;
end
end
end


%% Generate figure for plotting chosen species
if plotAllspecies~=1

color=hsv(size(id.output,2));

for idx_o = [5,6,7,9]%1:n_output2%n_output
    for idx_d = 1:length(doses(kkk))%length(doses)
        if idx_o==4 %IkBa
            
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+1:n_output2+3,:,idx_d),1);
        elseif idx_o==13%IkBb
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+4:n_output2+6,:,idx_d),1);
        elseif idx_o==15%IkBd
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+7:n_output2+9,:,idx_d),1);
        elseif idx_o==17%IkBe
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+10:n_output2+12,:,idx_d),1);
        else
            dataplot=wt_sim(idx_o,:,idx_d);
        end
        
    end
    figure('position', [-1600, 10, 600, 400])
    
    if idx_o==5 %NFkBn
    plot(0:id.DT:id.sim_time,dataplot,'linewidth',1.5,'color','g'); hold on;
    errorbar(data.NuclearNFkB(:,1)*60,data.NuclearNFkB(:,2)/15,data.NuclearNFkB(:,3)/15,'o','markersize',5,'Color','g','MarkerFaceColor','g'); hold on;
   % set(gca,'XLim',[0 id.sim_time],'Xtick',0:120*2:720*2,'Xticklabel',0:4:24)
    set(gca,'YLim',[0 2]/15,'Ytick',[0:0.5:2]/15,'Yticklabel',0:0.5:2)  
    xlabel('Time (hr)');ylabel1=ylabel('A.U.','FontSize',20);
    title(id.output{idx_o});
    end
    
    
    if idx_o==6 %A20mRNA
    plot(0:id.DT:id.sim_time,dataplot,'linewidth',1.5,'color','r'); hold on;
    errorbar(data.A20mRNA_new(:,1)*60,data.A20mRNA_new(:,2)/1e4/8,data.A20mRNA_new(:,3)/1e4/8,'o','markersize',5,'Color','r','MarkerFaceColor','r'); hold on;
   % set(gca,'XLim',[0 id.sim_time],'Xtick',0:120*2:720*2,'Xticklabel',0:4:24)
    set(gca,'YLim',[0 30]/1e4/6,'Ytick',[0:8:40]/1e4/6,'Yticklabel',0:2:10)  
    xlabel('Time (hr)');ylabel1=ylabel('A.U.','FontSize',20);
    title(id.output{idx_o});
    end
    
    if idx_o==7 %data.A20Protein_new
    plot(0:id.DT:id.sim_time,dataplot,'linewidth',1.5,'color','m'); hold on;
    %plot(data.A20Protein_new(:,1)*60,data.A20Protein_new(:,2)/10,'o','markersize',5,'Color','m','MarkerFaceColor','m');hold on;                        
    errorbar(data.A20Protein_new(:,1)*60,data.A20Protein_new(:,2)/10,data.A20Protein_new(:,3)/10,'o','markersize',5,'Color','m','MarkerFaceColor','m');hold on;                        
    %set(gca,'XLim',[0 id.sim_time],'Xtick',0:120*2:720*2,'Xticklabel',0:4:24)
    set(gca,'YLim',[0 0.2],'Ytick',[0:0.1:0.2],'Yticklabel',[0:0.1:0.2])  
    xlabel('Time (hr)');ylabel1=ylabel('A.U.','FontSize',20);
    title(id.output{idx_o});
    end
    
    
    if idx_o==9 %pMLKL
    plot(0:id.DT:id.sim_time,dataplot,'linewidth',1.5,'color','k'); hold on;
    ylim([0 1.5]);
    plot(data.pMLKL_insoluble(:,1)*60,data.pMLKL_insoluble(:,2)/1,'o','markersize',5,'Color','k','MarkerFaceColor','k');hold on;             
    
    set(gca,'YLim',[0 1.5],'Ytick',[0:0.5:1.5],'Yticklabel',[0:0.5:1.5])  
    xlabel('Time (hr)');ylabel1=ylabel('A.U.','FontSize',20);
    title(id.output{idx_o});
    end
    
set(gca,'XLim',[0 id.sim_time*3/4],'Xtick',0:120*2:720*2,'Xticklabel',0:4:24)
%set(gca,'XLim',[0 id.sim_time*3/4],'Xtick',0:240*2:720*2,'Xticklabel',0:4:24)    
set(gca,'linewidth',2);
set(gca,'fontweight','bold','FontSize',30);
    
% if Type==1
% print(['./WT_',num2str(idx_o),'.png'],'-dpng')
% elseif Type==2
% print(['./RelA_',num2str(idx_o),'.png'],'-dpng')
% elseif Type==3
% print(['./IkBaIkBeKO_',num2str(idx_o),'.png'],'-dpng')
% elseif Type==4
% print('./New_IkBaKO.png','-dpng')
% end
end


close all;


else
%% Plot time course of chosen species
for idx_d = 1:length(doses)
    run_id = id;
    run_id.dose = doses(idx_d);
    [wt_sim(:,:,idx_d),BASAL,v] = getSimData(run_id,Type);
end

% plot
n_output = length(id.output);
n_output2 = 17;%length(output2);
if scanA20Constituive==0
figure('Position',[ -1900         27        1815         933]);
end
dataplot=[];
for idx_o = 1:45%n_output2%n_output
    subplot(6,8,idx_o)
    hold on
    for idx_d = 1:length(doses)
        if idx_o==4 %IkBa
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+1:n_output2+3,:,idx_d),1);
        elseif idx_o==13%IkBb
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+4:n_output2+6,:,idx_d),1);
        elseif idx_o==15%IkBd
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+7:n_output2+9,:,idx_d),1);
        elseif idx_o==17%IkBe
            dataplot=wt_sim(idx_o,:,idx_d)+sum(wt_sim(n_output2+10:n_output2+12,:,idx_d),1);
        else
            dataplot=wt_sim(idx_o,:,idx_d);
        end
            
        
        if scanA20Constituive==0
            plot(0:id.DT:id.sim_time,dataplot,'linewidth',1.5,'color',mod_colormap(idx_d,:)); hold on;
            ytickformat('%.2f');
            limsy=get(gca,'YLim');
            set(gca,'Ylim',[0 limsy(2)]);
            dataplotTemp(idx_o,:)=dataplot;
            
            Plotdata.Basal(idx_o,1)=BASAL(idx_o);
            Plotdata.Induced(idx_o,1)=max(dataplot);
        
        elseif scanA20Constituive==2
            plot(0:id.DT:id.sim_time,dataplotTemp(idx_o,:),'linewidth',1.5,'color',mod_colormap(idx_d,:)); hold on;
            plot(0:id.DT:id.sim_time,dataplot,'linewidth',1.5,'color','r'); hold on;
            
            Plotdata.Basal(idx_o,2)=BASAL(idx_o);
            Plotdata.Induced(idx_o,2)=max(dataplot);
            
        end
    end
    set(gca,'XLim',[0 id.sim_time],'Xtick',0:120*2:720*2,'Xticklabel',0:2*2:12*2)  
    xlabel('Time (hr)')
    ylabel(' (\muM)')
    title(id.output{idx_o})
end
%subplot(2,4,8) 
%colormap(mod_colormap)
%hcb=colorbar;
%set(hcb,'YTick',log10(doses),'YTicklabel',{'.001','.01','.1','1','10','100','1000'})
%ylabel(hcb,'TNF doses(ng/ml)') 



%% Plot time course of chosen species for 2IkB KO
color=hsv(size(id.output,2));

%2KO
if Type==3
for kk=1:n_output2
            subplot(6,8,kk)
         if kk==5
                errorbar(data.NuclearNFkB(:,1)*60,data.NuclearNFkB(:,2)/15,data.NuclearNFkB(:,3)/15,'o','markersize',5,'Color','b','MarkerFaceColor','b'); hold on;
                ylim([0 0.2]);              
         elseif kk==6
             ylim([0 6e-4]);
             errorbar(data.A20mRNA_new(:,1)*60,data.A20mRNA_new(:,2)/1e4/5,data.A20mRNA_new(:,3)/1e4/5,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
         elseif kk==10
          elseif kk==2
         elseif kk==15
             errorbar(data.data.IkBdProtein_2KO(:,1)*60,data.data.IkBdProtein_2KO(:,2)/50,data.data.IkBdProtein_2KO(:,3)/50,'o','markersize',5,'Color','b','MarkerFaceColor','b'); hold on;              
          elseif kk==13
             errorbar(data.data.IkBdProtein_2KO(:,1)*60,data.IkBbProtein_2KO(:,1)/40,data.IkBbProtein_2KO(:,2)/40,'o','markersize',5,'Color','b','MarkerFaceColor','b'); hold on;              
         elseif kk==3          
            % plot(IkBamRNA(:,1)*60,IkBamRNA(:,2)/max(IkBamRNA(:,2))/35,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:));hold on;             
            ylim([0 0.04]);
         elseif kk==4
             ylim([0 0.04]);
             %errorbar(data.IkBaProtein(:,1)*60,data.IkBaProtein(:,2)/max(data.IkBaProtein(:,2))/35,data.IkBaProtein(:,3)/max(data.IkBaProtein(:,2))/35,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
         elseif kk==15
             ylim([0 4e-3]);
           % errorbar(data.IkBTime(:,1)*60,data.IkBdProtein(:,1)/400,data.IkBdProtein(:,2)/400,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
         elseif kk==17
             ylim([0 4e-3]);
            %errorbar(data.IkBTime(:,1)*60,data.IkBeProtein(:,1)/400,data.IkBeProtein(:,2)/400,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
          elseif kk==13
              ylim([0 10e-3]);
            %errorbar(data.IkBTime(:,1)*60,data.IkBbProtein(:,1)/200,data.IkBbProtein(:,2)/200,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;

         elseif kk==7
             ylim([0 0.2]);
             plot(data.A20Protein_new(:,1)*60,data.A20Protein_new(:,2)/10,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:));hold on;                        
        elseif kk==8
            ylim([0 0.15]);     
             %plot(data.RIP3_insoluble(:,1)*60,data.RIP3_insoluble(:,2)/max(data.RIP3_insoluble(:,2))/1e1,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:));hold on;             
         elseif kk==9
             ylim([0 1.5]);
            % plot(data.pMLKL_insoluble(:,1)*60,data.pMLKL_insoluble(:,2)/1,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:));hold on;             
         end
end
end

if Type==6 || Type==4
for kk=1:n_output2
            subplot(6,8,kk)
%          if kk==5
%                 errorbar(data.NuclearNFkB(:,1)*60,data.NuclearNFkB(:,2)/15,data.NuclearNFkB(:,3)/15,'o','markersize',5,'Color','b','MarkerFaceColor','b'); hold on;
%                 ylim([0 0.2]);              
%          elseif kk==6
%              ylim([0 6e-4]);
%              errorbar(data.A20mRNA_new(:,1)*60,data.A20mRNA_new(:,2)/1e4/5,data.A20mRNA_new(:,3)/1e4/5,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
%          elseif kk==10
%           elseif kk==2
%          elseif kk==1
%          elseif kk==3          
%              plot(IkBamRNA(:,1)*60,IkBamRNA(:,2)/max(IkBamRNA(:,2))/35,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:));hold on;             
%             ylim([0 0.04]);
%          elseif kk==4
%             % ylim([0 0.04]);
%              errorbar(data.IkBaProtein(:,1)*60,data.IkBaProtein(:,2)/7,data.IkBaProtein(:,3)/7,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
%          elseif kk==15
%              %ylim([0 4e-3]);
%             errorbar(data.IkBTime(:,1)*60,data.IkBdProtein(:,1)/50,data.IkBdProtein(:,2)/50,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
%          elseif kk==17
%              %ylim([0 4e-3]);
%             errorbar(data.IkBTime(:,1)*60,data.IkBeProtein(:,1)/70,data.IkBeProtein(:,2)/70,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
%           elseif kk==13
%               % ylim([0 10e-3]);
%             errorbar(data.IkBTime(:,1)*60,data.IkBbProtein(:,1)/40,data.IkBbProtein(:,2)/40,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
% 
%          elseif kk==7
%              ylim([0 0.2]);
%              plot(data.A20Protein_new(:,1)*60,data.A20Protein_new(:,2)/10,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:));hold on;                        
%         elseif kk==8
%             ylim([0 0.15]);     
%              plot(data.RIP3_insoluble(:,1)*60,data.RIP3_insoluble(:,2)/max(data.RIP3_insoluble(:,2))/1e1,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:));hold on;             
%          elseif kk==9
%              ylim([0 1.5]);
%              plot(data.pMLKL_insoluble(:,1)*60,data.pMLKL_insoluble(:,2)/1,'o','markersize',5,'Color',color(kk,:),'MarkerFaceColor',color(kk,:));hold on;             
%  end     


if kk==5
    if Type~=2
     errorbar(data.NuclearNFkB(:,1)*60,data.NuclearNFkB(:,2)/15*1.2,data.NuclearNFkB(:,3)/15*1.2,'o','markersize',5,'Color','g','MarkerFaceColor','g'); hold on;
    end
end
if kk==6
    errorbar(data.A20mRNA_new(:,1)*60,data.A20mRNA_new(:,2)/1e4/8,data.A20mRNA_new(:,3)/1e4/8,'o','markersize',5,'Color',[1 0.4 0.6],'MarkerFaceColor',[1 0.4 0.6]); hold on; 
end
if Type==1 % data for WT only
if kk==4
    errorbar(data.IkBaProtein(:,1)*60,data.IkBaProtein(:,2)/6,data.IkBaProtein(:,3)/6,'o','markersize',5,'Color','g','MarkerFaceColor','g'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
elseif kk==9
   %plot(data.pMLKL_insoluble(:,1),data.pMLKL_insoluble(:,2)/2.5,'o','markersize',5,'Color','c','MarkerFaceColor','c');hold on;             
    errorbar(data.pMLKL_insoluble(:,1)*60,data.pMLKL_insoluble(:,2)/2.5,data.pMLKL_insoluble(:,3)/2.5,'o','markersize',5,'Color','c','MarkerFaceColor','c');hold on;             
elseif kk==15
    errorbar(data.IkBTime(:,1)*60,data.IkBdProtein(:,1)/100*4*2,data.IkBdProtein(:,2)/100*4*2,'o','markersize',5,'Color','g','MarkerFaceColor','g'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
 elseif kk==17
    errorbar(data.IkBTime(:,1)*60,data.IkBeProtein(:,1)/70,data.IkBeProtein(:,2)/70,'o','markersize',5,'Color','g','MarkerFaceColor','g'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
  elseif kk==13
    errorbar(data.IkBTime(:,1)*60,data.IkBbProtein(:,1)/40*0.9,data.IkBbProtein(:,2)/40*0.9,'o','markersize',5,'Color','g','MarkerFaceColor','g'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
  elseif kk==7
      errorbar(data.A20Protein_new(:,1)*60,data.A20Protein_new(:,2)/10,data.A20Protein_new(:,3)/10,'o','markersize',5,'Color',[1 0.4 0.6],'MarkerFaceColor',[1 0.4 0.6]);hold on;                        
    %plot(data.A20Protein_new(:,1)*60,data.A20Protein_new(:,2)/10,'o','markersize',5,'Color',[1 0.4 0.6],'MarkerFaceColor',[1 0.4 0.6]);hold on;                        
  elseif kk==9 %RIP1

  elseif kk==8 %RIP3
      errorbar(data.RIP3_insoluble(:,1)*60,data.RIP3_insoluble(:,2)/max(data.RIP3_insoluble(:,2))/1e1,data.RIP3_insoluble(:,3)/max(data.RIP3_insoluble(:,2))/1e1,'o','markersize',5,'Color','c','MarkerFaceColor','c');hold on;                       
end
end


end
end
% xlabel3=xlabel('Time  (h)','FontSize',26, 'FontWeight', 'bold');
%        
% xlim([0 20]);xticks([0 4 8 12 16 20]);xticklabels({'0','4','8','12','16','20'});
% ylabel1=ylabel('A. U.','FontSize',26, 'FontWeight', 'bold');
% 
% set(gca,'FontSize',26, 'FontWeight', 'bold');
% tt=title(Name(kk));
% set(gca,'linewidth',2);
% set(gca,'FontSize',26, 'FontWeight', 'bold');
% figurename2=[filefolder,'\WT_',num2str(kk),'.jpg'];
%  saveas(gcf,figurename2); 


%%
% save 
if scanA20Constituive==0
if Type==1
print('./New_WT.png','-dpng')
elseif Type==2
print('./New_RelAKO.png','-dpng')
elseif Type==3
print('./New_IkBaIkBeKO.png','-dpng')
elseif Type==4
print('./New_IkBaKO.png','-dpng')
end
end
end

close all;

%hold on;
%end

% for ii=1:2
% figure('Position',[ -1900         27        1815         933]);
% for idx_o = 1:45
%     subplot(6,8,idx_o)
%     hold on;     
%      bar([Plotdata.Basal(idx_o,ii),Plotdata.Induced(idx_o,ii)])
%      %set(gca,'XLim',[0 id.sim_time],'Xtick',0:120*2:720*2,'Xticklabel',0:2*2:12*2)  
%     %xlabel('Time (hr)')
%     ylabel(' (\muM)')
%     title(id.output{idx_o})
%     
% end
% if ii==1
%     print('./figs/New_WT_Basal.png','-dpng')
% elseif ii==2
% print('./figs/New_A20KO_Basal.png','-dpng')
% end
% end
%close all;