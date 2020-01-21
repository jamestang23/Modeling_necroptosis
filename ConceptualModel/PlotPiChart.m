%(c) 2020 Signaling Systems Lab UCLA
%All rights reserved. 
%This MATLAB code plot the pi-chart figure for the bimodal or unimodal cell
%death time distribution. It loads the output file from
%DeathConcept3_PiChartDist.m. Users should check it load existing mat file
%in the corresponding folder. 
%A detailed application of this package is given in the main text. 


%% Load the file
clear;
filefolder=pwd;
filefolder=[filefolder,'\DeathConceptNew'];
Uni=2;
if Uni==1
    filename=[filefolder,'\DeathConceptModel1_0.1_1.mat'];%UniModalityCount
else
    filename=[filefolder,'\DeathConceptModel2_0.1_1.mat'];%BiModalityCount
end
load(filename);


%% Count the unimodal and bimodal death time distribution
Total=length(Output.SurvivalRate);
SurvivalThreshold=40;UniModalThreshold=0.5;
SurvivalCount=0;
ModalityCount=0;
BiModalityCount=0;UniModalityCount=0;
for i=1:Total
    if Output.SurvivalRate(i)<SurvivalThreshold  
        if Output.Unimodality(i)~=0% delete outlier
        ModalityCount=ModalityCount+1;
        UnimodalityPro(ModalityCount)=Output.Unimodality(i);
        if Output.Unimodality(i)<UniModalThreshold
            BiModalityCount=BiModalityCount+1;
        else
            UniModalityCount=UniModalityCount+1;
        end
        end
    else
    SurvivalCount=SurvivalCount+1;
    end
end

DeathCount=Total-SurvivalCount;
X = [SurvivalCount UniModalityCount BiModalityCount];
figure ('position', [00, 10, 700, 500]);
labels = {'Survival ratio>0.4','Unimodal','Bimodal'};
p=pie(X);
legend(labels,'location','southeastoutside')
t = p(6);t.FontSize = 15;
t = p(4);t.FontSize = 15;
t = p(2);t.FontSize = 15;
 set(gca,'linewidth',2);
 set(gca,'FontSize',20);
if Uni==1
    figurename2=[filefolder,'\PieChart.fig'];

else
    figurename2=[filefolder,'\PieChartBimodal.fig'];

end
saveas(gcf,figurename2);  





