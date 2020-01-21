%(c) 2020 Signaling Systems Lab UCLA
%All rights reserved. 
%This MATLAB code implements the mathematical modeling on the samll 
%network motifs for different types of cell death. It plots the time course
%of the signaling molecules and the death time distribution
%A detailed application of this package is given in the main text. 


%% Initialize the iteration
clear;
filefolder=pwd;
filefolder=[filefolder,'\DeathConceptNew'];
mkdir(filefolder)
Name={'Ligand','Adaptor','Effector','Inhibitor','Regulator','Death Time'};
%Name={'TNF','RIP1/3','pMLKL','NFkB','IkB','Death Time'};
Name={'RIP1/3','pMLKL','X','IkB','NFkB','Death Time'};
global Para 

fignum=1;
tf=24*3600; % hour
SpecieNum=5;
TimeResolution=200;
tt=linspace(0,tf,TimeResolution);
scale=3600;
index2=1;

% Para.k=[.1,.1,5,5,1]/scale;
% Para.Km=[1,0.05,0.5,.2,0.3,.3,0.03];
% Para.kd=[.01,.1,1,1,1]/scale;
Para.k=[.1,.1,5,5,1]/scale;
Para.Km=[1,0.15,0.5,.2,0.3,.3,0.03];
Para.kd=[.01,.1,1,1,1]/scale;
Para.n=[3,3,3,3,3,3,3];
Para.kTemp=Para.k;
Para.Km(5)=0.07;
Threshold=0.1;%should be <=0.3, >0.1

%Initial value and input
Para.InitialValue=zeros(1,SpecieNum);
Para.LigandStrength=1;
Para.HalfLife=10*scale;
TNFPlot=Para.LigandStrength*exp(-tt/Para.HalfLife);

FreePara=0;%0.05;%0 or 0.05
Para.InitialValue(3)=FreePara; %0.05 is proper effective inhibition strength
Para.InitialValue(4)=0.05;
sizeScan=300;


scale2=1;


DataSave=cell(sizeScan,1);

%% Start simulation on each cell
 for kscan=1:sizeScan
     
%Distribute the synthesis rates
% Para.k(1)=Para.kTemp(1)*lognrnd(0,.5);%have adjusted the Var to visually match distribution width; %abs(rand(1)-1);
% Para.k(2)=Para.kTemp(2)*lognrnd(0,.5);

options = odeset('RelTol',1e-5);%,'Stats','on');%,'OutputFcn',@odeplot);
if FreePara~=0
    index=1;
    Para.InitialValue(3)=FreePara*lognrnd(0,.3);
    
[T,solution] = ode45(@ODESimulation_1, tt, Para.InitialValue);%,options);
else    
    
    %Para.InitialValue(3)=0.1*lognrnd(0,.5);
    Para.InitialValue(5)=0.01;
    %Para.k(3)=Para.kTemp(3)*10;
    %Para.k(3)=Para.kTemp(3)*lognrnd(0,0.5);
    if rand(1)<0.5
        Para.k(3)=.5/scale*abs(randn(1)*8+7)*scale2;%1/scale*abs(randn(1)*1.5-3);%lognrnd(0,NoiseStrength)+5;%distribute a20 mRNA systhesis
    else
        Para.k(3)=.5/scale*abs(randn(1)*2+14)*scale2;%1/scale*abs(randn(1)*8-14);%lognrnd(0,NoiseStrength);%distribute a20 mRNA systhesis
    end

    DistributePara(kscan)=Para.k(3);
    
[T,solution] = ode45(@ODESimulation_2, tt, Para.InitialValue);%,options);
    index=2;
end
solution=solution';
DataSave{kscan,1}=solution;
 end
 if index==2
 figure 
 histogram(DistributePara)
 end
 %dd
%% Thrshold for effector

FirstPassage=[];ThresholdIndex=[];
for kscan=1:sizeScan
FreePara=find((DataSave{kscan,1}(2,:)-Threshold)>0);
if size(FreePara,2)>0
FirstPassage=[FirstPassage,tt(FreePara(1))/3600];
ThresholdIndex=[ThresholdIndex,FreePara(1)];
else
 ThresholdIndex=[ThresholdIndex,size(DataSave{kscan,1},2)];   
end
end



%% Plot figures of the simulated time coures and death time distribution

c=jet(sizeScan);
SpeciesToPlot=[2,3,5,6];

for kk=1:4
%figure(1);clf;
figure('position', [-1600, 10, 600, 400])
    


if kk<=3
for kscan=1:sizeScan
    if kk==1
        plot(tt(1:ThresholdIndex(kscan))./3600,max(DataSave{kscan,1}(SpeciesToPlot(kk),1:ThresholdIndex(kscan)),0),'-','Color',c(kscan,:),'linewidth',2);hold on;
        h2=plot(tt(ThresholdIndex(kscan):end)./3600,max(DataSave{kscan,1}(SpeciesToPlot(kk),ThresholdIndex(kscan):end),0),'-','Color',c(kscan,:),'linewidth',2);hold on;
        h2.Color(4)=0.1;  
    else
plot(tt./3600,max(DataSave{kscan,1}(SpeciesToPlot(kk),:),0),'-','Color',c(kscan,:),'linewidth',2);hold on;
    end
end
    if kk==1
        plot(tt./3600,Threshold*ones(1,TimeResolution),'-','Color','k','linewidth',2);hold on;
        
        %plot(tt./3600,max(DataSave{kscan,1}(SpeciesToPlot(kk),:),0),'-','Color',c(kscan,:),'linewidth',2);hold on;
        
        
        ylim([0 Threshold*2]);yticks([0:0.05:0.2]);  yticklabels({'0','0.05','0.1','0.15','0.2'});%yticklabels({'0','0.25','0.5','0.75','1'});
    end
end
xlabel3=xlabel('Time  (h)','FontSize',20);
xlim([0 tf/3600])
ylabel1=ylabel('A.U.','FontSize',20);
if kk==3
    ylim([-.01 .2]);
    yticks([0:0.05:0.2]);
    yticklabels({'0','0.25','0.5','0.75','1'});
end
if kk==2
    ylim([0 0.5]);
    yticks([0 0.25 0.5]);
   % yticklabels({'0','0.25','0.5','0.75','1'});
    
end
if kk==4
    if index==1
        histogram(FirstPassage,24,'BinLimits',[0,24],'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5],'Normalization','probability'); 
    else
        histogram(FirstPassage,24,'BinLimits',[0,24],'FaceColor','b','EdgeColor','b','Normalization','probability');
    end
   xlabel3=xlabel('Time  (h)','FontSize',20);
    xlim([0 tf/3600])
    %ylabel1=ylabel('Cell count','FontSize',20);ylim([0 150]);
%     ylabel1=ylabel('Cell fraction','FontSize',20);
%     %ylim([0 .5]);
%     ylim([0 150]);
    
   % h1 = histogram(FirstPassage,edges,'FaceColor','b','EdgeColor','b','Normalization','probability');  
    ylabel1=ylabel('death event frequency','FontSize',20);
    ylim([0 0.5]);   
    
    
end

%if kk<3
set(gca,'XLim',[0 24],'Xtick',0:4:24,'Xticklabel',0:4:24);
%set(gca,'XLim',[0 id.sim_time*3/4],'Xtick',0:240*2:720*2,'Xticklabel',0:4:24)    
set(gca,'linewidth',2);
set(gca,'fontweight','bold','FontSize',20);
%set(gca,'FontSize',20);
if kk==4
    h=legend(['Survival ',num2str(round(100-length(FirstPassage)/sizeScan*100)),'%'],'Location','northeast');
else
h=legend(Name(SpeciesToPlot(kk)),'Location','northeast');
end
set(h,'fontweight','bold','FontSize',20);
%set(gca,'linewidth',2);
%set(gca,'FontSize',20);
figurename2=[filefolder,'\DeathConceptModel',num2str(index),'_',num2str(Threshold),'_',num2str(kk),'.fig'];
 saveas(gcf,figurename2);   
end
   
 

 

c=jet(sizeScan);
%figure(1);clf;
figure('position', [-1600, 10, 1600, 800])
    
for kk=1:6
    %figure(kk),clf,
subplot(2,3,kk), hold on
%if kk==1
   % plot(tt./3600,max(TNFPlot,0),'-','Color',c(kk,:),'linewidth',2);hold on;
%elseif
if kk<6
for kscan=1:sizeScan
plot(tt./3600,max(DataSave{kscan,1}(kk,:),0),'-','Color',c(kscan,:),'linewidth',2);hold on;
end
    if kk==2
        plot(tt./3600,Threshold*ones(1,TimeResolution),'-','Color','k','linewidth',2);hold on;
        ylim([0 Threshold*2])
    end
    if kk==3 %&& Para.InitialValue(3)~=0
        ylim([0 1]);%Para.InitialValue(3)*2])
    end
end
xlabel3=xlabel('Time  (h)','FontSize',20);
xlim([0 tf/3600])
ylabel1=ylabel('Activity','FontSize',20);

if kk==6
   hist(FirstPassage);
   xlabel3=xlabel('Time  (h)','FontSize',20);
    xlim([0 tf/3600])
    ylabel1=ylabel('Cell count','FontSize',20);
    ylim([0 150])
end
set(gca,'FontSize',20);
h=legend(Name(kk),'Location','northeast');
set(h,'FontSize',20);
set(gca,'linewidth',2);
set(gca,'FontSize',20);

end
figurename2=[filefolder,'\DeathConceptModel',num2str(index),'_',num2str(Threshold),'_',num2str(index2),'.fig'];
while exist(figurename2)
    index2=index2+1;
    figurename2=[filefolder,'\DeathConceptModel',num2str(index),'_',num2str(Threshold),'_',num2str(index2),'.fig'];
end
 saveas(gcf,figurename2);      
 
filename=[filefolder,'\DeathConceptModel',num2str(index),'_',num2str(Threshold),'_',num2str(index2),'.mat'];
%save(filename,'Para','DataSave','tt','tf','Name','TNFPlot','Threshold','FirstPassage');
save(filename,'Para','tt','tf','Name','TNFPlot','Threshold','FirstPassage');
 
close all;
%% Function of ODE
function xp = ODESimulation_2(tt, x)
global Para 
%Derivative of Logistic function
%InputA=LA*(K*exp(-K*(t-Duration)))./(exp(-K*(t-Duration))+1).^2+epsilon; %Derivative of Logistic
%InputI=LI*(K*exp(-K*(t-Duration)))./(exp(-K*(t-Duration))+1).^2+epsilon; %Derivative of Logistic

%Exponential:
%TRIALExpRate=-log(1/2)/(10*3600);
%InputA=kA*exp(-TRIALExpRate*(t-Duration))*heaviside(t-Duration);
%InputI=kI*exp(-TRIALExpRate*(t-Duration))*heaviside(t-Duration);

k=Para.k;
Km=Para.Km;
n=Para.n;
kd=Para.kd;
TNF=Para.LigandStrength*exp(-tt/Para.HalfLife);

xp(1) = k(1).*(TNF/Km(1)).^n(1)./((TNF/Km(1)).^n(1)+1).*1./((x(3)/Km(5)).^n(5)+1)-kd(1)*x(1);% RIP1

xp(2) =  k(2).*(x(1)/Km(2)).^n(2)./((x(1)/Km(2)).^n(2)+1)-kd(2).*x(2);%pMLKL

xp(3) =  k(3).*(x(5)/Km(3)).^n(3)./((x(5)/Km(3)).^n(3)+1)-kd(3).*x(3);%Now, X

xp(4) = k(4).*(x(5)/Km(4)).^n(4)./((x(5)/Km(4)).^n(4)+1).*1./((TNF/Km(6)).^n(6)+1)-kd(4)*x(4);%IkB


xp(5) =  k(5)*1./((x(4)/Km(7)).^n(7)+1)-kd(5).*x(5);%NFkB.


xp=xp';

end


function xp = ODESimulation_1(tt, x)
global Para 
%Derivative of Logistic function
%InputA=LA*(K*exp(-K*(t-Duration)))./(exp(-K*(t-Duration))+1).^2+epsilon; %Derivative of Logistic
%InputI=LI*(K*exp(-K*(t-Duration)))./(exp(-K*(t-Duration))+1).^2+epsilon; %Derivative of Logistic

%Exponential:
%TRIALExpRate=-log(1/2)/(10*3600);
%InputA=kA*exp(-TRIALExpRate*(t-Duration))*heaviside(t-Duration);
%InputI=kI*exp(-TRIALExpRate*(t-Duration))*heaviside(t-Duration);

k=Para.k;
Km=Para.Km;
n=Para.n;
kd=Para.kd;
TNF=Para.LigandStrength*exp(-tt/Para.HalfLife);

xp(1) = k(1).*(TNF/Km(1)).^n(1)./((TNF/Km(1)).^n(1)+1).*1./((x(3)/Km(5)).^n(5)+1)-kd(1)*x(1);% RIP1

xp(2) =  k(2).*(x(1)/Km(2)).^n(2)./((x(1)/Km(2)).^n(2)+1)-kd(2).*x(2);% pMLKL

xp(3) = 0;% -kd(3).*x(3)/3;%0 %+sqrt(2*D(2))*randn(1); % % NFkB. Now, X

xp(4) = 0;%+sqrt(2*D(1))*randn(1); % IkB
xp(5) =0;
xp=xp';

end


    