%This script generates the figures of 
%time course of the modeling species and 
%death time distribution

function [Threshold,FirstPassage]=DeathTime(DataSave,IndexToShow,sizeScan,IDs,Threshold,tt,index2,fig,ToSave,data,Type,id)
filefolder=pwd;
filefolder=[filefolder,'\DeathConcept2'];
mkdir(filefolder);
Name=IDs;%{'Ligand','Adaptor','Effector','Inhibitor','Regulator','death Time'};
global Para 

%figurename2=[filefolder,'\DeathTime','_',num2str(Threshold),'_',num2str(index2),'.jpg'];
% if exist(figurename2)
%     display('Have it');
%     return;
% end


fignum=1;
if ToSave==1
tt=tt/60;
end
tf=tt(end); % hour
SpecieNum=4;
TimeResolution=1000;
%tt=linspace(0,tf,TimeResolution);
scale=3600;
%index2=1;
index=1;

Para.k=[.1,.1,5,2]/scale;
Para.Km=[1,0.1,0.5,.15,0.05,0.1];
Para.n=[3,3,3,3,3,3];
Para.kd=[.01,.1,1,1]/scale;
Para.kTemp=Para.k;

%Initial value and input
Para.InitialValue=zeros(1,SpecieNum);
Para.LigandStrength=1;
Para.HalfLife=10*scale;
TNFPlot=Para.LigandStrength*exp(-tt/Para.HalfLife);

Para.InitialValue(3)=0; %0.05 is proper effective inhibition strength

%sizeScan=300;
% DataSave=cell(sizeScan,1);
% 
%  for kscan=1:sizeScan
%      
% %Distribute the synthesis rates
% Para.k(1)=Para.kTemp(1)*lognrnd(0,.5);%have adjusted the Var to visually match distribution width; %abs(rand(1)-1);
% Para.k(2)=Para.kTemp(2)*lognrnd(0,.5);
% 
% options = odeset('RelTol',1e-5);%,'Stats','on');%,'OutputFcn',@odeplot);
% if Para.InitialValue(3)~=0
%     index=1;
% [T,solution] = ode45(@ODESimulation_1, tt, Para.InitialValue);%,options);
% else    
% [T,solution] = ode45(@ODESimulation_2, tt, Para.InitialValue);%,options);
%     index=2;
% end
% solution=solution';
% DataSave{kscan,1}=solution;
%  end
%%
%Thrshold for effector
%should be <=0.3, >0.1
FirstPassage=[];ThresholdIndex=[];
for kscan=1:sizeScan
Temp=find((DataSave{kscan,1}(3,:)-Threshold)>0);
if size(Temp,2)>0
    assignin('base', 'Temp', Temp);
FirstPassage=[FirstPassage,tt(Temp(1))];
ThresholdIndex=[ThresholdIndex,Temp(1)];
else
 ThresholdIndex=[ThresholdIndex,size(DataSave{kscan,1},2)];   
end
end

%%
%Thrshold for A20: 100 threshold, 10 time points;
% for mRNA now

%ThresholdA20=0:0.01:1;
%ThresholdA20Chose=[0.05,0.1,0.2,0.5];
ThresholdA20=linspace(0,1e-3,100); % for mRNA now
Chose=[1,10,20,40,60,100];
ThresholdA20Chose=ThresholdA20(Chose);%[0.05,0.1,0.2,0.5]*10*1e-4;% for mRNA now
TimePointsA20=[0, 0.5, 1, 2, 4, 8, 12,16]; %hours
OverA20Fraction=zeros(size(TimePointsA20,2),size(ThresholdA20,2));
for kscan=1:sizeScan
    for jj=1:size(TimePointsA20,2)
        Points=min(round(TimePointsA20(jj)/tf*size(tt,2))+1,size(tt,2));
        for kk=1:size(ThresholdA20,2)
            if DataSave{kscan,1}(1,Points)>ThresholdA20(kk)
                OverA20Fraction(jj,kk)=OverA20Fraction(jj,kk)+1;
            end
        end
        
    end
end
OverA20Fraction=OverA20Fraction/sizeScan;

% figure('position', [100, 10, 800, 600])
% cc=jet(size(TimePointsA20,2));
% for jj=1:size(TimePointsA20,2)
% plot(ThresholdA20,OverA20Fraction(jj,:),'-','linewidth',2,'Color',cc(jj,:));hold on;
% end
% %ylim([0.001 1])
% %set(gca, 'YScale', 'log')
% h=legend(num2str(TimePointsA20'),'Location','best');
% set(h,'FontSize',20);
% set(gca,'linewidth',2);
% set(gca,'FontSize',20);
% saveas(gcf,fig.figurename3);  
% 
% 
% figure('position', [100, 10, 800, 600])
% cc=jet(size(ThresholdA20Chose,2));l=0;
% for kk=1:size(ThresholdA20,2)
%     if min(abs(ThresholdA20(kk)-ThresholdA20Chose))==0
%         l=l+1;
%          plot(TimePointsA20,OverA20Fraction(:,kk),'-o','markersize',10,'Color',cc(l,:),'linewidth',2);hold on;
%         % disp(TimePointsA20);disp(OverA20Fraction(:,kk));dd
%     end
% end
% h=legend(num2str(ThresholdA20Chose'),'Location','best');
% set(h,'FontSize',20);
% set(gca,'linewidth',2);
% set(gca,'FontSize',20);
% saveas(gcf,fig.figurename4);  
%%
%Plot mean from single-cell simulations
c=jet(sizeScan);
SpeciesToPlot=[2,3,5,6];

for kk=[1 2 3 4 5 6 7 8 9 10 11 12]
%figure(1);clf;
figure('position', [-1600, 10, 600, 400])

MeanTemp=[];
for kscan=1:sizeScan 
    MeanTemp(kscan,:)=DataSave{kscan,1}(kk,:);
end
%MeanTemp=MeanTemp/sizeScan;
if kk==4
    plot(tt,max(mean(MeanTemp,1),0),'-','Color','r','linewidth',2);hold on;
else
plot(tt,max(mean(MeanTemp,1),0),'-','Color','r','linewidth',2);hold on;
end
x2 = [tt, fliplr(tt)];
%inBetween = [mean(MeanTemp,1)-std(MeanTemp,1)*0.5, fliplr(mean(MeanTemp,1)+0.5*std(MeanTemp,1))];

MeanTempSort=sort(MeanTemp,1);
Percentile=0.3;
inBetween = [MeanTempSort(round(sizeScan*Percentile),:), fliplr(MeanTempSort(round(sizeScan*(1-Percentile)),:))];


h2=fill(x2, inBetween, [0.5 0.5 0.5]);
alpha(0.25);


%%
%normalize for y-axis
if kk==3%pMLKL_insoluble
    set(gca,'YLim',[0 0.6]/1.08,'Ytick',[0:0.1:0.6]/1.08,'Yticklabel',[0:0.25:1.5]) ;    
        %set(gca,'YLim',[0 .2],'Ytick',[0:0.05:0.2],'Yticklabel',[0:0.05:0.2]) ;
   %     set(gca,'YLim',[0 .5],'Ytick',[0:0.1:.5],'Yticklabel',[0:0.1:.5]) ;       
elseif kk==4%NuclearNFkB
    NFkBSummary1=tt;NFkBSummary2=max(mean(MeanTemp,1),0);
        %ylim([0 0.2]);
        set(gca,'YLim',[0 0.08]*1.2,'Ytick',[0:0.02:0.08]*1.2,'Yticklabel',[0:0.3:1.2]);
elseif kk==1%A20mRNA_new
     %set(gca,'YLim',[0 100]/1e4/8/3,'Ytick',[0:20:100]/1e4/8/3,'Yticklabel',0:2:10);
   % set(gca,'YLim',[0 100]/1e4/8/3,'Ytick',[0:20:100]/1e4/8/3,'Yticklabel',0:2:10);
   if id.Perturb(19)==1
       ylabel('\muM');
       ylim([0 3e-3]);
   else
    set(gca,'YLim',[0 200]/1e4/8/3*1.2,'Ytick',[0:40:200]/1e4/8/3*1.2,'Yticklabel',[0:1:5]);
   end
elseif kk==5%IkBaProtein
    set(gca,'YLim',[0 0.24]/6*5,'Ytick',[0:0.04:0.24]/6*5,'Yticklabel',[0:0.04:0.24]*6/6*5);
     % ylim([0 0.04]);
elseif kk==7%IkBdProtein
      set(gca,'YLim',[0 0.012]*4,'Ytick',[0:0.002:0.02]*4,'Yticklabel',[0:0.2:2]);
elseif kk==8%IkBeProtein
     set(gca,'YLim',[0 0.018],'Ytick',[0:0.003:0.018],'Yticklabel',[0:0.2:1.2]);
elseif kk==6%IkBbProtein
      set(gca,'YLim',[0 0.024]/4*5*0.9,'Ytick',[0:0.004:0.04]/4*5*0.9,'Yticklabel',[0:0.004:0.04]*40/4*5);
elseif kk==10%A20Protein_new
     if id.Perturb(19)==1
       ylabel('\muM');
       ylim([0 4]);
   else
      set(gca,'YLim',[0 0.2]*1.5/1.15,'Ytick',[0:0.04:0.2]*1.5/1.15,'Yticklabel',[0:0.04:0.2]*10);
     end
elseif kk==9 %RIP1
      set(gca,'YLim',[0 0.04],'Ytick',[0:0.008:0.04],'Yticklabel',[0:0.008:0.04]*25);
elseif kk==2 %RIP3
      set(gca,'YLim',[0 0.2],'Ytick',[0:0.04:0.2],'Yticklabel',[0:0.04:0.2]*10);
end
 



%%
%plot data
if kk==4
    if Type~=2
        if id.Perturb(19)==1
        else
        errorbar(data.NuclearNFkB(:,1),data.NuclearNFkB(:,2)/15*1.2,data.NuclearNFkB(:,3)/15*1.2,'o','markersize',10,'Color','g','MarkerFaceColor','g'); hold on;
        end
    end
end
if kk==1
    if id.Perturb(19)==1
       
   else
    errorbar(data.A20mRNA_new(:,1),data.A20mRNA_new(:,2)/1e4/8*1.2,data.A20mRNA_new(:,3)/1e4/8*1.2,'o','markersize',10,'Color',[1 0.4 0.6],'MarkerFaceColor',[1 0.4 0.6]); hold on; 
    end
    
end


if Type==1 %&& id.Perturb(14)~=.6 % data for WT only
if kk==5
    errorbar(data.IkBaProtein(:,1),data.IkBaProtein(:,2)/6,data.IkBaProtein(:,3)/6,'o','markersize',10,'Color','g','MarkerFaceColor','g'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
elseif kk==3
   %plot(data.pMLKL_insoluble(:,1),data.pMLKL_insoluble(:,2)/2.5,'o','markersize',10,'Color','c','MarkerFaceColor','c');hold on;             
    errorbar(data.pMLKL_insoluble(:,1),data.pMLKL_insoluble(:,2)/2.5,data.pMLKL_insoluble(:,3)/2.5,'o','markersize',10,'Color','c','MarkerFaceColor','c');hold on;             
elseif kk==7
    errorbar(data.IkBTime(:,1),data.IkBdProtein(:,1)/100*4,data.IkBdProtein(:,2)/100*4,'o','markersize',10,'Color','g','MarkerFaceColor','g'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
 elseif kk==8
    errorbar(data.IkBTime(:,1),data.IkBeProtein(:,1)/70,data.IkBeProtein(:,2)/70,'o','markersize',10,'Color','g','MarkerFaceColor','g'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
  elseif kk==6
    errorbar(data.IkBTime(:,1),data.IkBbProtein(:,1)/40*0.9,data.IkBbProtein(:,2)/40*0.9,'o','markersize',10,'Color','g','MarkerFaceColor','g'); hold on;%,'Color',color(kk,:),'MarkerFaceColor',color(kk,:)); hold on;
  elseif kk==10
      if id.Perturb(19)==1  
        else
      errorbar(data.A20Protein_new(:,1),data.A20Protein_new(:,2)/10*1.3,data.A20Protein_new(:,3)/10*1.3,'o','markersize',10,'Color',[1 0.4 0.6],'MarkerFaceColor',[1 0.4 0.6]);hold on;                        
      end%plot(data.A20Protein_new(:,1)*60,data.A20Protein_new(:,2)/10,'o','markersize',10,'Color',[1 0.4 0.6],'MarkerFaceColor',[1 0.4 0.6]);hold on;                        
  elseif kk==9 %RIP1

  elseif kk==2 %RIP3
      errorbar(data.RIP3_insoluble(:,1),data.RIP3_insoluble(:,2)/max(data.RIP3_insoluble(:,2))/1e1,data.RIP3_insoluble(:,3)/max(data.RIP3_insoluble(:,2))/1e1,'o','markersize',10,'Color','c','MarkerFaceColor','c');hold on;                       
end
end

    
xlabel3=xlabel('time (h)','FontSize',20);
    xlim([0 tf/3600])
%if kk<3
set(gca,'XLim',[0 24],'Xtick',0:4:24,'Xticklabel',0:4:24);
%set(gca,'XLim',[0 id.sim_time*3/4],'Xtick',0:240*2:720*2,'Xticklabel',0:4:24)    
set(gca,'linewidth',2);
% yname=[Name(kk),'(A.U.)'];
% ylabel([string(Name(kk)),'(A.U.)'],'FontSize',20)
ylabel(sprintf('%s ',char(Name(kk)),'(A.U.)'),'FontSize',20)
 if id.Perturb(19)==1
      if  kk==10 || kk==1
          ylabel(sprintf('%s ',char(Name(kk)),'(\muM)'),'FontSize',20)
      end
 end
       
set(gca,'fontweight','bold','FontSize',20);
%set(gca,'FontSize',20);

% h=legend(Name(kk),'Location','northeast');
% set(h,'fontweight','bold','FontSize',20);



if kk==1
    figurename2=fig.figurename13;
elseif kk==4
    figurename2=fig.figurename14;
elseif kk==3
    figurename2=fig.figurename15;
    elseif kk==5
    figurename2=fig.figurename21;
elseif kk==6
    figurename2=fig.figurename22;
    elseif kk==7
    figurename2=fig.figurename23;
    elseif kk==8
    figurename2=fig.figurename24;
    elseif kk==9
    figurename2=fig.figurename25;
    elseif kk==10
    figurename2=fig.figurename26;
    elseif kk==11
    figurename2=fig.figurename27;
    elseif kk==2
    figurename2=fig.figurename28;
    elseif kk==12
    figurename2=fig.figurename29;
end
 saveas(gcf,figurename2);   
end




%%
%Plot single-cell simulations
c=jet(sizeScan);
SpeciesToPlot=[2,3,5,6];

for kk=[1 3 4 6]
%figure(1);clf;
figure('position', [-1600, 10, 600, 400])
if kk<=5
for kscan=1:sizeScan
    if kk==3
        plot(tt(1:ThresholdIndex(kscan)),max(DataSave{kscan,1}(kk,1:ThresholdIndex(kscan)),0),'-','Color',c(kscan,:),'linewidth',2);hold on;
        h2=plot(tt(ThresholdIndex(kscan):end),max(DataSave{kscan,1}(kk,ThresholdIndex(kscan):end),0),'-','Color',c(kscan,:),'linewidth',2);hold on;
        h2.Color(4)=0.1;  
    else
    plot(tt,max(DataSave{kscan,1}(kk,:),0),'-','Color',c(kscan,:),'linewidth',2);hold on;
    end
end
end

 if kk==3
         %set(gca,'YLim',[0 0.5],'Ytick',[0:0.1:0.5],'Yticklabel',[0:0.2:1]) ;
      set(gca,'YLim',[0 1],'Ytick',[0:0.2:1],'Yticklabel',[0:0.2:1]) ;
       %set(gca,'YLim',[0 .2],'Ytick',[0:0.05:0.2],'Yticklabel',[0:0.05:0.2]) ;
       % set(gca,'YLim',[0 .5],'Ytick',[0:0.1:.5],'Yticklabel',[0:0.1:.5]) ;
 elseif kk==4
        ylim([0 0.2]); 
 elseif kk==1
  % set(gca,'YLim',[0 100]/1e4/6*10,'Ytick',[0:20:100]/1e4/6*10,'Yticklabel',0:2:10);
     %set(gca,'YLim',[0 100]/1e4/8/3,'Ytick',[0:20:100]/1e4/8/3,'Yticklabel',0:2:10);
    set(gca,'YLim',[0 200]/1e4/8/3,'Ytick',[0:40:200]/1e4/8/3,'Yticklabel',0:4:20);
 end      

 
 
if kk==3
    plot(tt,Threshold*ones(1,length(tt)),'-','Color','k','linewidth',2);hold on;
end
if kk==4
    ylim([0 0.2]); 
end

if kk==6
   %histogram(FirstPassage,15,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5]); 
   % histogram(FirstPassage,15,'FaceColor','b','EdgeColor','b'); 
   edges = [0:2:24];
%     h = histogram(FirstPassage,edges,'FaceColor','b','EdgeColor','b');  
%     ylabel1=ylabel('cell count','FontSize',20);
%     ylim([0 150]);
if id.Perturb(19)==1
    h1 = histogram(FirstPassage,edges,'FaceColor','b','EdgeColor','b');  
    ylabel1=ylabel('death cell count','FontSize',20);
    ylim([0 50]);   
else
     h1 = histogram(FirstPassage,edges,'FaceColor','b','EdgeColor','b','Normalization','probability');  
    ylabel1=ylabel('death event frequency','FontSize',20);
    ylim([0 0.5]);   
end
end
xlabel3=xlabel('time (h)','FontSize',20);
    xlim([0 tf/3600])
%if kk<3
set(gca,'XLim',[0 24],'Xtick',0:4:24,'Xticklabel',0:4:24);
%set(gca,'XLim',[0 id.sim_time*3/4],'Xtick',0:240*2:720*2,'Xticklabel',0:4:24)    
set(gca,'linewidth',2);
set(gca,'fontweight','bold','FontSize',20);
%set(gca,'FontSize',20);
if kk==6
   % disp(round(length(FirstPassage)/sizeScan*100));
    h=text(0.6,0.9,['survival ',num2str(round(100-length(FirstPassage)/sizeScan*100)),'%'],'Units','normalized');set(h,'fontweight','bold','FontSize',20);
%h=legend(['survival ',num2str(round(100-length(FirstPassage)/sizeScan*100)),'%'],'Location','northeast');
else
%h=legend(Name(kk),'Location','northeast');
% yname=[Name(kk),'(A.U.)'];
% ylabel(yname,'FontSize',20)
ylabel(sprintf('%s ',char(Name(kk)),'(A.U.)'),'FontSize',20)


end

%set(gca,'linewidth',2);
%set(gca,'FontSize',20);
if kk==6
    figurename2=fig.figurename5;
elseif kk==1
    figurename2=fig.figurename6;
elseif kk==4
    figurename2=fig.figurename8;
elseif kk==3
    figurename2=fig.figurename7;
%figurename2=[filefolder,'\Bimodal',num2str(index),'_',num2str(Threshold),'_',num2str(kk),'.jpg'];
end
 saveas(gcf,figurename2);   
end

DeathTimeTemp=sort(FirstPassage);
n=1;
tTemp(n)=0;DeathRatesPlot(n)=0;
Toplot=1;
for tt3=0:0.1:24-5
    tt2=tt3+5;
    DeathInWindow=find(DeathTimeTemp>tt3 & DeathTimeTemp<tt2);
      assignin('base', 'deathInWindow', DeathInWindow);
    if ~isempty(DeathInWindow)
    if (sizeScan-DeathInWindow(1))>100
    n=n+1;
    DeathRatesPlot(n)=length(DeathInWindow)/(sizeScan-DeathInWindow(1))/5;
    tTemp(n)=tt3;
    
    end
    else
        Toplot=0;
    end
end
%if Toplot==1 
figure('position', [-1600, 10, 600, 400])
plot(tTemp,DeathRatesPlot,'Color','k','linewidth',2);hold on;
ylim([0 0.2]);
xlim([0 24]);
xlabel('time (h)','FontSize',20);
ylabel('death rates','FontSize',20);
set(gca,'linewidth',2);
set(gca,'fontweight','bold','FontSize',20);
figurename2=fig.figurename10;
saveas(gcf,figurename2);   
%end


%%
%plot heat map of WT
%Need to revise...
if fig.PlotHeatMapWT==1
    figure('position', [-1600, 10, 800, 200])
    h1 = histogram(FirstPassage,edges,'FaceColor','b','EdgeColor','b');  
        surf(edges,[0 1],[[h1.Values,0];[h1.Values,0]]);
        assignin('base', 'hh1', h1);
        h=colorbar; 
        ylabel(h, 'cell count') 
        %ylabel(h, 'Frequency') 
        shading interp;
        %shading flat;
        view(2);
        if fig.scanA20Constituive==0
        title('WT');
        elseif fig.scanA20Constituive==2
             title('A20 KO');
        end
        caxis([0 50]);
        colormap jet;
        %ylabel('constitutive A20 (2^X)')
        set(gca,'ytick',[]);
         xlabel('time (h)')
        %ylim([min(FoldExponent), max(FoldExponent)]);
        xlim([min(edges(1:end)), max(edges(1:end))]);
         set(gca,'FontSize',20);
        saveas(gcf,fig.figurename123);
        
        figure('position', [-1600, 10, 800, 200])
        %surf(NFkBSummary1(1:100:end),[0 1],[NFkBSummary2(1:100:end);NFkBSummary2(1:100:end)]);
        surf(NFkBSummary1,[0 1],[NFkBSummary2;NFkBSummary2]);
        assignin('base', 'NFkBSummary2', NFkBSummary2);
        %surf(edges(1:end-1),[0 1],[DeathTimeSummary(end,:);DeathTimeSummary(end,:)]);
        h=colorbar; 
        ylabel(h, 'NF\kappaB');
        colormap jet;
        caxis([0 0.08]);
        %h.Ticks = linspace(0, max(NFkBSummary2), 3) ; %Create 8 ticks from zero to 1
        h.Ticks = linspace(0, 0.08, 5) ; %Create 8 ticks from zero to 1
       % set(gca,'YLim',[0 0.08],'Ytick',[0:0.02:0.08],'Yticklabel',[0:0.3:1.2]);
        h.TickLabels = num2cell([0:0.3:1.2]) ;
        assignin('base', 'h', h);
        %ylabel(h, 'Frequency') 
        shading interp;
       % shading flat;
        view(2);
        if fig.scanA20Constituive==0
            title('WT');
        elseif fig.scanA20Constituive==2
             title('A20 KO');
        end
        %ylabel('constitutive A20 (2^X)')
        set(gca,'ytick',[]);
         xlabel('time (h)')
        %ylim([min(FoldExponent), max(FoldExponent)]);
        xlim([min(edges(1:end)), max(edges(1:end))]);
         set(gca,'FontSize',20);
        saveas(gcf,fig.figurename124);
end

%%
% c=jet(sizeScan);
% %figure(1);clf;
% figure('position', [100, 10, 1600, 800])
%     
% for kk=1:6
%     %figure(kk),clf,
% subplot(2,3,kk), hold on
% if  kk<6
% for kscan=1:sizeScan
% %     disp(size(DataSave{kscan,1}(kk,:)));disp(size(tt));
% plot(tt,max(DataSave{kscan,1}(kk,:),0),'-','Color',c(kscan,:),'linewidth',2);hold on;
% end
%     if kk==3
%         plot(tt,Threshold*ones(1,length(tt)),'-','Color','k','linewidth',2);hold on;
%         %ylim([0 Threshold*2])
%     end
% %     if kk==4 && Para.InitialValue(3)~=0
% %         ylim([0 Para.InitialValue(3)*2])
% %     end
% end
% xlabel3=xlabel('time (h)','FontSize',20);
% xlim([0 tf])
% ylabel1=ylabel('Activity','FontSize',20);
% 
% if kk==6
%    hist(FirstPassage);
%    xlabel3=xlabel('time (h)','FontSize',20);
%     xlim([0 tf])
%     ylabel1=ylabel('cell count','FontSize',20);
%     ylim([0 150])
% end
% set(gca,'FontSize',20);
% h=legend(Name(kk),'Location','northeast');
% set(h,'FontSize',20);
% set(gca,'linewidth',2);
% set(gca,'FontSize',20);
% 
% end
% 
%  saveas(gcf,fig.figurename2);      
 
filename=[filefolder,'\DeathTime',num2str(Threshold),'_',num2str(index2),'.mat'];
%save(fig.figurename9,'Para','DataSave','tt','tf','Name','TNFPlot','Threshold','FirstPassage');
%dddd
%'edges','FirstPassage','tTemp','DeathRatesPlot'
if ToSave==1
save(fig.figurename9,'Para','data','DataSave','tt','tf','IndexToShow','sizeScan','IDs','index2','fig','Name','TNFPlot','Threshold','h1','edges','FirstPassage','tTemp','DeathRatesPlot','ThresholdIndex');
end
end