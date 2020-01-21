
nboot = 500;
plotIndex=1;
figure;
p_value_matrix=zeros(2,3);
gmmodelArray={};
for i=1:3
    fileToLoad=strcat('arrayOfDeathTimes_T',num2str(i),'.mat');
    load(fileToLoad);
    %arrayOfDeathTimes=(arrayOfDeathTimes/1.5)-30;
    tnfSample=sort(arrayOfDeathTimes);
    
    [xpdf, n, b] = compute_xpdf(tnfSample);
    [dip, p_value, xlow, xup] = HartigansDipSignifTest(xpdf, nboot); 
    p_value_matrix(1,i)=p_value;
    subplot(3,2,plotIndex);
    
    
    bar(b,n,'FaceColor',[0.5,0.5,0.5],'EdgeColor',[0,0,0])
    set(gca,'xticklabel',[2:2:24])
    set(gcf,'color','w');
    xlim([0,1550]); 
    %xlim('auto'); 
    options = statset('MaxIter',1000000);
    scoresTNF=[0,0,0,0];
    for terms=1:4
        gmmodel=fitgmdist(tnfSample',terms,'Options',options);
        AICScore=gmmodel.BIC;
        scoresTNF(terms)=AICScore;
    end
    [~,bestModes]=min(scoresTNF);
    title({sprintf('Probability of unimodal %.2f', p_value),sprintf('lowest AIC:%s',num2str(bestModes))})
    gmmodel=fitgmdist(tnfSample',2,'Options',options)
    gmmodelArray{end+1}=gmmodel;
    hold on;%plot([60:960+59],pdf(gmmodel,[1:960]'),'k','linewidth',2);hold on;
    plot([60+30:1440+59+30],gmmodel.ComponentProportion(1)*normpdf([1:1440],gmmodel.mu(1),sqrt(gmmodel.Sigma(1))),'r','linewidth',2);
    plot([60+30:1440+59+30],gmmodel.ComponentProportion(2)*normpdf([1:1440],gmmodel.mu(2),sqrt(gmmodel.Sigma(2))),'b','linewidth',2);
    
    
    
    
    
    plotIndex=plotIndex+1;
    fileToLoad=strcat('arrayOfDeathTimes_CT',num2str(i),'.mat');
    load(fileToLoad);
    %arrayOfDeathTimes=arrayOfDeathTimes*1.5;

    CTSample=sort(arrayOfDeathTimes);

    [xpdf, n, b] = compute_xpdf(CTSample);
    [dip, p_value, xlow, xup] = HartigansDipSignifTest(xpdf, nboot); 
    p_value_matrix(2,i)=p_value;
    subplot(3,2,plotIndex);
    bar(b,n,'FaceColor',[0.5,0.5,0.5],'EdgeColor',[0,0,0])
    title(sprintf('Probability of unimodal %.2f', p_value))
    set(gca,'xticklabel',[2:2:24]);
    set(gcf,'color','w');

    xlim([0,1550]); 

        options = statset('MaxIter',1000000);
    scoresCT=[0,0,0,0];
        plotIndex=plotIndex+1;
    %let's fit 2X gaussian distrbution to these!
    %print -dpng modality.png
end
TNFPVals=p_value_matrix(1,:);
CTPVals=p_value_matrix(2,:);
meanTNF=mean(TNFPVals);
stdTNF=std(TNFPVals);
meanCT=mean(CTPVals);
stdCT=std(CTPVals);
titleString={strcat('TNF Unimodel Prob:',num2str(round(meanTNF,2)),'+/-',num2str(round(stdTNF,2))),strcat('Cx+TNF Unimodel Prob:',num2str(round(meanCT,2)),'+/-',num2str(round(stdCT,2)))};
suptitle(titleString);
print -dpng histsOutput.png

figure;
meanUnimodalStatsTNF=meanTNF;
stdUnimodalStatsTNF=stdTNF;
meanUnimodalStatsCx=meanCT;
stdUnimodalStatsCx=stdCT;

bar([1,2],[meanUnimodalStatsTNF,meanUnimodalStatsCx],'FaceColor',[0.5,0.5,0.5],'EdgeColor',[0,0,0],'LineWidth',2);
hold on;
capsize=0.3;
linewidth=2;
set(gcf,'color','w');

plot([1,1],[meanUnimodalStatsTNF,meanUnimodalStatsTNF+stdUnimodalStatsTNF],'k','linewidth',linewidth);
plot([1-capsize,1+capsize],[meanUnimodalStatsTNF+stdUnimodalStatsTNF,meanUnimodalStatsTNF+stdUnimodalStatsTNF],'k','linewidth',linewidth);
plot([2,2],[meanUnimodalStatsCx,meanUnimodalStatsCx+stdUnimodalStatsCx],'k','linewidth',linewidth);
plot([2-capsize,2+capsize],[meanUnimodalStatsCx+stdUnimodalStatsCx,meanUnimodalStatsCx+stdUnimodalStatsCx],'k','linewidth',linewidth);
xlim([0.5,2.5]);
ylim([0,1.1]);

[~,pval]=ttest2([0.44,0.34,0.13],[0.91,0.96,0.98]);
title(sprintf('pVal:%f',pval));


mode1=[gmmodelArray{1}.mu(2), gmmodelArray{2}.mu(2),gmmodelArray{3}.mu(2)];
mode2=[gmmodelArray{1}.mu(1), gmmodelArray{2}.mu(1),gmmodelArray{3}.mu(1)];
mode1mean=mean(mode1)
mode2mean=mean(mode2)
mode1STD=std(mode1)
mode2STD=std(mode2);
stdArray=[gmmodelArray{1}.Sigma(2),gmmodelArray{1}.Sigma(1);...
    gmmodelArray{2}.Sigma(2),gmmodelArray{2}.Sigma(1);...
    gmmodelArray{3}.Sigma(2),gmmodelArray{3}.Sigma(1)];
mode1stdmean=mean(stdArray(:,1))
mode1stdstd=std(stdArray(:,1))
mode2stdmean=mean(stdArray(:,2))
mode2stdstd=std(stdArray(:,2))

mixingProportionMode1=[gmmodelArray{1}.ComponentProportion(2),gmmodelArray{2}.ComponentProportion(2),gmmodelArray{3}.ComponentProportion(2)];
mixingProportionMode2=[gmmodelArray{1}.ComponentProportion(1),gmmodelArray{2}.ComponentProportion(1),gmmodelArray{3}.ComponentProportion(1)];
mode1mean=mean(mixingProportionMode1)
mode2mean=mean(mixingProportionMode2)
mode1STD=std(mixingProportionMode1)
mode2STD=std(mixingProportionMode2)