function [t,TNF,IL6,IL10,IL1,t2,Temp,BP,HR,color,corrTNF,corrIL6,corrIL10] = data_model(i,data)
% Assign data

if i == 1
    t = data.Dillingh2014bolus05_cytokines(:,1);
    TNF = data.Dillingh2014bolus05_cytokines(:,2);
    IL6 = data.Dillingh2014bolus05_cytokines(:,3);
    IL10 = zeros(length(t),1);
    IL1 = zeros(length(t),1);
    t2= data.Dillingh2014bolus05_vitals(:,1);
    Temp = data.Dillingh2014bolus05_vitals(:,2)+36.5;
    HR = data.Dillingh2014bolus05_vitals(:,4)+70;
    BP = zeros(length(t2),1);
    color = [255 183 5]/256;
    corrTNF=1;
    corrIL6=1;
    corrIL10=0;
elseif i == 2
    t = data.Dillingh2014bolus1_cytokines(:,1);
    TNF = data.Dillingh2014bolus1_cytokines(:,2);
    IL6 = data.Dillingh2014bolus1_cytokines(:,3);
    IL10 = zeros(length(t),1);
    IL1 = zeros(length(t),1);
    t2= data.Dillingh2014bolus1_vitals(:,1);
    Temp = data.Dillingh2014bolus1_vitals(:,2)+36.5;
    HR = data.Dillingh2014bolus1_vitals(:,4)+70;
    BP = zeros(length(t2),1)+90;
    color = [150 199 217]/256;
    corrTNF=1;
    corrIL6=1;
    corrIL10=0;
elseif i == 3
    t = data.Dillingh2014bolus2_cytokines(:,1);
    TNF = data.Dillingh2014bolus2_cytokines(:,2);
    IL6 = data.Dillingh2014bolus2_cytokines(:,3);
    IL10 = zeros(length(t),1);
    IL1 = zeros(length(t),1);
    t2= data.Dillingh2014bolus2_vitals(:,1);
    Temp = data.Dillingh2014bolus2_vitals(:,2)+36.5;
    HR = data.Dillingh2014bolus2_vitals(:,4)+70;
    BP = zeros(length(t2),1)+90;
    color = [28 72 109]/256;
    corrTNF=1;
    corrIL6=1;
    corrIL10=0;
elseif i == 4
    t = data.Taudorf2007bolus03_cytokines(:,1);
    TNF = data.Taudorf2007bolus03_cytokines(:,2);
    IL6 = data.Taudorf2007bolus03_cytokines(:,3);
    IL10 = zeros(length(t),1);
    IL1 = zeros(length(t),1);
    t2= data.Taudorf2007bolus03_vitals(:,1);
    Temp = data.Taudorf2007bolus03_vitals(:,2);
    HR = data.Taudorf2007bolus03_vitals(:,4);
    BP = zeros(length(t2),1)+90;
    color = [205 14 5]/256;
    corrTNF=1;
    corrIL6=1;
    corrIL10=0;
elseif i == 5
    t = data.Kiers2017bolus1_cytokines(:,1);
    TNF = data.Kiers2017bolus1_cytokines(:,2);
    IL6 = data.Kiers2017bolus1_cytokines(:,3);
    IL10 = data.Kiers2017bolus1_cytokines(:,4);
    IL1 = data.Kiers2017bolus1_cytokines(:,5);
    t2 = data.Kiers2017bolus1_vitals(:,1);
    Temp = data.Kiers2017bolus1_vitals(:,2);
    BP = data.Kiers2017bolus1_vitals(:,3);
    HR = data.Kiers2017bolus1_vitals(:,4);
    color = [1 31 75]/256;
    corrTNF=1;
    corrIL6=1;
    corrIL10=1;
elseif i == 6
    t = data.Kiers2017bolus2_cytokines(:,1);
    TNF = data.Kiers2017bolus2_cytokines(:,2);
    IL6 = data.Kiers2017bolus2_cytokines(:,3);
    IL10 = data.Kiers2017bolus2_cytokines(:,4);
    IL1 = data.Kiers2017bolus2_cytokines(:,5);
    t2 = data.Kiers2017bolus1_vitals(:,1);
    Temp = data.Kiers2017bolus2_vitals(:,2);
    BP = data.Kiers2017bolus2_vitals(:,3);
    HR = data.Kiers2017bolus2_vitals(:,4);
    color = [0.9290 0.6940 0.1250];
    corrTNF=1;
    corrIL6=1;
    corrIL10=1;
elseif i == 7
    t = data.Kiers2017threehour1_cytokines(:,1);
    TNF = data.Kiers2017threehour1_cytokines(:,2);
    IL6 = data.Kiers2017threehour1_cytokines(:,3);
    IL10 = data.Kiers2017threehour1_cytokines(:,4);
    IL1 = data.Kiers2017threehour1_cytokines(:,5);
    t2 = data.Kiers2017threehour1_vitals(:,1);
    Temp = data.Kiers2017threehour1_vitals(:,2);
    BP = data.Kiers2017threehour1_vitals(:,3);
    HR = data.Kiers2017threehour1_vitals(:,4);
    color = [255 183 5]/256;
    corrTNF=1;
    corrIL6=1;
    corrIL10=1;
else
    t = 0;
    TNF = 0;
    IL6 = 0;
    IL10 = 0;
    IL1 = zeros(length(t),1);
    t2=t;
    Temp = [36.5];
    BP = [90];
    HR = [65];
    color = [0 0.4470 0.7410];
    corrTNF=0;
    corrIL6=0;
    corrIL10=0;
end
end