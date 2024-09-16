%clc; close all;

%% Run PLA over all estimated parameters (Nothing changes)
load('p_opt.mat');

%%
thetaN = length(p_opt);
upPLval = [];
upPLres = [];
downPLval = [];
downPLres = [];

nr_parameters = 6;
lsq_options = optimset('Algorithm','trust-region-reflective');

lb = [ 0   0   0    0    0   0   ];
ub = [1e6 1e6 1e6 2 2 2];

%%
for i = 1:nr_parameters
    if i <4
        [upPLval(i,:), upPLres(i,:)] = PLA(@(p_opt)PLA_cost_func(p_opt),p_opt',i,ub(i),[],lb,ub,lsq_options,0.0005,0.002,[],[],75);
        [downPLval(i,:), downPLres(i,:)] = PLA(@(p_opt)PLA_cost_func(p_opt),p_opt',i,lb(i),[],lb,ub,lsq_options,0.015,0.02,[],[],50);
    else
        [upPLval(i,:), upPLres(i,:)] = PLA(@(p_opt)PLA_cost_func(p_opt),p_opt',i,ub(i),[],lb,ub,lsq_options,0.0125,0.0175,[],[],75);
        [downPLval(i,:), downPLres(i,:)] = PLA(@(p_opt)PLA_cost_func(p_opt),p_opt',i,lb(i),[],lb,ub,lsq_options,0.015,0.02,[],[],50);
    end
   
end
figure
%%


labels = {'sTNF','sIL6','sIL10','kTNFmRNA','kIL6mRNA','kIL10mRNA'};

for i = 1:nr_parameters
    subplot(2,3,i)
    semilogx(downPLval(i,:),downPLres(i,:),'k')
    xlabel('Value')
    ylabel('Residual error')
    hold on
    semilogx(upPLval(i,:),upPLres(i,:),'k')
    hold on
    semilogx([p_opt(i) p_opt(i)],[min([upPLres(i,:) downPLres(i,:)]) max([upPLres(i,:) downPLres(i,:)])],'k--')
    hold on
    semilogx(upPLval(i,1),upPLres(i,1),'x','LineWidth',2,'Color',[8 144 153]/256)
    ylim([0 0.025])
    title(labels(i))
    grid on
end