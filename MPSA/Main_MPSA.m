%Template for Multi-Parametric Sensitivity Analysis with Try-Catch.
%Try-Catch can be used to prevent that a Matlab error stops execution of the program.
%Sampling can generate unfavorable combinations of values for the model
%parameter, which can result in a model that numerically unstable
%(difficult to simulate; a case of ill-posedness).
%To prevent that a numerical issue generates a Matlab error which will
%terminate execution of the program it is recommended to use a try ? catch
%implementation.
%Also see: https://nl.mathworks.com/help/matlab/ref/try.html .

%Natal van Riel, TU/e
%Freek Relouw, TU/e

clear all
close all
tic

%Nominal parameter values
% Cytokine scaling (4)
sTNF = 14217.8909603617;    %estimated from data
sIL6 = 4879.44319226363;    %estimated from data
sIL10 = 23718.6860292186;   %estimated from data
sIL1 = 450;

% mRNA elimination (4)
kTNFmRNA = 1.0479;      %estimated from data
kIL6mRNA = 1.6005;      %estimated from data
kIL10mRNA= 1.5038;      %estimated from data
kIL1mRNA = 0.38*0.693;

% Cytokine elimination (4)
kTNF = 2*0.693;
kIL6 = 6*0.693;
kIL10= 30*0.693;
kIL1 = 17.85*0.693;

% Hill equation parameters (18)
xLPS_TNF = 13.8;
nLPS_TNF = 1;
xLPS_IL6 = 0.0006;
nLPS_IL6 = 0.4;
xLPS_IL10 = 29.95;
nLPS_IL10 = 0.64;
xLPS_IL1 = 0.85;
nLPS_IL1 = 0.98;
xIL10_TNF = 0.6;
nIL10_TNF = 0.9;
xIL10_IL6 = 15;
nIL10_IL6 = 0.3;
xIL10_IL1 = 0.01;
nIL10_IL1 = 2;
xIL10_IL10 = 100;
nIL10_IL10 = 2;
xIL6_Temp = 3000;
nIL6_Temp = 2;

% Heart rate (2)
kHRBP = 0.5;
kHRTemp = 15;

% Differential equation parameters
LPS_dose = 165;
kM = 0.33;
kM1 = 0.1;
nTemp = 2.7;
sTempIL1 = 0.01;
kTemp = 0.9;
sTempIL6 = 5;
TempMax = 42;
sD = 0.0001;
kD = 0.01;

tmaxIL6 = 2;
ntIL6 = 2;
tmaxIL10 = 2.5;
ntIL10 = 3;
tmaxIL1 = 1;
ntIL1 = 1;

dummy1 = 0.0001;
dummy2 = 0.001;
dummy3 = 0.01;
dummy4 = 0.1;
dummy5 = 1;
dummy6 = 10;
dummy7 = 100;
dummy8 = 1000;
dummy9 = 10000;
dummy10 = 100000;
dummy11 = 1000000;
dummy12 = 10000000;


%Nominal parameter array
par0 = [sTNF, sIL6, sIL10, sIL1,...
    kTNFmRNA, kIL6mRNA, kIL10mRNA, kIL1mRNA, ...
    kTNF, kIL6, kIL10, kIL1, ...
    xLPS_TNF, nLPS_TNF, xLPS_IL6, nLPS_IL6, xLPS_IL10, nLPS_IL10, xLPS_IL1, nLPS_IL1, ...
    xIL10_TNF, nIL10_TNF, xIL10_IL6, nIL10_IL6, xIL10_IL1, nIL10_IL1, xIL10_IL10, nIL10_IL10, ...
    xIL6_Temp, nIL6_Temp, ...
    kHRBP, kHRTemp, ...
    LPS_dose, kM, kM1, nTemp, sTempIL1, kTemp, sTempIL6, TempMax, sD, kD, ...
    tmaxIL6, ntIL6, tmaxIL10, ntIL10, tmaxIL1, ntIL1, dummy1,dummy2,dummy3,dummy4,dummy5,...
    dummy6,dummy7,dummy8,dummy9,dummy10,dummy11,dummy12
    ];



% Simulate default model for nominal parameter values par0
% and calculate feature from selected output(s)
feature0 = PSA_model_code(5, par0, 0, 0, 0, 37, 90, 70)';  %Simulate is a function you have to write yourself

Np = 48;         % no. of parameters (rate constants, initial conditions) included in analysis
Nd = 12;          % no. of dummies
n_MPSA = 5000;    % no. of parameter ensembles / no. of Monte Carlo simulations


% Make the V matrix
V = NaN(length(feature0),n_MPSA);

% Latin Hypercube Sampling of parameter space
scale = lhsdesign(n_MPSA, Np+Nd);   % random uniform distributed parameter sets (values [0,1])
% (scale is an n_MPSA x Np+Nd matrix)
% Use the scale matrix to generate MC parameter sets relative to reference values par0:
% e.g. 0.1*par0(i) < par_MC(i,j) < 10*par0(i)
% The Matlab function 'rescale' could be used.

par_MC = zeros(n_MPSA, Np+Nd);
for i = 1:Np+Nd
    par_MC(:,i) = rescale(scale(:,i),par0(i)-.5*par0(i),par0(i)+.5*par0(i));
end

%% Monte Carlo simulations
for j = 1 : n_MPSA
    partemp=par_MC(j,1:Np); %select the Np model parameters from set j
    % Simulate model for partemp and calculate feature of interest from
    % selected model outputs
    try %if simulation works, the sensitivity criterion for parameter set j gets a value
        feature = PSA_model_code(5,partemp,0,0,0,37,90,70)';
        % Calculate sensitivity criterion
        % Sum of squared differences between the perturbed and reference output as example
        V(:,j) = (feature0-feature).^2;
    catch %if simulation fails, the sensitivity criterion for parameter set j is marked as NaN
        V(:,j)= NaN;
    end
    fprintf('iteration %d completed \n',j)
end
toc

%% Preprocessing
V = real(V);

%%
var_interest = 4:15;
nr_rep = 50;

S = NaN(length(var_interest),Np+Nd);

%% Bootstrap
n_temp_array = [50 100 200 300 400 500 600 700 800 900 1000 2000 3000 4000];% 5000 7500 10000 12500];
collectAverageDummy = [];
collectStdDummy = [];
x = 1:Np;
d = Np+1:Np+Nd;

for u = var_interest 
    for m = 1:length(n_temp_array)
        n_temp = n_temp_array(m);
        K_S = [];

        for k = 1:nr_rep
            bootstrap = randperm(n_MPSA,n_temp);
            V_temp = V(u,bootstrap);
            par_MC_temp = par_MC(bootstrap,:);


            %% Classification of acceptable vs unacceptable simulations
            flag  = zeros(n_temp, 1);
            Sa    = zeros(n_temp, Np);
            Su    = zeros(n_temp, Np);
            value = zeros(n_temp, Np);

            threshold = mean(V_temp);    %mean as threshold
            acc       = find(V_temp <= threshold);
            unacc     = find(V_temp  > threshold);
            flag(acc) = 1;

            %% Cumulative distributions (for model parameters and dummies)
            for i = 1 : Np+Nd
                temp = [par_MC_temp(:,i), flag];     %associate 1 to acceptable cases and 0 to unacceptable parameter values
                temp = sortrows(temp,1);        %sorts temp based on column 1

                value(:,i) = temp(:,1);
                Sa(:,i) = cumsum(temp(:,end));
                Su(:,i) = cumsum(-1*temp(:,end)+ones(n_temp,1));
                Sa(:,i) = Sa(:,i)/max(Sa(:,i));
                Su(:,i) = Su(:,i)/max(Su(:,i));
            end

            %% Kolmogorov-Smirnov statistic
            K_S(k,:) = max(abs(Sa-Su));

        end
        averageKS = mean(K_S);
        stdKS = std(K_S);

        [M,I] = max(averageKS(:,d));

        collectAverageDummy = [collectAverageDummy M];
        collectStdDummy = [collectStdDummy stdKS(:,Np+I)];
    end
    S(u-3,:) = averageKS;
end


%% Display K_S
% h = figure;
% h = bar(x,diag(averageKS(:,x)),'stacked');
% grid on
% 
% %Perform Z test to find parameters that are significantly  bigger than
% %biggest dummy variable.
% KSmean = averageKS(:,x);
% KSstd = stdKS(:,x);
% Dummymean = averageKS(:,Np+I);
% Dummystd = stdKS(:,Np+I);
% Z = (KSmean-Dummymean)./sqrt( (KSstd.^2./nr_rep) + (Dummystd^2/nr_rep) );
% index = find(Z <= 2.33);
% 
% 
% set(h(x),'facecolor',[0 0.4470 0.7410])
% set(h(index),'facecolor',[1 1 1])
% hold on
% 
% er = errorbar(x,averageKS(:,x),2*stdKS(:,x),2*stdKS(:,x));
% er.Color= [0 0 0];
% er.LineStyle = 'none';
% 
% ylabel('K_S distance')
% title('Multi parametric sensitivity analysis')
% xticks(1:length(par0))
% xticklabels({'sTNF','sIL6','sIL10','sIL1','kTNFmRNA','kIL6mRNA','kIL10mRNA','kIL1mRNA,'...
%     'kTNF','kIL6','kIL10','kIL1','xLPS_{TNF}','nLPS_{TNF}','xLPS_{IL6}','nLPS_{IL6}','xLPS_{IL10}',...
%     'nLPS_{IL10}','xLPS_{IL1}','nLPS_{IL1}','xIL10_{TNF}','nIL10_{TNF}','xIL10_{IL6}','nIL10_{IL6}',...
%     'xIL10_{IL1}','nIL10_{IL1}','kHR_{BP}','kHR_{Temp}'})
% xtickangle(45)
% 
% % plot line with biggest dummy variable value + std
% hold on
% plot([0 Np+1],[averageKS(:,Np+I) averageKS(:,Np+I)],'--r',[0 Np+1],[averageKS(:,Np+I)-2*stdKS(:,Np+I) averageKS(:,Np+I)-2*stdKS(:,Np+I)],':r',[0 Np+1],[averageKS(:,Np+I)+2*stdKS(:,Np+I) averageKS(:,Np+I)+2*stdKS(:,Np+I)],':r')


%% Heatmap representation
N ={'sTNF','sIL6','sIL10','sIL1','kTNFmRNA','kIL6mRNA','kIL10mRNA','kIL1mRNA,'...
    'kTNF','kIL6','kIL10','kIL1','xLPS_{TNF}','nLPS_{TNF}','xLPS_{IL6}','nLPS_{IL6}','xLPS_{IL10}',...
    'nLPS_{IL10}','xLPS_{IL1}','nLPS_{IL1}','xIL10_{TNF}','nIL10_{TNF}','xIL10_{IL6}','nIL10_{IL6}',...
    'xIL10_{IL1}','nIL10_{IL1}','kHR_{BP}','kHR_{Temp}'};
n = length(N);
M ={'TNFmRNA','TNF','IL6mRNA','IL6','IL10mRNA','IL10','IL1mRNA','IL1','Temp','BP','HR','D'};
m = length(M);
imagesc(S(1:end,1:n)); % plot the matrix
xticks(1:n);
xticklabels(N);
yticks(1:m);
yticklabels(M);
colorbar