clear all
close all

%% LPSA
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


%Nominal parameter array
param0 = [sTNF, sIL6, sIL10, sIL1,... 
    kTNFmRNA, kIL6mRNA, kIL10mRNA, kIL1mRNA, ...
    kTNF, kIL6, kIL10, kIL1, ...
    xLPS_TNF, nLPS_TNF, xLPS_IL6, nLPS_IL6, xLPS_IL10, nLPS_IL10, xLPS_IL1, nLPS_IL1, ...
    xIL10_TNF, nIL10_TNF, xIL10_IL6, nIL10_IL6, xIL10_IL1, nIL10_IL1, xIL10_IL10, nIL10_IL10, ...
    xIL6_Temp, nIL6_Temp, ...
    kHRBP, kHRTemp, ...
    LPS_dose, kM, kM1, nTemp, sTempIL1, kTemp, sTempIL6, TempMax, sD, kD, ...
    tmaxIL6, ntIL6, tmaxIL10, ntIL10, tmaxIL1, ntIL1
];

% Factor change from nominal value
factor = 1.5;

% Initialize sensitivity matrix S
S = zeros(15, length(param0));  % Assuming there are 15 variables

for variable = 1:15

    % Simulation with nominal values
    AUC = PSA_model_code(5, param0, 0, 0, 0, 37, 90, 70);
    SV_nominal = real(AUC(:, variable));  % SV_nominal is 801x1

    % Sensitivity analysis with changed parameter values
    for i = 1:length(param0)
        % Create a copy of the nominal parameters and modify the i-th parameter
        param1 = param0;
        param1(i) = param0(i) * factor;  % Change the i-th parameter by the factor
        
        % Run the simulation with modified parameter set
        AUC = PSA_model_code(5, param1, 0, 0, 0, 37, 90, 70);
        
        % Get the sensitivity result for the i-th parameter
        SV_sensitivity = real(AUC(:, variable));  % SV_sensitivity should be 801x1
        
        % Calculate sensitivity (element-wise difference, normalized by SV_nominal)
        S(variable, i) = mean(abs(SV_sensitivity - SV_nominal) ./ SV_nominal,'omitnan');  % Mean of relative changes
    end
end

%% Heatmap representation
N = {'sTNF', 'sIL6', 'sIL10', 'sIL1', 'kTNFmRNA', 'kIL6mRNA', 'kIL10mRNA', 'kIL1mRNA', ...
    'kTNF', 'kIL6', 'kIL10', 'kIL1', 'xLPS_{TNF}', 'nLPS_{TNF}', 'xLPS_{IL6}', 'nLPS_{IL6}', ...
    'xLPS_{IL10}', 'nLPS_{IL10}', 'xLPS_{IL1}', 'nLPS_{IL1}', 'xIL10_{TNF}', 'nIL10_{TNF}', ...
    'xIL10_{IL6}', 'nIL10_{IL6}', 'xIL10_{IL1}', 'nIL10_{IL1}', 'kHR_{BP}', 'kHR_{Temp}'};

M = {'TNFmRNA', 'TNF', 'IL6mRNA', 'IL6', 'IL10mRNA', 'IL10', 'IL1mRNA', 'IL1\beta', 'Temp', 'BP', 'HR', 'D'};

% Create the heatmap
imagesc(S(4:15,1:length(N)));  % Plot the matrix for variables 4 to 15
xticks(1:length(N));
xticklabels(N);
yticks(1:length(M));
yticklabels(M);
colorbar;

%% Individual variable representation
% TNFmRNA   = 4     TNF  = 5
% IL6mRNA   = 6     IL6  = 7
% IL10mRNA  = 8     IL10 = 9
% IL1mRNA   = 10    IL1  = 10
% Temp      = 12    BP   = 13
% HR        = 14    D    = 15

% var_interest = 7;
% 
% 
% figure
% b = bar(S(var_interest,:));
% grid on
% ylabel('Sensitivity')
% %xlabel('Parameter')
% title('Local parametric sensitivity analysis')
% xticks(1:length(param0))
% xticklabels({'sTNF','sIL6','sIL10','sIL1','kTNFmRNA','kIL6mRNA','kIL10mRNA','kIL1mRNA,'...
%     'kTNF','kIL6','kIL10','kIL1','xLPS_{TNF}','nLPS_{TNF}','xLPS_{IL6}','nLPS_{IL6}','xLPS_{IL10}',...
%     'nLPS_{IL10}','xLPS_{IL1}','nLPS_{IL1}','xIL10_{TNF}','nIL10_{TNF}','xIL10_{IL6}','nIL10_{IL6}',...
%     'xIL10_{IL1}','nIL10_{IL1}','kHR_{BP}','kHR_{Temp}'})
% xtickangle(45)