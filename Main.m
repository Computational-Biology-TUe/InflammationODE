clear all
close all
warning('off')

%% Load data
% Load the data from the 'Data.mat' file
data = load('Data.mat');

%% Prepare optimization
% Initial parameter array
param0 = ones(6,1);

% Set options for the optimization algorithm
lsq_options = optimset('Display','iter', 'Algorithm','trust-region-reflective', 'TolFun',2e-10, 'MaxFunEvals',2e4);

% Select the dataset to fit
% [1 2 3] --> Dillingh 2014 (0.5 - 1 - 2 ng/kg bolus)
% [4]    --> Taudorf   2007 (0.3 ng/kg bolus)
% [5 6 7] --> Kiers    2017 (1 - 2 ng/kg bolus - continuous 3h)

mode = 5;     % 1 ng/kg bolus Kiers
%mode = 6;    % 2 ng/kg bolus Kiers

% Lower and upper bounds for the parameters
lb = [0 0 0 0 1.25 0];
ub = [1e6 1e6 1e6 2 2 2];

%% Estimate parameters
% Perform parameter estimation using least squares optimization
[p_opt,~,residual,~,~,~,J] = lsqnonlin(@cost_func, param0, lb, ub, lsq_options, mode, data);

%% Plot results
% Choose validation dataset for plotting

if mode == 5
    plot_mode = 5;            % 1 ng/kg bolus
    %plot_mode = 7;           % Continuous Infusion
elseif mode == 6
    plot_mode = [1 2 3 4];     % Different Dosages
end

figure
for i = plot_mode
    % Get the model data
    [t, TNF, IL6, IL10, IL1, t2, Temp, BP, HR, color, corrTNF, corrIL6, corrIL10] = data_model(i, data);

    % Simulate the system with the optimized parameters
    [X, Y] = model_code(i, p_opt, 0, 0, 0, Temp(1), BP(1), HR(1));

    % Define the cytokine subplots and their properties
    subplot_info = {
        struct('index', 1, 'data', TNF, 'model_data', Y(:,5), 'title', 'TNF', 'ylabel', 'Concentration [pg/ml]', 'ylim', [0.1 1000]),
        struct('index', 2, 'data', IL6, 'model_data', Y(:,7), 'title', 'IL6', 'ylabel', 'Concentration [pg/ml]', 'ylim', [0.1 2000]),
        struct('index', 3, 'data', IL10, 'model_data', Y(:,9), 'title', 'IL10', 'ylabel', 'Concentration [pg/ml]', 'ylim', [0.1 1000], 'condition', corrIL10),
        struct('index', 4, 'data', IL1, 'model_data', Y(:,11), 'title', 'IL1-b', 'ylabel', 'Concentration [pg/ml]', 'ylim', [0.01 1000], 'condition', ismember(i, [5 6 7]))
        };

    % Loop through the cytokines
    for k = 1:length(subplot_info)
        info = subplot_info{k};
        subplot(2, 3, info.index)
        plot(X, info.model_data, 'Color', [color 1], 'LineWidth', 3)
        hold on
        if ~isfield(info, 'condition') || info.condition
            plot(t, info.data, 'kx', 'LineWidth', 2, 'Color', color, 'MarkerSize', 7)
        end
        grid on
        title(info.title)
        ylabel(info.ylabel)
        xlabel('Time [h]')
        set(gca, 'YScale', 'log')
        ylim(info.ylim)
    end

    % Plot temperature change
    subplot(2, 3, 5)
    plot(X, Y(:,12)-Y(1,12), 'Color', [color 1], 'LineWidth', 3)
    hold on
    plot(t2, Temp-Temp(1), 'kx', 'LineWidth', 2, 'Color', color, 'MarkerSize', 7)
    grid on
    title('Temperature change')
    set(gca, 'YScale', 'log')
    ylim([0.1 10])
    ylabel('degrees [C]')
    xlabel('Time [h]')

    % Plot heart rate change
    subplot(2, 3, 6)
    plot(X, Y(:,14)-Y(1,14), 'Color', [color 1], 'LineWidth', 3)
    hold on
    plot(t2, HR-HR(1), 'kx', 'LineWidth', 2, 'Color', color, 'MarkerSize', 7)
    grid on
    title('Heart rate change')
    ylabel('Rate [bpm]')
    xlabel('Time [h]')
    set(gca, 'YScale', 'log')
    ylim([1 100])
end
