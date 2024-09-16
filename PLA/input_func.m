function [tstim1, tstim2, dose, simtime] = input_func(mode)
%INPUT_FUNC Provides stimulation times, dose, and simulation time based on mode
%   Detailed explanation goes here

% Default simulation time
simtime = 8;

% Define parameters for each mode using cell arrays
params = {...
    {1, struct('tstim1', 0.1, 'tstim2', 0.1, 'dose', 0.5, 'simtime', 24)}, ...
    {2, struct('tstim1', 0.1, 'tstim2', 0.1, 'dose', 1, 'simtime', 24)}, ...
    {3, struct('tstim1', 0.1, 'tstim2', 0.1, 'dose', 2, 'simtime', 24)}, ...
    {4, struct('tstim1', 0.1, 'tstim2', 0.1, 'dose', 0.3, 'simtime', 24)}, ...
    {5, struct('tstim1', 0.1, 'tstim2', 0.1, 'dose', 1)}, ...
    {6, struct('tstim1', 0.1, 'tstim2', 0.1, 'dose', 2)}, ...
    {7, struct('tstim1', 0.1, 'tstim2', 3, 'dose', 1)} ...
};

% Get parameters based on mode, use default values if mode is not specified
found = false;
for i = 1:numel(params)
    if params{i}{1} == mode
        mode_params = params{i}{2};
        tstim1 = mode_params.tstim1;
        tstim2 = mode_params.tstim2;
        dose = mode_params.dose;
        if isfield(mode_params, 'simtime')
            simtime = mode_params.simtime;
        end
        found = true;
        break;
    end
end

if ~found
    tstim1 = 0;
    tstim2 = 0;
    dose = 0;
end

end
