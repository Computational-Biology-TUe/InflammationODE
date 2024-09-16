function e = PLA_cost_func(param)
    data = load('Data.mat');
    mode = 5;

    e = [];
    for i = 1:length(mode)
        x = mode(i);

        % Retrieve experimental data
        [t, TNF, IL6, IL10, IL1, ~, Temp, BP, HR, ~, corrTNF, corrIL6, corrIL10] = data_model(x, data);

        % Simulate cytokine concentrations
        [tsim, ysim] = model_code(x, param, 0, 0, 0, Temp(1), BP(1), HR(1));

        % Interpolate simulated cytokine concentrations to match experimental time points
        TNF_simulated = interp1(tsim, real(ysim(:,5)), t, 'spline');
        IL6_simulated = interp1(tsim, real(ysim(:,7)), t, 'spline');

        % Calculate squared error for TNF and IL6
        TNFerror = corrTNF * ((TNF - TNF_simulated) ./ max(TNF)).^2;
        IL6error = corrIL6 * ((IL6(1:7) - IL6_simulated(1:7)) ./ max(IL6)).^2;

        % Calculate error for IL10 if available
        if corrIL10 == 1
            IL10_simulated = interp1(tsim, real(ysim(:,9)), t, 'spline');
            IL10error = corrIL10 * ((IL10 - IL10_simulated) ./ max(IL10)).^2;
            e = [e; TNFerror; IL6error; IL10error];
        else
            e = [e; TNFerror; IL6error];
        end
    end


end
