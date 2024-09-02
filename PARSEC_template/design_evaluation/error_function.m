function WSSE = error_function(prm, data_file, prm_sample_ind)
    load(data_file)
    time_pts = Data_set(:, 1, prm_sample_ind);
    en = max(time_pts);
    tspan = [0:1:en]';
    dyn_data = Data_set(:, 2:size(Data_set, 2), prm_sample_ind);
    
    %% simulation
    % disp(10.^prm)
    [t,y] = ode45(@(t,y) model_system(t, y, prm), tspan, y0);
    dyn_pred = y(time_pts' + 1, specs_var_interest');
    %% error calculation
    error = sum(sum((dyn_data - dyn_pred).^2));
    WSSE = sqrt(error);
end