function gen_data(n_cluster, specs_tim_d, D_file, show_designs)
    
    %% the folder where the designs will be saved
    expt_folder = ['data_set_', num2str(specs_tim_d), 'unit/C', num2str(n_cluster)];
    mkdir(expt_folder);

    %% load D_file (the design file) and parameter/initial condition information
    load([D_file, '.mat'], 'TIM_PARSEC', 'TIM_ANTIPARSEC', 'TIM_ANTIPARSEC2', 'TIM_RANDOM', 'specs_var_interest');
    prm_test = [];
    load('GT_val.mat', 'prm', 'tspan', 'y0', 'prm_test');
    if size(prm_test, 2) > 0
        test_set = 1;
    else test_set = 0;
    end
    %% feasible time points of observtions
    st = 0; en = tspan(1, 2); clear tspan; tspan = [0:1:en];

    %% 
    for train_test = 1:(test_set + 1)
        clear prm_temp pre_pf;
        if train_test == 1
            prm_temp = prm; pre_pf = 'train';
        else
            prm_temp = prm; pre_pf = 'test';
        end

        for prm_sample = 1:size(prm_temp, 2)
            clear t y;
            [t, y] = ode45(@(t, y) model_system(t, y, prm_temp(:, prm_sample)), tspan, y0);
        %     y_aug = abs(y + 0.05*rand(size(y))); clear y; y = y_aug; clear y_aug;
            y_3d(:, :, prm_sample) = [t, y(:, specs_var_interest')];
        end

        %% Designs
        for prefix_ind = 1:2
            clear pf TIM_MAT;
            if prefix_ind == 1
                TM_MAT = TIM_PARSEC; pf = [pre_pf, ' PARSEC'];
            elseif prefix_ind == 2
                TM_MAT = TIM_RANDOM; pf = [pre_pf, ' RANDOM'];
            elseif prefix_ind == 3
                TM_MAT = TIM_ANTIPARSEC; pf = [pre_pf, ' ANTIPARSEC'];
            elseif prefix_ind == 4
                TM_MAT = TIM_ANTIPARSEC2; pf = [pre_pf, ' ANTIPARSEC2'];
            end

            n_ED = size(TM_MAT, 2);
            for ED_ind = 1:n_ED
                clear Data_set t_fin t_indices;
                t_fin = sort(TM_MAT(:, ED_ind));
                t_indices = t_fin' + 1; % disp(t_indices - 1);
                Data_set = y_3d(t_indices, :, :);
                save([expt_folder, '/', pf, '_ED', num2str(ED_ind), '.mat'], ...
                    'y0', 'tspan', 'Data_set', 'specs_var_interest');
                
                if rand < show_designs
                    figure; hold on;
                    title([num2str(n_cluster), ' samples ', pf, ' #', num2str(ED_ind)])
                    clear t y;
                    for ind = 1:( size(Data_set, 2) - 1 )
                        plot(Data_set(:, 1, 1), Data_set(:, ind + 1, 1), '.', 'Markersize', 25);
                        legend_text{1, ind} = ['Data: V ', num2str(specs_var_interest(ind, 1))];
                    end
                    for ind = 1:( size(Data_set, 2) - 1 )
                        plot(tspan, y_3d(tspan + 1, ind + 1, 1), '--');
                        legend_text{1, ind + size(Data_set, 2) - 1} = ['Actual: V ', num2str(specs_var_interest(ind, 1))];
                    end
                    legend(legend_text); legend boxoff;
                    xlabel('Time'); ylabel('Molecules');
                    hold off;
                end

            end

        end

        %% Designs
        clear Data_set; Data_set = y_3d;
        save([expt_folder, '/data0.mat'], ...
            'y0', 'tspan', 'Data_set', 'specs_var_interest');
    end
end