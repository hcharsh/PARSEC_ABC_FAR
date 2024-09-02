% clear all;
% clc;
% close all;
% warning off;
% format shortG; % % housekeeping
% n_cluster = 5;
% specs_tim_d = 3;
% show_rows = 3;
% plot_cd = 0; plot_cd_rank = 0;
% 
function display_summary(n_cluster, specs_tim_d, show_rows, plot_cd, plot_cd_rank)
    res_folder = ['../RES_', num2str(specs_tim_d), 'unit/C', num2str(n_cluster)];
    eval_file = [res_folder, '/eval_file']; load([eval_file, '.mat']);
    
    %% show the good designs
    if show_rows > 0
        clear rank srt_Sep_expt_lbl srt_TimPts srt_estm_err;
        [vecs_Sort_TrainPARSEC, inds_Sort_TrainPARSEC] = sort(prm_err0(1:n_expt, 4));
        for ind = 1:min(show_rows, size(inds_Sort_TrainPARSEC, 1))
            ed_ind = inds_Sort_TrainPARSEC(ind, 1);
            rank(ind, 1) = ind;
            srt_Sep_expt_lbl{ind, 1} = Sep_tot_expt_lbl{ed_ind, 1};
            srt_TimPts{ind, 1}       = measure_tim_tot{ed_ind, 1};
            srt_estm_err(ind, 1)     = round(vecs_Sort_TrainPARSEC(ind, 1), 2);
            if test_set == 1
                srt_estm_err_test(ind, 1) = prm_err_full2d(ed_ind, 1);
                rank_test(ind, 1) = inds_Sort_TrainPARSEC(ed_ind, 1);
            end
    
        end
        st.Rank = rank; st.Expt_id = srt_Sep_expt_lbl; st.ToM = srt_TimPts; st.TrainEstmErro = srt_estm_err;
        % if test_set == 1
        %     st.Rank_test = rank_test; st.TestEstmError = srt_estm_err_test;
        % end
        Table_srt_train  = struct2table(st);
        disp(Table_srt_train);
    end
    
    %% plot the cummulative designs
    if plot_cd == 1
    
        x_rank = [0; SSort_Train.rank];
        x_err = [vecs_Sort_Train(1, 1); vecs_Sort_Train];
        y_PARSEC = [0; SSort_Train.good_PARSECs];
        y_RANDOM = [0; SSort_Train.good_RANDOMs];
    
        x = x_err;  x_lbl = 'Estm. error (arb. units)';
        if plot_cd_rank == 1
            x = x_rank; x_lbl = 'Rank';
        end
    
        f1 = figure('Units', 'normalized', 'Position', [0.35, 0.1, 0.3, 0.3]);
        ax1 = axes('Parent', f1); hold on;
        plot(x, y_PARSEC, '-k', 'Linewidth', 1.5);
        plot(x, y_RANDOM, '-r', 'Linewidth', 1.5);
        legend({'PARSEC', 'RANDOM'}, 'Location', 'best'); legend boxoff;
        xlabel(x_lbl); ylabel('Fraction of better designs');
        title(['Evaluation of ', num2str(n_cluster), '-sample designs']);
        set(ax1, 'FontSize', 12, 'Linewidth', 1);
    end
end
