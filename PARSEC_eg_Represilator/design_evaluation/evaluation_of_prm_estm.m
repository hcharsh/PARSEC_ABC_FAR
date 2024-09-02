% clear all;
% clc;
% close all;
% warning off;
% format shortG; % % housekeeping
% n_cluster = 7;
% specs_tim_d = 2;

function ar = evaluation_of_prm_estm(n_cluster, specs_tim_d, n_expt)

    res_folder = ['../RES_', num2str(specs_tim_d), 'unit/C', num2str(n_cluster)];
    eval_file = [res_folder, '/eval_file'];
    
    prm_test = [];
    load('../GT_val.mat');
    if size(prm_test, 2) > 0
        test_set = 1;
    else test_set = 0;
    end
    
    %%
    
    expt_lbl{1, 1} = 'GT';
    prefix_list = {'train PARSEC_ED', 'train RANDOM_ED', 'test PARSEC_ED', 'test RANDOM_ED'};
    n_expt_vec = n_expt * [1, 1, 1, 1]';
    prm_sample_vec = [size(prm, 2), size(prm, 2), size(prm_test, 2), size(prm_test, 2)]';
    
    total_ed_ind = 0; index = 0;
    
    for prefix_ind = 1:(2 + 2*test_set)
        clear prefix; pf = prefix_list{1, prefix_ind}; 
%         clear n_expt; n_expt = n_expt_vec(prefix_ind, 1);
        for ed_index = 1:n_expt_vec(prefix_ind, 1)
            total_ed_ind = total_ed_ind + 1;
            prm_err_full2d(total_ed_ind, 1) = 0;
            ed_pf_mrkr(total_ed_ind, 1:(2 + 2*test_set)) = zeros(1, (2 + 2*test_set)); ed_pf_mrkr(total_ed_ind, prefix_ind)  = 1;
            Sep_tot_expt_lbl{total_ed_ind, 1}    = [pf, ' #', num2str(ed_index)];
            for prm_sample_ind = 1:prm_sample_vec(prefix_ind, 1)
                if prefix_ind <= 2
                    prm_log_GT = round(prm(:, prm_sample_ind), 2);
                else
                    prm_log_GT = round(prm_test(:, prm_sample_ind), 2);
                end
                clear time_pts FAR_file cmf_est PRM_nest;
                FAR_file = [res_folder, '/full_', pf, num2str(ed_index), '_s', num2str(prm_sample_ind)];
                index = index + 1;
                pf_mrkr(index, 1:(2 + 2*test_set)) = zeros(1, (2 + 2*test_set)); pf_mrkr(index, prefix_ind)  = 1;
                load(FAR_file, 'info', 'cmf_est', 'PRM_nest');
                load(info.data_file, 'Data_set'); time_pts = Data_set(:, 1, prm_sample_ind);
                cmf_est = round(cmf_est, 2);
                sample_size(index, 1) = (size(Data_set, 2) - 1) * size(time_pts, 1);
                clear info Data_set;
                GT{index, :}          = mat2str(prm_log_GT);

                Sep_expt_lbl{index, 1}    = [pf, ' #', num2str(ed_index), 's #', num2str(prm_sample_ind)];
                expt_lbl{index + 1, 1}    = [pf, ' #', num2str(ed_index), 's #', num2str(prm_sample_ind)];
                measure_tim{index, 1} = mat2str(time_pts);
                blank{index, 1} = '';

                n_nest = size(PRM_nest, 3); n_prm = size(PRM_nest, 2);
                lst = (n_nest - 1) * n_prm + [1 : n_prm];

                A_25(index, :)        = cmf_est(lst, 2)';
                A_median(index, :)    = cmf_est(lst, 3)';
                A_75(index, :)        = cmf_est(lst, 4)';
                A_diff(index, :)      = cmf_est(lst, 4)' - cmf_est(lst, 2)';
                clear cmf_est;
                clear dyn_err prm_err sort_dyn_err sort_prm_err d50 p50 dyn_mean prm_mean;
                prm_err = fn_dev_prm(prm_log_GT, PRM_nest);
                sort_prm_err = sort(prm_err); co = size(sort_prm_err, 1);

                p50 = mean(sort_prm_err([co/2, 1 + co/2], 1));
                prm_mean_err = mean(sort_prm_err);
                prm_err0(index, :) = [sort_prm_err(co/4, 1), p50, ...
                    sort_prm_err(3*co/4, 1), prm_mean_err];
                
                prm_err_full(prefix_ind, ed_index, prm_sample_ind) = prm_mean_err;
                prm_err_full2d(total_ed_ind, 1) = prm_err_full2d(total_ed_ind, 1) + prm_mean_err/prm_sample_vec(prefix_ind, 1);
                measure_tim_tot{total_ed_ind, 1} = mat2str(time_pts);
            end
        end
    end
    PARSEC = mean(prm_err0(1:n_expt, 4));
    RANDOM = mean(prm_err0(n_expt + 1:n_expt + n_expt, 4));
    ar = RANDOM/PARSEC;
    disp(['For ', num2str(n_cluster), ' clusters:'])
    disp(['Mean estm. error for ', num2str(n_expt), ' PARSEC design(s) is ', num2str(PARSEC), ' arbitary units'])
    disp(['Mean estm. error for ', num2str(n_expt), ' random design(s) is ', num2str(RANDOM), ' arbitary units'])
    disp(['The accuracy ratio is ', num2str(RANDOM/PARSEC)])
    for prm_ind = 1:n_prm
        A_sep(:, prm_ind) = abs(A_median(:, prm_ind) - prm_log_GT(prm_ind, 1));
    end

    A_Sep_lbl = Sep_expt_lbl;

    %% Time of measurements for all the designs
    S_all_dsn_info.expt_ind     = Sep_expt_lbl;
    S_all_dsn_info.GT           = GT;
    S_all_dsn_info.TimPts       = measure_tim;
    S_all_dsn_info.sample_size  = sample_size;
    S_all_dsn_info.estm_err     = round(prm_err0(:, 4), 2);

    T_all_dsn_info              = struct2table(S_all_dsn_info);

    %% summary of error in parameter estimation

    Errs.expt_ind      = A_Sep_lbl;
    Errs.errPRM        = blank;

    Errs.estm_err_q25   = round(prm_err0(:, 1), 2);
    Errs.estm_err_q50   = round(prm_err0(:, 2), 2);
    Errs.estm_err_q75   = round(prm_err0(:, 3), 2);
    Errs.errPRM         = blank;
    Errs.estm_err_avg   = round(prm_err0(:, 4), 2);

    T_Errs = struct2table(Errs);


    %% Train
    Train_ind = n_expt_vec(1, 1) + n_expt_vec(2, 1);
    Train.expt_ind    = Sep_tot_expt_lbl(1:Train_ind, 1);
    Train.estm_error  = prm_err_full2d(1:Train_ind, 1);
    TTrain = struct2table(Train);

    [vecs_Sort_Train, inds_Sort_Train]      = sort(prm_err_full2d(1:Train_ind, 1));
    Sort_Train                              = ed_pf_mrkr(inds_Sort_Train, :);
    cSort_Train(1, :)                       = Sort_Train(1, :);
    for nr = 2:size(Sort_Train, 1)
        cSort_Train(nr, :)   = cSort_Train(nr - 1, :) + Sort_Train(nr, :);
    end
    
    SSort_Train.rank             = [1:Train_ind]';
    SSort_Train.sort_expt_ind    = Sep_tot_expt_lbl(inds_Sort_Train, 1);
    SSort_Train.sort_estm_error  = vecs_Sort_Train;

    SSort_Train.good_PARSECs        = cSort_Train(:, 1)/n_expt_vec(1, 1);
    SSort_Train.good_RANDOMs        = cSort_Train(:, 2)/n_expt_vec(2, 1);

    TSort_Train = struct2table(SSort_Train);

    %%
    excel_file = [eval_file, '.xlsx'];
    warning('off','MATLAB:xlswrite:AddSheet'); %optional
    writetable(T_all_dsn_info, excel_file, 'Sheet', 1); %, 'All Designs');
    writetable(T_Errs, excel_file, 'Sheet', 2); %, 'Estimate error');
    writetable(TTrain, excel_file, 'Sheet', 3); %, 'Error in Estimate');
    writetable(TSort_Train, excel_file, 'Sheet', 4); %, 'Error in Estimate');




    if test_set == 1
        %% Test
        Test_ind = sum(n_expt_vec);
        Test.expt_ind    = Sep_tot_expt_lbl(Train_ind + 1:Test_ind, 1);
        Test.estm_error  = prm_err_full2d(Train_ind + 1:Test_ind, 1);
        TTest = struct2table(Test);

        [vecs_Sort_Test, inds_Sort_Test]    = sort(prm_err_full2d(Train_ind + 1:Test_ind, 1));
        Sort_Test                           = ed_pf_mrkr(inds_Sort_Test, :);
        cSort_Test(1, :)                    = Sort_Test(1, :);
        for nr = 2:size(Sort_Test, 1)
            cSort_Test(nr, :)   = cSort_Test(nr - 1, :) + Sort_Test(nr, :);
        end

        SSort_Test.rank             = [1:Test_ind - Train_ind]';
        SSort_Test.sort_expt_ind    = Sep_tot_expt_lbl(inds_Sort_Test + Train_ind, 1);
        SSort_Test.sort_estm_error  = vecs_Sort_Test;

        SSort_Test.good_PARSECs        = cSort_Test(:, 1)/n_expt_vec(1, 1);
        SSort_Test.good_RANDOMs        = cSort_Test(:, 2)/n_expt_vec(2, 1);

        TSort_Test = struct2table(SSort_Test);
        writetable(TTest, excel_file, 'Sheet', 5); %, 'Error in Estimate');
        writetable(TSort_Test, excel_file, 'Sheet', 6); %, 'Error in Estimate');
    end

    save([eval_file, '.mat']);

end
