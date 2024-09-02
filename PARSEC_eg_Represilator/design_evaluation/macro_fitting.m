function macro_fitting(specs_tim_d, selected_clusters, n_expt)

addpath '../';

for n_cluster_ind = 1:size(selected_clusters, 1)
    n_cluster = selected_clusters(n_cluster_ind, 1); % disp(n_cluster)
    expt_folder = ['../data_set_', num2str(specs_tim_d), 'unit/C', num2str(n_cluster)];
    res_folder = ['../RES_', num2str(specs_tim_d), 'unit/C', num2str(n_cluster)];    
    lb = log10(0.3)*ones(6, 1); lb(1, 1) = -2; lb(2, 1) = -3;
    ub =  log10(30)*ones(6, 1); ub(1, 1) =  0; ub(2, 1) = -1;

    prm_test = []; load('../GT_val.mat', 'prm_test');
    if size(prm_test, 2) > 0
        test_set = 1;
    else
        test_set = 0;
    end

    for ED_ind = 1:n_expt
        fit_setup_file_UD(expt_folder, res_folder, 'train PARSEC_ED', ED_ind, lb, ub);
        fit_setup_file_UD(expt_folder, res_folder, 'train RANDOM_ED', 1 * ED_ind, lb, ub);
        if test_set == 1
            fit_setup_file_UD(expt_folder, res_folder, 'test PARSEC_ED', ED_ind, lb, ub);
            fit_setup_file_UD(expt_folder, res_folder, 'test RANDOM_ED', 1 * ED_ind, lb, ub);
        end
    end
    
    % evaluation_of_prm_estm(n_cluster, specs_tim_d, n_expt) % this is
    % already being called in the live script
    
end

end