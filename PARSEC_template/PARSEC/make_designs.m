% clear all; clc; close all; warning off; format shortG; % % housekeeping
% n_expt = 100;
% n_cluster = 6;
% tim_d = 3;

function [D_file, intra_d_mean, intra_d_stdev] = make_designs(n_cluster, n_expt, full_Z, specs_tim_d, specs_var_interest)

    folder = ['Results/DES_', num2str(specs_tim_d), 'unit']; mkdir(folder);
    C_file = [folder, '/C',  num2str(n_cluster)];
    D_file = [folder, '/D',  num2str(n_cluster)];

    %% clustering
    if ~exist([C_file, '.mat'], 'file')
         clustering_sync(n_cluster, n_expt, full_Z, specs_tim_d, C_file);
    end
    load([C_file, '.mat'],...
        'k_cluster', 'MT', 'k2');
    intra_d         = k2.sum_intraDist; 
    intra_d_mean    = mean(intra_d);
    intra_d_stdev   = std(intra_d);
    %% selecting one candidates from each clusters
    for expt_ind = 1:n_expt
        clear bag;
        bag = k_cluster.sample(:, expt_ind);
        for cl_ind = 1:n_cluster % max(bag)
            clear int_vec size_vec;
            int_vec = find(bag==cl_ind);
            size_vec = size(int_vec, 1);
            clear rand_val rand_ind select_measure_ind;
            rand_val = rand;
            rand_ind = floor(rand_val * size_vec) + 1;
            select_measure_ind = int_vec(rand_ind, 1);
            PARSEC_ind(cl_ind, expt_ind) = select_measure_ind;
            PARSEC_CH(cl_ind, expt_ind) = bag(select_measure_ind, 1);
            TIM_PARSEC(cl_ind, expt_ind) = MT(select_measure_ind, 1);
        end
    end
    VAR_PARSEC = specs_var_interest;

    %% selecting multiple candidates from a few clusters
    for expt_ind = 1:n_expt
        clear bag size_cl sort_sc sort_ind;
        bag = k_cluster.sample(:, expt_ind);
        for cl_ind = 1:n_cluster % max(bag)
            clear int_vec;
            int_vec = find(bag==cl_ind);
            size_cl(cl_ind, 1) = size(int_vec, 1);
        end
        [sort_sc, sort_ind] = sort(size_cl);
        fill_m = 0; cl_ind = 1;
        while 1
            nm_m = sort_sc(n_cluster + 1 - cl_ind, 1);
            cl_m = sort_ind(n_cluster + 1 - cl_ind, 1);
            int_vec = find(bag==cl_m);
            if nm_m < (n_cluster - fill_m)
                TIM_ANTIPARSEC(fill_m + 1 : fill_m + nm_m, expt_ind) = MT(int_vec', 1);
                fill_m = fill_m + nm_m;
            else
                clear temp_perm;
                temp_perm = int_vec(randperm(nm_m), 1);
                TIM_ANTIPARSEC(fill_m + 1 : n_cluster, expt_ind) = ...
                    MT(temp_perm(1 : n_cluster - fill_m, 1)', 1);
                CLS_ANTIPARSEC(fill_m + 1 : n_cluster, expt_ind) = ...
                    cl_m;
                fill_m = n_cluster;
            end

            if (fill_m >= n_cluster)
                break;
            end
            if cl_ind >= n_cluster
                break;
            end
            cl_ind = cl_ind + 1;
        end
    end
    VAR_ANTIPARSEC = specs_var_interest;

    %% not the worst
    a = 0.3; % bigger this better it is
    for expt_ind = 1:n_expt
        clear bag size_cl sort_sc sort_ind;
        bag = k_cluster.sample(:, expt_ind);
        for cl_ind = 1:n_cluster % max(bag)
            clear int_vec;
            int_vec = find(bag==cl_ind);
            size_cl(cl_ind, 1) = size(int_vec, 1);
        end
        [sort_sc, sort_ind] = sort(size_cl);
        tot_nm = size(MT, 1);
        tW2_nm = tot_nm * a + n_cluster * (1 - a);
        small_set = []; cls_set = [];
        fill_m = 0; cl_ind = 1;
        while 1
            clear nm_m cl_m int_vec;
            nm_m = sort_sc(n_cluster + 1 - cl_ind, 1);
            cl_m = sort_ind(n_cluster + 1 - cl_ind, 1);
            int_vec = find(bag==cl_m);
            cls_vec = cl_m * ones(size(int_vec));
            fill_m = fill_m + nm_m;

            small_set   = [small_set; int_vec];
            cls_set     = [cls_set; int_vec];

            if (fill_m >= tW2_nm) && (cl_ind >= 2)
                break;
            end
            if cl_ind >= n_cluster
                break;
            end
            cl_ind = cl_ind + 1;
        end
        nmc = size(small_set, 1); temp_perm = randperm(nmc);
        reqd_perm = temp_perm(1, 1:n_cluster);
        temp_tim_ind = small_set(reqd_perm, 1)';

        TIM_ANTIPARSEC2(:, expt_ind) = MT(temp_tim_ind, 1);
        CLS_ANTIPARSEC2(:, expt_ind) = cls_set(reqd_perm, 1);
    end
    VAR_ANTIPARSEC2 = specs_var_interest;

    %% Random Synchronous
    nmc = size(MT, 1); % n_measurement_candidates
    for expt_ind = 1:n_expt
        clear temp_perm;
        temp_perm = randperm(nmc)';
        for cl_ind = 1:n_cluster % max(bag)
            clear temp_id var_id tim_id;
            tim_id = temp_perm(cl_ind, 1);
            TIM_RANDOM(cl_ind, expt_ind) = MT(tim_id, 1);
        end
    end
    VAR_RANDOM = specs_var_interest;

    clear rand_val rand_ind select_measure_ind...
        int_vec size_vec cl_ind expt_ind...
        temp_perm temp_id var_id tim_id;
    save(D_file);
end
