% clear all; clc; close all; warning off; format shortG; % % housekeeping
% n_cluster = 5; n_iter = 5; % n_samples = 2;
% specs_tim_d = 2;
% load('TSA_full_2unit.mat', 'full_Z');

function clustering_sync(n_cluster, n_iter, full_Z, specs_tim_d, C_file)
    TSA_file = ['Results/TSA_prm1_',  num2str(specs_tim_d), 'unit.mat'];
    load(TSA_file, 'MT');
    
    z_prm_mat = zscore(full_Z);
    k1.time_point = MT;
    k_cluster.time_point = MT;
    %% k_means

    % disp(z_prm_mat); disp(n_cluster)
    iter = 0;
    while iter < n_iter
        clear idx centroids sumd DIST;
        [idx, centroids, sumd, DIST] = kmeans(z_prm_mat, n_cluster);
        iter = iter + 1;
        k_cluster.iteration{iter, 1} = ['Iteration #', num2str(iter)];
        k2.lbl{iter, 1} = ['Iteration #', num2str(iter)];
        k_cluster.sample(:, iter) = idx;
        for ind = 1:size(idx, 1)
            k_cluster.cluster_id{ind, iter} = ['in cluster ', num2str(idx(ind, 1))];
        end
        k2.intraDist_cluster(iter, :) = round(sumd, 4);
        k2.sum_intraDist(iter, :) = round(sum(sumd), 4);
        k2.totalDist(iter, :) = sum(sum(DIST));
        % disp(k1); disp(k_cluster);
    end
    k1.iteration = k_cluster.cluster_id;
    % disp(k1); disp(k2);
    k1_Table  = struct2table(k1);
    k2_Table  = struct2table(k2);
    xl_file = [C_file, '.xlsx'];
    warning('off','MATLAB:xlswrite:AddSheet'); %optional
    writetable(k1_Table, xl_file, 'Sheet', 1);
    writetable(k2_Table, xl_file, 'Sheet', 2);
    clear k1_Table k2_Table;
    save([C_file, '.mat'])

end