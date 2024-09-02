function cmf_est = fn_estimate_cmf(info, nest_ind)

    %% load the relevant ABC_FAR results
    load([info.info_fit_folder, '/', info.ABC_FAR_file, '.mat'], 'Pmf_all');
    lb = info.lb;
    ub = info.ub;

    %% distribution info
    n_bin = size(Pmf_all, 1);
    n_prm = size(Pmf_all, 2);
    % nest_ind = size(Pmf_all, 3);
    PMF = Pmf_all(:, :, nest_ind);

    %% simulation

    CMF = zeros(n_bin, n_prm);
    CMF(1, :) = PMF(1, :);
    for bin_iter = 2:n_bin
        CMF(bin_iter, :) = PMF(bin_iter, :) + CMF(bin_iter - 1, :);
    end
    
    DELTA_PRM = ub - lb;
    for ind_prm = 1:n_prm
        for bin_iter = 1:n_bin
            PRM_VEC(bin_iter, ind_prm) = lb(ind_prm, 1) + (bin_iter - 0.5)*DELTA_PRM(ind_prm, 1)/(n_bin - 0.5);
        end
    end

    %% find median, 0.25, 0.75, 0.1 and 0.9
    q_vec = [0.1; 0.25; 0.5; 0.75; 0.9];
    cmf_est = zeros(n_prm, size(q_vec, 1));

    for ind_prm = 1:n_prm
        for ind = 1:size(q_vec, 1)
            f_val = q_vec(ind, 1);
            bin_iter = 1;
            while bin_iter < 1+n_bin
                if CMF(bin_iter, ind_prm) > f_val
                    break;
                end
                bin_iter = bin_iter + 1;
            end
            if bin_iter == 1
                prev_val = 0;
            else
                prev_val = CMF(bin_iter - 1, ind_prm);
            end
            slop = CMF(bin_iter, ind_prm) - prev_val;
            bin_num = (f_val - prev_val)/slop + bin_iter - 1;
            est = lb(ind_prm, 1) + (bin_num)*DELTA_PRM(ind_prm, 1)/(n_bin);
            cmf_est(ind_prm, ind) = est;
        end
    end    
end
