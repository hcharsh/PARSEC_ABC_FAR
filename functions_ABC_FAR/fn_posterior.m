function [PostDist] = fn_posterior(sort_PRM, lb, ub, n_bins)
    thres = size(sort_PRM, 1);
    num_prm = size(sort_PRM, 2);
    PostDist = zeros(n_bins, num_prm);
    
    for ind_iter = 1:thres
        clear ind_interest;
        for ind_prm = 1:num_prm
            clear abs_rlz rlz;
            abs_rlz = sort_PRM(ind_iter, ind_prm);
            rlz = (abs_rlz - lb(ind_prm, 1))/(ub(ind_prm, 1) - lb(ind_prm, 1));
            bin_rlz = floor(n_bins*rlz) + 1;
            PostDist(bin_rlz, ind_prm) = PostDist(bin_rlz, ind_prm) + (1/thres);
        end
    end
end

