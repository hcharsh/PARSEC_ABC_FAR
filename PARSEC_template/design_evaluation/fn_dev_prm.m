function prm_err = fn_dev_prm(prm_log_GT, PRM_nest)

    cut_off     = size(PRM_nest, 1);
    n_prm       = size(PRM_nest, 2);
    n_nest      = size(PRM_nest, 3);
    PARAMETER = PRM_nest(:, :, n_nest);
    clear PRM_nest;

    free_log_GT = prm_log_GT(1:n_prm, 1);
    for prm_ind = 1:cut_off
        clear prm dev sq_dev rssq;
        prm = PARAMETER(prm_ind, :)';
        dev = abs(prm - free_log_GT);
        sq_dev = dev.^2; ssq = sum(sq_dev);
        rssq = sqrt(ssq);
        prm_err(prm_ind, 1) = rssq;
    end

end