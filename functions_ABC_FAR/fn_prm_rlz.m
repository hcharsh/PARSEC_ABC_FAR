function PARAMETER = fn_prm_rlz(currDist, lb, ub, n_rlz)
    [n_bins, n_prm] = size(currDist);
    unifDist = linspace(0, 1, n_bins)';

    X1 = lhsdesign(n_rlz, n_prm,'iterations', 5); % uniform rand
    for ind_prm = 1:n_prm
        low = lb(ind_prm, 1); high = ub(ind_prm, 1); dif = high - low;
        for ind_rlz = 1:n_rlz

            rlz_unif = X1(ind_rlz, ind_prm);
            bin_dist = 1;
            while 1
                if rlz_unif<currDist(bin_dist, ind_prm)
                    break;
                else
                    bin_dist = bin_dist + 1;
                end
            end
%             disp([rlz_unif, bin_dist]);
            % % unif -> dist
            d1 = currDist(bin_dist, ind_prm);
            d0 = currDist(bin_dist - 1, ind_prm);
            x1 = unifDist(bin_dist, 1);
            x0 = unifDist(bin_dist - 1, 1);
            rlz_dist = x0 + (x1 - x0)*(rlz_unif - d0)/(d1 - d0);
            PARAMETER_rel(ind_rlz, ind_prm) = rlz_dist;
            val = low + rlz_dist*dif;
            PARAMETER(ind_rlz, ind_prm) = val;
            clear d0 d1 x0 x1 rlz_unif bin_dist rlz_dist val;
        end
        clear low high dif;
    end
    % % test
    pmf = zeros(n_bins-1, n_prm);
%     'same'
%     size(pmf);
    for ind_prm = 1:n_prm
        for ind_rlz = 1:n_rlz
            clear rlz bin_num
            rlz = PARAMETER_rel(ind_rlz, ind_prm);
            bin_num = floor(rlz*(n_bins - 1) + 1);
            pmf(bin_num, ind_prm) = pmf(bin_num, ind_prm) + 1/n_rlz;
        end
    end
end