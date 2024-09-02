function err_vec = fn_calc_error_UD(parameter_matrix, info)

    n_rlz           = size(parameter_matrix, 1);
    err_vec         = zeros(n_rlz, 1);
    % disp('fn_calc_error_UD')
    parfor iter = 1:n_rlz
        % simulation
        prm_arr = parameter_matrix(iter, :);
        prm     = prm_arr';        
        if info.work_progress_nest == 1
            if mod(iter-1, round(n_rlz/3))==0
                disp([num2str((iter-1)*100/n_rlz), '% done for the nest']);
            end
        end
        WSSE = error_function(prm, info.data_file, info.prm_sample_ind);
        err_vec(iter, 1) = WSSE;
    end    
end