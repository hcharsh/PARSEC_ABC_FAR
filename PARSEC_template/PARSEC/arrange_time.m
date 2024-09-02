function [mat_arrange] = arrange_time(mat_var, mat_tim, var_int)
n_expt = size(mat_var, 2);
for expt_ind = 1:n_expt
    clear temp_var temp_tim;
    temp_var = mat_var(:, expt_ind);
    temp_tim = mat_tim(:, expt_ind);
    for var_ind = 1:size(var_int, 2)
        clear var_indices;
        var_indices = find(temp_var == var_int(1, var_ind));
        temp_indices = temp_tim(var_indices', 1);
%         disp([temp_var(var_indices', 1), temp_indices]);
        mat_arrange{expt_ind, var_ind} = sort(temp_indices);
    end
end