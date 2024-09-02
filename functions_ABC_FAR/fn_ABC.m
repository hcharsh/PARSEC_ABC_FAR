function fn_ABC(info)

ABC_start_time = datetime;
info.ABC_start_time = ABC_start_time;

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% few checks regarding the fitting options
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% check for lower and upper bounds
if ~isfield(info, 'lb')
    error('Please provive lower bound in field, lb, of the input argument of this function')
end

if ~isfield(info, 'ub')
    error('Please provive upper bound in field, ub, of the input argument of this function')
end

if size(info.lb, 1) ~= size(info.ub, 1)
    warning(['Lower bound and upper bound vectors have unequal sizes. ',...
        'The number of parameters will correspond to the ', ...
        'smaller of the vectors.']);
    temp_lb = info.lb;
    temp_ub = info.ub;
    n_prm = min([size(temp_lb, 1), size(temp_ub, 1)]);
    info.lb = temp_lb(1:n_prm, 1);
    info.ub = temp_ub(1:n_prm, 1);
end
lb = info.lb;
ub = info.ub;

bound_check = (lb>=ub);
if sum(bound_check) ~= 0
    error(['Atleast one of the lower bound is higher than or equal', ...
        'to the corresponding upper bound.', ...
        'Please ensure upper bounds are strictly higher than the',...
        'correponding lower bounds.']);
end
clear bound_check;

%% fitting options
if ~isfield(info, 'far')
    info.far = 0.025;
    info.default.far = 1;
end
acc_rate = info.far;

if ~isfield(info, 'cut_off')
    info.cut_off = 250;
    info.default.cut_off = 1;
end
cut_off = info.cut_off;

if ~isfield(info, 'n_nest')
    info.n_nest = 8;
    info.default.n_nest = 1;
end
n_nest = info.n_nest;

if ~isfield(info, 'n_bins')
    if isfield(info, 'PriorDist')
        nbins = size(info.PriorDist, 1);
        info.n_bins = nbins - 1;
        warning('n_bins is assumed to correspond to the size of PriorDist supplied.')
        clear nbins;
    else
        info.n_bins = 30;
        info.default.n_bins = 1;
    end
end
n_bins = info.n_bins;

if ~isfield(info, 'alpha')
    info.alpha = 0.1;
    info.default.alpha = 1;
end
alpha = info.alpha;

if ~isfield(info, 'use_prev')
    info.use_prev = 1;
    info.default.use_prev = 1;
end
use_prev = info.use_prev;

if ~isfield(info, 'PriorDist')
    n_prm = min([size(lb, 1), size(ub, 1)]);
    for ind_prm = 1:n_prm
        unif_Dist(:, ind_prm) = linspace(0, 1, n_bins + 1)'; % uniform distribution
    end
    info.PriorDist = unif_Dist;
    clear n_prm unif_Dist;
    info.default.PriorDist = 1;
end
PriorDist = info.PriorDist;

if ~isfield(info, 'sampling_noise_Dist')
    n_prm = min([size(lb, 1), size(ub, 1)]);
    for ind_prm = 1:n_prm
        unif_Dist(:, ind_prm) = linspace(0, 1, n_bins + 1)'; % uniform distribution
    end
    info.sampling_noise_Dist = unif_Dist;
    clear n_prm unif_Dist;
    info.default.sampling_noise_Dist = 1;
end
sampling_noise_Dist = info.sampling_noise_Dist;

PriorDist_bins = size(PriorDist, 1) - 1;
noise_bins = size(sampling_noise_Dist, 1) - 1;
min_nbins = min([n_bins, noise_bins, PriorDist_bins]);
if n_bins ~= min_nbins
    clear n_bins;
    info.n_bins = min_nbins;
    n_bins = info.n_bins;
    warning('Bining of PriorDist and/or sampling_noise_Dist is inconsistent with n_bins. n_bins is updated');
end

if PriorDist_bins ~= min_nbins
    current_res = linspace(0, 1, PriorDist_bins + 1)';
    reduced_res = linspace(0, 1, min_nbins + 1)';
    new_PriorDist = interp1(current_res, PriorDist, reduced_res);
    clear PriorDist current_res reduced_res;
    info.PriorDist = new_PriorDist;
    clear new_PriorDist;
    PriorDist = info.PriorDist;
    warning('Bining inconsistency. Resolution of PriorDist is reduced accordingly.');
end

if noise_bins ~= min_nbins
    current_res = linspace(0, 1, noise_bins + 1)';
    reduced_res = linspace(0, 1, min_nbins + 1)';
    new_sampling_noise_Dist = interp1(current_res, sampling_noise_Dist, reduced_res);
    clear sampling_noise_Dist current_res reduced_res;
    info.sampling_noise_Dist = new_sampling_noise_Dist;
    clear new_sampling_noise_Dist;
    sampling_noise_Dist = info.sampling_noise_Dist;
    warning('Bining inconsistency. Resolution of sampling_noise_Dist is reduced accordingly.');
end

clear PriorDist_bins noise_bins min_nbins;

%% folder and file names
if ~isfield(info, 'info_fit_folder')
    info.info_fit_folder = 'fit_folder';
    info.default.info_fit_folder = 1;
    warning(['Results will be saved in the folder: ', info.info_fit_folder])
end
mkdir(info.info_fit_folder);
Filename = sprintf('ABC_%s', datestr(now,'mm-dd-yyyy HH-MM'));
if ~isfield(info, 'ABC_FAR_file')
    info.ABC_FAR_file = Filename;
    info.default.ABC_FAR_file = 1;
    warning(['Results will be saved in the file: ', info.ABC_FAR_file, '.mat'])
end
clear Filename;
Filename = info.ABC_FAR_file;
if ~isfield(info, 'ABC_FAR_results_file')
    info.ABC_FAR_results_file = ['p_', Filename];
    info.default.ABC_FAR_results_file = 1;
    warning(['Selected statistics of results will be saved in the file: ', info.ABC_FAR_results_file, '.mat'])
end

if isfield(info, 'default')
    warning('Default values for certain fitting options/information were used.')
end

%% additional details
if ~isfield(info, 'work_progress_time')
    info.work_progress_time = 0;
    info.default.work_progress_time = 1;
end
if ~isfield(info, 'work_progress_nest')
    info.work_progress_nest = 1;
    info.default.work_progress_nest = 1;
end
if ~isfield(info, 'diary_reqd')
    info.diary_reqd = 1;
    info.default.diary_reqd = 1;
end
diary_reqd = info.diary_reqd;
diary_file = ['diary_', Filename];

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% the algorithm starts here 
%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

n_prm = min([size(lb, 1), size(ub, 1)]);
for ind_bin = 1:n_bins
    Pmf_all(ind_bin, :, 1) = PriorDist(ind_bin + 1, :) - PriorDist(ind_bin, :);
end

% rv_vec = zeros(n_bins, n_prm);
% for ind_prm = 1:n_prm
%     clear arr;
%     arr = linspace(lb(ind_prm, 1), ub(ind_prm, 1), n_bins)';
%     rv_vec(:, ind_prm) = arr;
% end

%% distribution biased by error/deviation from observations
PostDist = PriorDist; % for the first iteration
cmf_est = [];
alpha_CURRENT(1, :) = [0, 0, 0];

for nest_ind = 1:n_nest
    tic;
    info.current_nest = nest_ind;
    if info.work_progress_nest == 1
        disp(['nest: ', num2str(nest_ind)]);
    end
    
    %% evaluate the new alpha
    clear alpha_curr;
    alpha_curr = alpha/(nest_ind+1);
    alpha_CURRENT(nest_ind + 1, 1) = alpha_curr;

    %% noisy sampling
    currDist = (1 - alpha_curr)*PostDist + alpha_curr*sampling_noise_Dist;
    n_rlz = floor(cut_off/acc_rate);
    clear PostDist;
    clear PARAMETER err_vec;
    tic;
    PARAMETER = fn_prm_rlz(currDist, lb, ub, n_rlz);
    sampling_toc = toc;
    tic;
    err_vec = fn_calc_error_UD(PARAMETER, info);
    fn_eval_toc = toc;

    if nest_ind > 1
        if use_prev == 1
            clear err_prev prm_prev;
            err_prev = sort_ERR;
            prm_prev = sort_PRM;
        end
    end
    
    clear Sort_err index sort_P sort_PRM sort_ERR;
    tic;
    [Sort_err, index] = sort(err_vec);
    sort_P = PARAMETER(index, :);
    sort_PRM = sort_P(1:cut_off, :);
    sort_ERR = Sort_err(1:cut_off, 1);
    
    if nest_ind > 1
        if use_prev == 1
            q_err = [sort_ERR; err_prev];
            q_prm = [sort_PRM; prm_prev];
            clear Sort_err index sort_P sort_PRM sort_ERR;
            [Sort_err, index] = sort(q_err);
            sort_P = q_prm(index, :);
            sort_PRM = sort_P(1:cut_off, :);
            sort_ERR = Sort_err(1:cut_off, 1);
        end
    end
    sort_update_toc = toc;

    ERROR_CUTS(nest_ind, :) = [Sort_err(1, 1), Sort_err(cut_off, 1)];
    PRM_nest(:, :, nest_ind) = sort_PRM;
    ERROR_MAT(:, nest_ind) = Sort_err(1:cut_off, 1);
    [Pmf_all(:, :, nest_ind + 1)] = fn_posterior(sort_PRM, lb, ub, n_bins);
    clear Post_pmf;
    Post_pmf = Pmf_all(:, :, nest_ind + 1);
    PostDist = zeros(n_bins + 1, n_prm);
    for ind_bin = 2:(n_bins + 1)
        PostDist(ind_bin, :) = PostDist(ind_bin - 1, :) + Post_pmf(ind_bin - 1, :);
    end
    
    
    %%
    clear cmf_est_temp;
    save([info.info_fit_folder, '/', info.ABC_FAR_file, '.mat']);
    cmf_est_temp = fn_estimate_cmf(info, nest_ind);
    cmf_est = [cmf_est; cmf_est_temp];
    save([info.info_fit_folder, '/', info.ABC_FAR_file, '.mat'],...
        'cmf_est_temp', 'cmf_est', '-append');
    
    time_taken_iteration(nest_ind, :) = [sampling_toc, fn_eval_toc, sort_update_toc];
    if info.work_progress_time == 1
        disp(['This iteration (#', nest_ind, ') was completed in '...
            time_taken, 's.'])
    end
    clear sampling_toc fn_eval_toc sort_update_toc;
end

info.time_taken_iteration = time_taken_iteration;

ABC_end_time = datetime;
info.ABC_end_time = ABC_end_time;

total_time_taken = ABC_end_time - ABC_start_time;
info.total_time_taken = total_time_taken;

save([info.info_fit_folder, '/', info.ABC_FAR_file, '.mat']);
save([info.info_fit_folder, '/', info.ABC_FAR_results_file, '.mat'],...
    'info', 'total_time_taken', 'time_taken_iteration',...
    'Pmf_all', 'PRM_nest', 'cmf_est_temp', 'cmf_est',...
    'ERROR_CUTS', 'ERROR_MAT');


if diary_reqd == 1
    disp('-----');
    disp(['Diary is saved: ', info.info_fit_folder, '\', diary_file, '.txt']);
    diary([info.info_fit_folder, '\', diary_file, '.txt']);
    disp('-----');
    disp('Information on fitting analysis')
    disp(info); disp('-----');
    disp('Range of initial guess: ');
    disp('Lower bound'); disp(info.lb');
    disp('Upper bound'); disp(info.ub');
    disp('-----');
    if isfield(info, 'default')
        disp('Default values for the following were used:');
        fields_default = fieldnames(info.default);
        disp(fields_default);
    end
    disp('-----');
    disp('ABC fitting done using:');
    disp(['Release R' version('-release')]);

    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        poolsize = 0;
    else
        poolsize = p.NumWorkers;
    end
    disp('Parallel pool information');
    disp(['Pool connected: ', num2str(p.Connected)]);
    disp(['Number of workers: ', num2str(poolsize)]);
    disp('-----');
    diary off;
end

end
