function fit_setup_file_UD(expt_folder, res_folder, pf, ed_index, lb, ub)

% addpath '../functions_ABC_FAR';
addpath '../../functions_ABC_FAR';

%% options for the simulations
n_nest = 10;
acc_rate = 0.01;
n_bins = 25;
cut_off = 16*n_bins;
alpha = 0.25;
use_prev = 1;

%% data, id and folders
ED_file = [expt_folder, '/', pf, num2str(ed_index), '.mat'];

%% setting the options
info.far = acc_rate; % fixed acceptance rate (can be made into a vector)
% number of samples being used to get the updated distribution
info.cut_off = cut_off;
% include the best parameter vectors selected in prev. step to update
% parameter distribution in the current step. Default is 1.
info.use_prev = use_prev;
% number of nests or iterations
info.n_nest = n_nest;
% max. sampling error without dummy unbiasing
info.alpha = alpha;
% number of itervals to discretize the each parameter range
info.n_bins = n_bins;

%% parameter ranges, prior guess and sampling distribution
% lower (lb) and upper (ub) bounds of each parameter.
% lb and ub are vectors (nx1) following the indexing of info_full.prm_name
info.lb = lb;
info.ub = ub;

% unif_Dist = matrix with n_prm columns of uniform distributions
n_prm = size(lb, 1);
for ind_prm = 1:n_prm
    unif_Dist(:, ind_prm) = linspace(0, 1, n_bins + 1)'; % uniform distribution
end

% the initial guess (or prior distribution) and noise distribution to be
% added while sampling with range normalized between 0 and 1.
info.PriorDistName = 'Uniform prior';
info.PriorDist = unif_Dist;
% To input other distributions, input as a matrix with each columns
% specifying the distribution for each parameter.

info.sampling_noise_DistName = 'Uniform Noise';
info.sampling_noise_Dist = unif_Dist;

%% additional options
info.diary_reqd = 0;
info.work_progress_nest = 0;
info.work_progress_time = 0;

if exist(ED_file)
    load(ED_file);
    info.error_calc_fn = 'error_function(prm, data_file, prm_sample_ind)'; % to be changed in 'fn_calc_error_UD'
    info.data_file = ED_file;
    info.mainfolder = res_folder;

    % folder name in which the information and fit result files are to be stored
    info.info_fit_folder = res_folder; % [res_folder, '/info_fit'];
    mkdir(info.info_fit_folder);
    
    for prm_sample_ind = 1:size(Data_set, 3)
        clear info1; info1 = info;
        info1.prm_sample_ind = prm_sample_ind;
        info1.system_id = [pf, num2str(ed_index), '_s', num2str(prm_sample_ind)];
        
        % file name for the fit analysis
        info1.ABC_FAR_file           = ['full_', info1.system_id];
        info1.ABC_FAR_results_file   = ['sh_', info1.system_id];
        % don't put .mat extn(s).
        
        fn_ABC(info1);

        %% done
        % disp(['Results saved in ', info.info_fit_folder, '/', info.ABC_FAR_file, '.mat']);
        disp(['Results saved in ', info.mainfolder, '/', info1.ABC_FAR_file, '.mat']);

    end

else
    disp('Expt. design file does not exist!');
end

end