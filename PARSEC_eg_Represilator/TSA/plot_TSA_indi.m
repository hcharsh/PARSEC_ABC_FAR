% clear all; clc; close all; warning off; format shortG; % % housekeeping

function plot_TSA_indi(tim_d, prm_ind)

folder_name = 'Results';
TSA_file = [folder_name, '/TSA_prm',  num2str(prm_ind), '_',  num2str(tim_d), 'unit.mat'];
file_name = [folder_name, '/plot_TSA_prm',  num2str(prm_ind), '_',  num2str(tim_d), 'unit'];

%% load chosen parameter sets
load(TSA_file, 'Si_mat', 'Sti_mat', 'prm_name', 'MT', 'var_interest'); %, 'tim_d', 'system_ind');
n_tim = size(MT, 1); mat_Tm = Sti_mat';
n_Si  = size(mat_Tm, 1); n_prm = size(mat_Tm, 2); n_QI = size(var_interest, 1);

%% figure formatting
fx0 = 0.05; fy0 = 0.1; fx1 = 1; fy1 = n_QI*0.9; lw = 1.5;
fs = 24;
clear f ax s;
f = figure('Units', 'normalized', 'Position',[fx0 fy0 fx1 fy1]);
% f = figure;
n_rows = floor((n_QI - 1)/2) + 1;
for ind_QI = 1:n_QI
    clear ax;
    ax = subplot(n_rows,2,ind_QI);
    % disp(var_interest(ind_QI, 1));
    clear Z;
    st = (ind_QI - 1) * n_tim + 1; fn = ind_QI * n_tim;
    Z = [mat_Tm(st:fn, 1:n_prm - 1), mat_Tm(st:fn, n_prm - 1)];
    x1 = size(Z, 1); x2 = size(Z, 2);
    % close all;
    s = surf(tim_d*[1:x1], [1:x2], Z');
    xlim([tim_d - 1, tim_d*x1 + 1])
    ylim([1, x2])
    title(['V ', num2str(var_interest(ind_QI, 1))])
    s.EdgeColor = 'none';
    colorbar;
    set(ax, 'TickLength', [0.0125 0.0125], 'TickDir', 'out',...
        'YTick', [1:1:x2 + 0*n_prm]+0.5, 'YTickLabel', prm_name...
        , 'Linewidth', lw, 'Fontsize', fs...
        , 'XTick', tim_d*[1, x1/4, x1/2, 3*x1/4, x1]..., 'XTickLabel', {'3', '18', '36', '54', '72'}...
        );
    caxis(ax,[0, 1]);
    view(ax,[0 90]);
    grid off;
end
set(ax, 'TickLength', [0.0125 0.0125], 'TickDir', 'out',...
    'YTick', [1:1:x2 + 0*n_prm]+0.5, 'YTickLabel', prm_name...
    , 'Linewidth', lw, 'Fontsize', fs...
    , 'XTick', tim_d*[1, x1/4, x1/2, 3*x1/4, x1]..., 'XTickLabel', {'3', '18', '36', '54', '72'}...
    );
xlabel('Time (unit)')
save([file_name, '.mat']);
end