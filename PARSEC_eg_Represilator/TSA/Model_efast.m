% % First order and total effect indices for a given
% % model computed with Extended Fourier Amplitude
% % Sensitivity Test (EFAST).
% % Andrea Saltelli, Stefano Tarantola and Karen Chan.
% % 1999. % "A quantitative model-independent method for global
% % sensitivity analysis of model output". % Technometrics 41:39-56
% clear all;
% close all;
% warning off;
% % clc;

function [Z, TSA_file] = Model_efast(tim_d, var_interest, prm_ind)

load('GT_val.mat', 'prm', 'tspan', 'y0', 'prm_name', 'var_label');
dummy = 1;
prm_vec = [prm(:, prm_ind); dummy]; clear prm;
TSA_file = ['Results/TSA_prm',  num2str(prm_ind), '_',  num2str(tim_d), 'unit.mat'];
if ~exist('Results', 'dir')
       mkdir('Results')
end
%% INPUT
NR = 5; %: no. of search curves - RESAMPLING

k = size(prm_vec, 1); % # of input factors (parameters varied) + dummy parameter
% # of inout factors same as number of free parameters

NS = 500; % # of samples per search curve

wantedN=NS*k*NR; % wanted no. of sample points

% OUTPUT
% SI[] : first order sensitivity indices
% STI[] : total effect sensitivity indices
% Other used variables/constants:
% OM[] : vector of k frequencies
% OMi : frequency for the group of interest
% OMCI[] : set of freq. used for the compl. group
% X[] : parameter combination rank matrix
% AC[],BC[]: fourier coefficients
% FI[] : random phase shift
% V : total output variance (for each curve)
% VI : partial var. of par. i (for each curve)
% VCI : part. var. of the compl. set of par...
% AV : total variance in the time domain
% AVI : partial variance of par. i
% AVCI : part. var. of the compl. set of par.
% Y[] : model output

MI = 4; %: maximum number of fourier coefficients
% that may be retained in calculating the partial
% variances without interferences between the
% assigned frequencies

%% PARAMETERS AND ODE SETTINGS (they are included in the following file)
% Parameter_settings_EFAST

pmin = prm_vec' - log10(1.2);
pmax = prm_vec' + log10(1.2);
% disp(prm_vec);
clear prm_vec;
% efast_var = prm_name;
% y_var_label = var_label;% Variables Labels

% disp([prm_vec, pmin, pmax])
% TIME SPAN OF THE SIMULATION
t_bgn = tspan(1, 1); t_end = tspan(1, 2);
clear tspan; % disp(tim_d);
tspan = [t_bgn : tim_d : t_end]; % time points of model calculation
MT = [tim_d : tim_d : t_end]'; % time points of QT calculation

%% Computating the frequency
% Computation of the frequency for the group
% of interest OMi and the # of sample points NS (here N=NS)
OMi = floor(((wantedN/NR)-1)/(2*MI)/k);
NS = 2*MI*OMi+1;
if(NS*NR < 65)
    fprintf(['Error: sample size must be >= ' ...
    '65 per factor.\n']);
    return;
end


%% Pre-allocation of the output matrix Y
%% Y will save only the points of interest specified in
%% the vector time_points
% Y(NS,length(1),length(y0),length(pmin),NR)=0;  % pre-allocation

% % % % % % % % % % % % % % % % % % % % % % % tic;
% Loop over k parameters (input factors)
for i=1:k           % i=# of replications (or blocks)
    % Algorithm for selecting the set of frequencies.
    % OMci(i), i=1:k-1, contains the set of frequencies
    % to be used by the complementary group.
    OMci = SETFREQ(k,OMi/2/MI,i);   
    % Loop over the NR search curves.
    for L=1:NR
%         disp([i L]);
        % Setting the vector of frequencies OM
        % for the k parameters
        cj = 1;
        for j=1:k
            if(j==i)
                % For the parameter (factor) of interest
                OM(i) = OMi;
            else
                % For the complementary group.
                OM(j) = OMci(cj);
                cj = cj+1;
            end
        end
        % Setting the relation between the scalar
        % variable S and the coordinates
        % {X(1),X(2),...X(k)} of each sample point.
        FI = rand(1,k)*2*pi; % random phase shift
        S_VEC = pi*(2*(1:NS)-NS-1)/NS;
        OM_VEC = OM(1:k);
        FI_MAT = FI(ones(NS,1),1:k)';
        ANGLE = OM_VEC'*S_VEC+FI_MAT;
        
        X(:,:,i,L) = 0.5+asin(sin(ANGLE'))/pi; % between 0 and 1
        
        % Transform distributions from standard
        % uniform to general.
        X(:,:,i,L) = parameterdist(X(:,:,i,L),pmax,pmin,0,1,NS,'unif'); %%this is what assigns 'our' values rather than 0:1 dist
        % Do the NS model evaluations.
        for run_num=1:NS
            % [i run_num L] % keeps track of [parameter run NR]
            % ODE solver call
            clear prm_vec t y;
            prm_vec(:, 1) = X(run_num, :, i, L)';            
            [t,y]=ode23s(@(t, y) model_system(t, y, prm_vec), tspan, y0);
            
            MY = interp1(t, y, MT);
            QT_interest = analyze_this(MY, var_interest); % disp(QT_interest(1:4, 1)');
            for QT_ind = 1:size(QT_interest, 1)
                Y(run_num, 1, QT_ind, i, L) = QT_interest(QT_ind,1);
            end
            
        end %run_num=1:NS
    end % L=1:NR
end % i=1:k
% time_taken = toc;
% disp(['The analysis took ', num2str(time_taken), 's.'])
disp(TSA_file);

save(TSA_file);

sy3= size(Y, 3);
% CALCULATE Si AND STi for each resample (1,2,...,NR) [ranges]
[Si,Sti,rangeSi,rangeSti] = efast_sd(Y,OMi,MI,1,1:sy3);
save(TSA_file);

% clear all;
% load(['TSA_',  num2str(tim_d), 'hr.mat'], 'Si', 'Sti');

for a = 1:size(Si, 1)
    for b = 1:size(Si, 3)
        Si_mat(a, b) = Si(a, 1, b)/sum(Si(:, 1, b));
        Sti_mat(a, b) = Sti(a, 1, b)/sum(Sti(:, 1, b));
    end
end

save(TSA_file);


%%
n_tim = size(MT, 1); mat_Tm = Sti_mat';
n_Si  = size(mat_Tm, 1); n_prm = size(mat_Tm, 2); n_QI = 2; size(var_interest, 1);
Z = mat_Tm(1 : n_tim, 1:n_prm - 1);

for ind_QI = 2:n_QI
    st = (ind_QI - 1) * n_tim + 1; fn = ind_QI * n_tim;
    Z = [Z, mat_Tm(st:fn, 1:n_prm - 1)];
end
Z = [Z, mat_Tm(st:fn, n_prm - 1)];
clear n_tim mat_Tm n_Si n_prm n_QI ind_QI st fn;
save(TSA_file);

end