% This file estimates the parameter of a cell population model
% describing the second data set in Luzyanina et. al, 2007.

clear all;
close all;
clc;

%% PREPROCESSING OF DATA
% Load data
load data_set_2_Luzyanina;
D = data;
% Process data;
D.t = D.t([2,3,4,5,6]);
D.t = D.t - D.t(1);
D.cellcount = D.cellcount([2,3,4,5,6],:);
D.t_name = {'1st day','2nd day','3rd day','4th day','5th day'};
D.t_plot = [1,2,3,4,5];

%% DEFINITION OF MODEL
% model_1__alpha__beta__a0_delta;
% model_2__alpha_a__beta__a0_delta;
model_2_4__alpha_01a__beta__a0_delta;
% model_3__alpha_i__beta__a0_delta;
% model_4__alpha_ia__beta__a0_delta;
% model_5__alpha__beta_a__a0_delta
% model_6__alpha_a__beta_a__a0_delta
% model_7__alpha_i__beta_a__a0_delta
% model_8__alpha_ia__beta_a__a0_delta
% model_9__alpha__beta_i__a0_delta
% model_10__alpha_a__beta_i__a0_delta
% model_11__alpha_i__beta_i__a0_delta
% model_12__alpha_ia__beta_i__a0_delta
% model_13__alpha__beta_ia__a0_delta
% model_14__alpha_a__beta_ia__a0_delta
% model_15__alpha_i__beta_ia__a0_delta
% model_16__alpha_ia__beta_ia__a0_delta

%% SIMULATION AND OPTIMIZATION PARAMETER
% Simulation
a_sim = linspace(D.t(1),D.t(end),800); % age vector for simulation
t_sim = linspace(D.t(1),D.t(end),800); % time vector for simulation

% Likelihood function
options_logL.simulation.a_sim = a_sim;
options_logL.simulation.t_sim = t_sim;
options_logL.simulation.noise.flag = 'yes';
options_logL.costfuntype = 'logL';

% Optimization
options.n_starts = 10;
options.comp_type = 'sequential'; options.mode = 'visual';
% options.comp_type = 'sequential'; options.mode = 'text';
% options.comp_type = 'parallel'; options.mode = 'text'; n_workers = 1;
options.fmincon = optimset('algorithm','interior-point',...
                           'display','iter-detailed',...
                           'GradObj','on',...
                           'MaxIter',1000,...
                           'TolFun',1e-9,...
                           'TolX',1e-9,...
                           'MaxFunEvals',1000*parameters.number);

% Prepare visualization
if strcmp(options.mode,'visual')
    % Figure handle
    options.fh = figure;
    options_logL.fh = figure;
    % Timer for visulaization of fit
    global tp
    tp = clock;
end

% Open Matlabpool
if strcmp(options.comp_type,'parallel') && exist('n_workers')
    if n_workers >= 2
        matlabpool(n_workers);
    end
end

%% OPTIMIZATION
[parameters,fh_opt] = getMultiStarts(parameters,...
    @(theta) logLikelihood_proliferation(theta,M,D,options_logL),options);

%% PLOT BEST FIT
% Simulation
[Sim] = CPsimulateDALSP(M,parameters.MS.par(:,1),D.t,a_sim,...
                        [D(1).bins(1,1); D(1).bins(:,2)]',...
                        options_logL.simulation);
                    
% Plot
if strcmp(options.mode,'visual')
    plotProliferationAssay(Sim,D,options_logL.fh);
end

% Model selection criteria
parameters.AIC = -2*parameters.MS.logPost + 2*parameters.number;
parameters.BIC = -2*parameters.MS.logPost + parameters.number*log(sum(sum(D.cellcount)));

%% SAVE RESULTS
pathname = ['./project/results/' M.name '__grad_' options.fmincon.GradObj];
% Data
save(pathname);
% Figures
if strcmp(options.mode,'visual')
    saveas(options_logL.fh,[pathname '_fit'],'fig');
    saveas(options.fh  ,[pathname '_MS'] ,'fig');
end

%% VISUALIZE MEASUREMENT UNCERTAINTIES
if strcmp(options.mode,'visual')
plotProliferationAssay_w_MU(Sim,D,options_logL.fh);
end

% Close Matlabpool
if strcmp(options.comp_type,'parallel') && exist('n_workers')
    if n_workers >= 2
        matlabpool('close');
    end
end
