% This file estimates the parameter of a cell population model
% describing the second data set in Luzyanina et. al, 2007.
%
% Uncertainty analysis method: Markov chain Monte Carlo (MCMC)
%
%
%
% Example for input:
%   M_number = 16;
%   setting = 2;

function [parameters,M_number] = estDALSP_MCMC_batch(M_number,setting)

%% PREPROCESSING OF DATA
% Load data
load('./project/data/data_set_2_Luzyanina/data_set_2_Luzyanina.mat')
D = data;

% Process data;
D.t = D.t([2,3,4,5,6]);
D.t = D.t - D.t(1);
D.cellcount = D.cellcount([2,3,4,5,6],:);
D.t_name = {'1st day','2nd day','3rd day','4th day','5th day'};
D.t_plot = [1,2,3,4,5];

%% DEFINITION OF MODEL
switch M_number
    case 1
        model_1__alpha__beta__a0_delta;
    case 2
        model_2__alpha_a__beta__a0_delta;
    case 3
        model_3__alpha_i__beta__a0_delta;
    case 4
        model_4__alpha_ia__beta__a0_delta;
    case 5
        model_5__alpha__beta_a__a0_delta
    case 6
        model_6__alpha_a__beta_a__a0_delta
    case 7
        model_7__alpha_i__beta_a__a0_delta
    case 8
        model_8__alpha_ia__beta_a__a0_delta
    case 9
        model_9__alpha__beta_i__a0_delta
    case 10
        model_10__alpha_a__beta_i__a0_delta
    case 11
        model_11__alpha_i__beta_i__a0_delta
    case 12
        model_12__alpha_ia__beta_i__a0_delta
    case 13
        model_13__alpha__beta_ia__a0_delta
    case 14
        model_14__alpha_a__beta_ia__a0_delta
    case 15
        model_15__alpha_i__beta_ia__a0_delta
    case 16
        model_16__alpha_ia__beta_ia__a0_delta
end

%% SIMULATION AND OPTIMIZATION PARAMETER
% Simulation
a_sim = linspace(D.t(1),D.t(end),800); % age vector for simulation
t_sim = linspace(D.t(1),D.t(end),800); % time vector for simulation

% Likelihood function
options_logL.simulation.a_sim = a_sim;
options_logL.simulation.t_sim = t_sim;
options_logL.simulation.noise.flag = 'yes';
options_logL.costfuntype = 'logL';

%% LOAD RESULTS OF MULTI-START OPTIMIZATION
load(['./project/results/' M.name '__grad_on']);
options.mode = 'visual';
options.comp_type = comp_mode;

%% Single-chain Markov chain Monte-Carlo sampling -- Parameters

% Algorithm options
options.mode = 'silent';
switch setting
    case 1
        options.sampling_scheme = 'DRAM';
        options.AM.memory_length = inf;
        options.nsimu_warmup = 5e4;
        options.nsimu_run    = 5e5;
    case 2
        options.sampling_scheme = 'single-chain';
        options.proposal_scheme = 'AM';
        options.AM.adaption_scheme = 'difference';
        options.AM.memory_length = inf;
        options.thinning = 1;
        options.nsimu_warmup = 1e3;
        options.nsimu_run = 2e5;
end

% Sampling
parameters = getParameterSamples(parameters,@(theta) logLikelihood_proliferation(theta,M,D,options_logL),options);

% %%
% clc;
% options.theta_0 = parameters.S.par(:,end);
% options.Sigma_0 = cov(parameters.S.par');
% options.nsimu_warmup = 0;
% options.nsimu_run = 1e4;
%         
% parameters_old = parameters;
% parameters_new = getParameterSamples(parameters,@(theta) logLikelihood_proliferation(theta,M,D,options_logL),options);
% 
% %%
% parameters.S.par = [parameters_old.S.par,parameters_new.S.par];
% parameters.S.logPost = [parameters_old.S.logPost;parameters_new.S.logPost];
% 
% %%
% parameters_short = parameters;
% parameters_short.S.par = parameters_short.S.par(:,30000:end);
% parameters_short.S.logPost = parameters_short.S.logPost(30000:end);
% 
% % Diagnosis plots
% plotMCMCdiagnosis(parameters_short,'log-posterior');
% plotMCMCdiagnosis(parameters_short,'parameters');
% 
% % Parameter distribution
% plotParameterUncertainty(parameters_short,'1D');
% % plotParameterSamples(parameters,'2D');
% 
% % Chain statistics
% chainstats(parameters_short.S.par');

%% SAVE RESULTS
save(['./project/results/' M.name '__grad_' options.fmincon.GradObj '__MCMC__SampleSize_' num2str(options.nsimu_run)]);

end