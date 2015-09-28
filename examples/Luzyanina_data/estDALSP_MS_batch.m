% This file estimates the parameter of a cell population model
% describing the second data set in Luzyanina et. al, 2007.
%
% Estimation Methods: Multi-start local optimization
%
%
%
% Example for input:
%   comp_mode = 'sequential';
%   M_number = 1;

function [parameters,M_number] = estDALSP_MS_batch(comp_mode,M_number)

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

% Optimization
options.n_starts = 250;
options.comp_type = comp_mode;
options.mode = 'text';
options.fmincon = optimset('algorithm','interior-point',...
                           'display','iter-detailed',...
                           'GradObj','on',...
                           'MaxIter',1000,...
                           'TolFun',1e-9,...
                           'TolX',1e-9,...
                           'MaxFunEvals',1000*parameters.number);

%% OPTIMIZATION
parameters = getMultiStarts(parameters,@(theta) logLikelihood_proliferation(theta,M,D,options_logL),options);

% Model selection criteria
parameters.AIC = -2*parameters.MS.logPost + 2*parameters.number;
parameters.BIC = -2*parameters.MS.logPost + parameters.number*log(sum(sum(D.cellcount)));

%% SAVE RESULTS
save(['./project/results/' M.name '__grad_' options.fmincon.GradObj]);

end