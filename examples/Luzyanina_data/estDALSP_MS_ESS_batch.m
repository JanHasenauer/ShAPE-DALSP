% This file estimates the parameter of a cell population model
% describing the second data set in Luzyanina et. al, 2007.
%
% Estimation Methods: ESS
%
%
%
% Example for input:
%   M_number = 1;

function [parameters,M_number] = estDALSP_MS_ESS_batch(M_number)

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
options_logL_ESS.simulation.a_sim = a_sim;
options_logL_ESS.simulation.t_sim = t_sim;
options_logL_ESS.simulation.noise.flag = 'yes';
options_logL_ESS.costfuntype = 'logL';
options_logL_ESS.sign = 'negative';

% Optimization
options.n_starts = 2;

%% PARTICLE SWARM OPTIMIZATION
% Problem definition
Problem_ESS.f = 'logLikelihood_proliferation';
Problem_ESS.x_L = parameters.min;
Problem_ESS.x_U = parameters.max;
Problem_ESS.c_L = [];
Problem_ESS.c_U = [];

% Algorithm options
Options_ESS.maxeval = 3e2; %1000*parameters.number % default setting
Options_ESS.maxtime = 100;
Options_ESS.iterprint = 1; % 0

% Initialization
par = nan(parameters.number,options.n_starts);
logPost = nan(options.n_starts,1);
n_objfun = nan(options.n_starts,1);
n_iter = nan(options.n_starts,1);
t_cpu = nan(options.n_starts,1);
exitflag = nan(options.n_starts,1);

global counter
counter = 0;

% Loop: Mutli-starts
disp('start loop');
parfor i = 1:options.n_starts
    % Run the algorithm
    tcpu = cputime;
    Results = MEIGO(Problem_ESS,Options_ESS,'ESS',M,D,options_logL_ESS);
    tcpu = cputime - tcpu;

    % Assignment
    par(:,i) = Results.xbest';
	logPost(i) = -Results.fbest;
	n_objfun(i) = Results.numeval;
	n_iter(i) = length(Results.f);
    exitflag(i) = Results.end_crit;
	t_cpu(i) = tcpu;
end

% Assignment
parameters.MS.par = par;
parameters.MS.logPost = logPost;
parameters.MS.n_objfun = n_objfun;
parameters.MS.n_iter = n_iter;
parameters.MS.t_cpu = t_cpu;
parameters = sortMultiStarts(parameters);

%% SAVE RESULTS
save(['./project/results/' M.name '__ESS' ...
       '__MaxEval_' num2str(Options_ESS.maxeval,'%d')],...
       'parameters','Problem_ESS','Options_ESS','M_number');

end