% This file estimates the parameter of a cell population model
% describing the second data set in Luzyanina et. al, 2007.
%
% Estimation Methods: Particle Swarm Optimitaion
%
%
%
% Example for input:
%   M_number = 1;
%   setting = 3;

function [parameters,M_number] = estDALSP_MS_PS_batch(M_number,setting)

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
options_logL_PS.simulation.a_sim = a_sim;
options_logL_PS.simulation.t_sim = t_sim;
options_logL_PS.simulation.noise.flag = 'yes';
options_logL_PS.costfuntype = 'logL';
options_logL_PS.sign = 'negative';

% Optimization
options.n_starts = 250;

%% PARTICLE SWARM OPTIMIZATION
% Problem definition
Problem_PS.ObjFunction = @(theta) logLikelihood_proliferation(theta,M,D,options_logL_PS);
Problem_PS.LB = parameters.min;
Problem_PS.UB = parameters.max;

% Algorithm options
switch setting
    case 1
        Options_PS.Size    = 20;
        Options_PS.MaxIter = 1e3*parameters.number;
        Options_PS.MaxObj  = 1e3*parameters.number;
    case 2
        Options_PS.Size    = 20;
        Options_PS.MaxIter = 1e4*parameters.number;
        Options_PS.MaxObj  = 1e4*parameters.number;
    case 3
        Options_PS.Size    = 1e1*parameters.number;
        Options_PS.MaxIter = 1e3*parameters.number;
        Options_PS.MaxObj  = 1e3*parameters.number;
    case 4
        Options_PS.Size    = 1e1*parameters.number;
        Options_PS.MaxIter = 1e4*parameters.number;
        Options_PS.MaxObj  = 1e4*parameters.number;
end
Options_PS.IPrint = -1;

% Initialization
par = nan(parameters.number,options.n_starts);
logPost = nan(options.n_starts,1);
n_objfun = nan(options.n_starts,1);
n_iter = nan(options.n_starts,1);
t_cpu = nan(options.n_starts,1);

% Initial guess
x0 = bsxfun(@plus,parameters.min,bsxfun(@times,rand(parameters.number,options.n_starts),(parameters.max-parameters.min)));

% Loop: Mutli-starts
parfor i = 1:options.n_starts
    % Run the algorithm
    tcpu = cputime;
    [theta,J,RunData] = PSwarm(Problem_PS,struct('x',x0(:,i)),Options_PS);
    tcpu = cputime - tcpu;

    % Assignment
    par(:,i) = theta;
	logPost(i) = -J;
	n_objfun(i) = RunData.ObjFunCounter;
	n_iter(i) = RunData.IterCounter;
	t_cpu(i) = tcpu;
    
    % Save results in between
    path = [ './project/results/'  num2str(M_number,'%d') '__PS__Size_' num2str(Options_PS.Size,'%d') '__MaxObj_' num2str(Options_PS.MaxObj,'%d') '__Start_' num2str(i)];
    dlmwrite([path '_theta'],theta,'delimiter',',','precision',12);
    dlmwrite([path '_J'],J,'delimiter',',','precision',12);
    dlmwrite([path '_ObjFunCounter'],RunData.ObjFunCounter,'delimiter',',','precision',12);
    dlmwrite([path '_IterCounter'],RunData.IterCounter,'delimiter',',','precision',12);
    dlmwrite([path '_tcpu'],tcpu,'delimiter',',','precision',12);

end

% Assignment
parameters.MS.par = par;
parameters.MS.logPost = logPost;
parameters.MS.n_objfun = n_objfun;
parameters.MS.n_iter = n_iter;
parameters.MS.t_cpu = t_cpu;
parameters = sortMultiStarts(parameters);

%% SAVE RESULTS
save(['./project/results/' M.name '__PS' ...
       '__Size_' num2str(Options_PS.Size,'%d') ...
       '__MaxObj_' num2str(Options_PS.MaxObj,'%d')],...
       'parameters','Problem_PS','Options_PS','options_logL_PS','M_number');

end