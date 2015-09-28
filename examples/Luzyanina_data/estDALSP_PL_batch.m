% This file evaluationes the uncertainties of the parameter of a cell 
% population model describing the second data set in Luzyanina et. al, 2007.
%
% Uncertainty analysis method: Profile likelihood
%
%
% Example for input:
%   comp_mode = 'sequential';
%   M_number = 1;
%   PL_index = 1;

function [parameters,M_number] = estDALSP_PL_batch(comp_mode,M_number,PL_index)

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

%% PROFILE CALCULATION -- PARAMETERS
% Options
options_PL_par = options;
options_PL_par.R_min = exp(-chi2inv(0.99,1)/2-1); % single parameter CI
% options_PL.R_min = exp(-chi2inv(0.99,parameters.number)/2-1); % simultaneous CI
options_PL_par.dR_max = 0.15;
options_PL_par.options_getNextPoint.min = 1e-4;
options_PL_par.options_getNextPoint.max = 1e0;
options_PL_par.parameter_index = PL_index;
options_PL_par.fmincon = optimset('algorithm','active-set',...
                                  'display','iter-detailed',...
                                  'GradObj','on',...
                                  'MaxIter',1000,...
                                  'MaxFunEvals',1000*parameters.number);

% Calculation
parameters = getParameterProfiles(parameters,@(theta) logLikelihood_proliferation(theta,M,D,options_logL),options_PL_par);

% %% PROFILE CALCULATION -- PROPERTIES
% % Options
% options_PL_prop = options;
% options_PL_prop.R_min = exp(-chi2inv(0.99,1)/2-1); % single parameter CI
% options_PL_prop.dR_max = 0.05;
% options_PL_prop.dJ = 0.5;
% options_PL_prop.fmincon = optimset('algorithm','active-set',...
%                                    'display','iter-detailed',...
%                                    'GradObj','on',...
%                                    'GradConstr','on',...
%                                    'MaxIter',100,...
%                                    'TolCon',1e-4,...
%                                    'MaxSQPIter',100*parameters.number,...
%                                    'MaxFunEvals',1000*parameters.number);
% options_PL_prop.property_index = PL_index;
% 
% % Definition of property struct
% if exist('properties') ~= 1
%     for i = 1:parameters.number
%         properties.name{i} = parameters.name{i};
%         properties.function{i} = @(theta) propertyFunction_theta(theta,i);
%         properties.min(i) = parameters.min(i);
%         properties.max(i) = parameters.max(i);
%     end
%     properties.number = parameters.number;
% end
% 
% % Evaluate properties for multi-start results
% properties = getPropertyMultiStarts(properties,parameters,options);
% 
% % Calculation
% properties = getPropertyProfiles(properties,parameters,@(theta) logLikelihood_proliferation(theta,M,D,options_logL),options_PL_prop);

% SAVE RESULTS
save(['./project/results/' M.name '__grad_' options.fmincon.GradObj ...
            '__PL_par__' strrep(num2str(PL_index),' ','_')]);

end