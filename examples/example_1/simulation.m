clear all;
close all;
clc;

%% DEFINITION OF MODEL
model__w__time_dependence;

%% SIMULATION
% Simulation parameters
x_sim = logspace(1,4,2^8+1); % concentration vector for simulation
a_sim = linspace(0,4,1000);  % age vector for simulation
t_sim = linspace(0,4,1000);  % time vector for simulation
t = linspace(0,4,7); % time vector for observation

options_simulation.a_sim = a_sim;
options_simulation.t_sim = t_sim;
options_simulation.noise.flag = 'yes';

% Simulation
[Sim] = CPsimulateDALSP(M,parameters.guess,t,a_sim,x_sim,options_simulation);

%% VISUALIZATION
options_plot.plot_simulation = 'true';
options_plot.plot_data = 'false';
plotProliferationAssay(Sim,[],[],options_plot);

%%
D = data;
% Process data;
D.t = D.t([2,3,4,5,6]);
D.t = D.t - D.t(1);
D.cellcount = D.cellcount([2,3,4,5,6],:);
D.t_name = {'1st day','2nd day','3rd day','4th day','5th day'};
D.t_plot = [1,2,3,4,5];

