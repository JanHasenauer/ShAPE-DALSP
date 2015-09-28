clear all;
close all;
clc;

%% Load data
load('./project/results/model_16__alpha_ia__beta_ia__a0_delta__grad_on');

%% Plot
%options_plot.lincol = colormap(gray(6));
%plotProliferationAssay_w_MU(Sim,D,[],options_plot);


%% Computation
% Calculation of bin probability
x = [D(1).bins(1,1); D(1).bins(:,2)]';
Sim.H = bsxfun(@times,Sim.n_y(:,1:end-1) + Sim.n_y(:,2:end), 0.5*diff(x(:)'));
Sim.h = bsxfun(@times,Sim.p_y(:,1:end-1) + Sim.p_y(:,2:end), 0.5*diff(x(:)'));

% Calculation of lower and upper bounds
lb = zeros(size(Sim.H));
ub = zeros(size(Sim.H));
perc = [0.005 0.995];
for j = 1:length(D.t)
    for i = 1:length(lb)
        b = binoinv(perc,sum(D.cellcount(j,:)),Sim.h(j,i));
        lb(j,i) = b(1);
        ub(j,i) = b(2);
    end
end

% Calculation of population size on fine time scale
t_fine = linspace(0,4,100);
[Sim2] = CPsimulateDALSP(M,parameters.MS.MAP.par,t_fine,a_sim,x,options_logL.simulation);

%% Visualization
figure;
for j = 1:length(D.t)
    subplot(2,3,j);
    % Measurement
    lh(1)=stairs(x,(D.cellcount(j,[1:end,end])/sum(D.cellcount(j,:))),'-','Color',0.6*[1,1,1],'linewidth',1); hold on;
    % Simulation
    lh(2)=stairs(x,Sim.h(j,[1:end,end]),'-','Color',[0,0,0],'linewidth',1.5);
    lh(3)=stairs(x,lb(j,[1:end,end])/sum(D.cellcount(j,:)),'-','Color',[0,0,0],'linewidth',1);
    lh(3)=stairs(x,ub(j,[1:end,end])/sum(D.cellcount(j,:)),'-','Color',[0,0,0],'linewidth',1);

    % Set legende label
    set(gca,'Xscale','log');
    xlim([5,x(end)]);
    ylabel('frequency');
    set(gca,'xtick',[10,100,1000])
    xlabel('flourescence intensity, {\ity} [UI]');
    if j == 1
        legend(lh,{'measurement','simulation','confidence interval (99%)'},'location','NorthWest');
    end
    switch j
        case 1
            ylim([0,0.012]);
        case 2
            ylim([0,0.007]);
        case 3
            ylim([0,0.006]);
        case 4
            ylim([0,0.005]);
        case 5
            ylim([0,0.004]);
    end
            yl = get(gca,'ylim');
text(900,yl(2)*0.9,['{\bf day ' num2str(D.t_plot(j),'%d') '}']);
end

subplot(2,3,6);
plot(t_fine+D.t_plot(1),Sim2.N,'-','Color',[0,0,0],'linewidth',1.5); hold on;
plot(D.t_plot,sum(D.cellcount,2),'o','linewidth',1.5,'Color',0.6*[1,1,1]);
f_perc = icdf('norm',perc,0,options_logL.noise.sigma_mean);
plot(t_fine+D.t_plot(1),exp(log(Sim2.N)+f_perc(1)),'-','Color',[0,0,0],'linewidth',1);
plot(t_fine+D.t_plot(1),exp(log(Sim2.N)+f_perc(2)),'-','Color',[0,0,0],'linewidth',1);

set(gca,'xtick',1:5,'xlim',[1,5])
xlabel('time, {\itt} [days]');
ylabel('# cells, {\itN} [-]');